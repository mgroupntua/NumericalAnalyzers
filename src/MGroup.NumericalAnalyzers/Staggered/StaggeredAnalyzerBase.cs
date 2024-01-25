using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using MGroup.LinearAlgebra.Iterative;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers.Staggered
{
	public delegate void CreateNewModelDelegate(IParentAnalyzer[] analyzers, ISolver[] solvers); 
	
	public abstract class StaggeredAnalyzerBase : IAnalyzer
	{
		protected int maxStaggeredSteps;
		protected double tolerance;
		protected IterativeStatistics analysisStatistics = new IterativeStatistics() { AlgorithmName = "Base staggered analyzer" };
		protected List<IList<IList<IterativeStatistics>>> nestedAnalysisStatistics = new List<IList<IList<IterativeStatistics>>>();
		protected IParentAnalyzer[] analyzers;
		protected ISolver[] solvers;
		protected CreateNewModelDelegate CreateNewModel;
		protected IGlobalVector[] currentSolutions;
		protected GenericAnalyzerState[] analyzerStates;

		public IAnalysisWorkflowLog[] Logs { get; set; }

		public IterativeStatistics AnalysisStatistics => analysisStatistics;

		public IList<IList<IList<IterativeStatistics>>> NestedAnalysisStatistics => nestedAnalysisStatistics;

		public IGlobalVector CurrentAnalysisResult { get => throw new NotSupportedException("Staggered analyzer has more than one nested analyzers. Use CurrentAnalysisResult of each individual nested analyzer"); }

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => CreateState();
			set
			{
			}
		}

		GenericAnalyzerState CreateState() => new GenericAnalyzerState(this, new[]
		{
			(String.Empty, (IGlobalVector)null)
		});

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();

		public void BuildMatrices()
		{
			foreach (var analyzer in analyzers)
			{
				analyzer.BuildMatrices();
			}
		}

		/// <summary>
		/// Initializes the values of the system to be solved.
		/// </summary>
		public void Initialize(bool isFirstAnalysis = true)
		{
			foreach (var analyzer in analyzers)
			{
				analyzer.Initialize(isFirstAnalysis);
			}
		}

		protected virtual void Solve(Func<Action[]> solveMethods)
		{
			int staggeredStep = 0;
			var solutionNorm = 0d;
			double previousSolutionNorm, error;

			do
			{
				Debug.WriteLine("\n\nStaggered step: {0}", staggeredStep);
				previousSolutionNorm = solutionNorm;
				var currentStatistics = new IList<IterativeStatistics>[analyzers.Length];

				for (int i = 0; i < solvers.Length; i++)
				{
					if (solvers[i].LinearSystem.Solution != null)
					{
						currentSolutions[i] = solvers[i].LinearSystem.Solution.Copy();
					}

					solveMethods()[i]();

					currentStatistics[i] = analyzers[i].AnalysisStatistics.ToArray();
				}

				nestedAnalysisStatistics.Add(currentStatistics);
				solutionNorm = 0;
				for (int i = 0; i < solvers.Length; i++)
				{
					solutionNorm += solvers[i].LinearSystem.Solution.Norm2();
				}

				error = solutionNorm != 0 ? Math.Abs(solutionNorm - previousSolutionNorm) / solutionNorm : 0;
				Debug.WriteLine("Staggered step: {0} - error {1}", staggeredStep, error);
				staggeredStep++;

				if (staggeredStep < maxStaggeredSteps && error > tolerance)
				{
					CreateNewModel(analyzers, solvers);
				}
			}
			while (staggeredStep < maxStaggeredSteps && error > tolerance);

			analysisStatistics.NumIterationsRequired = staggeredStep;
			analysisStatistics.ResidualNormRatioEstimation = error;
			analysisStatistics.HasConverged = staggeredStep < maxStaggeredSteps;
		}

		public abstract void Solve();
	}
}
