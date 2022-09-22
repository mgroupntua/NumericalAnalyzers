using System;
using System.Diagnostics;

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
		protected IParentAnalyzer[] analyzers;
		protected ISolver[] solvers;
		protected CreateNewModelDelegate CreateNewModel;
		protected IGlobalVector[] currentSolutions;
		protected GenericAnalyzerState[] analyzerStates;

		public IAnalysisWorkflowLog[] Logs { get; set; }

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
				Console.WriteLine("\n\nStaggered step: {0}", staggeredStep);
				previousSolutionNorm = solutionNorm;

				for (int i = 0; i < solvers.Length; i++)
				{
					if (solvers[i].LinearSystem.Solution != null)
					{
						currentSolutions[i] = solvers[i].LinearSystem.Solution.Copy();
					}

					solveMethods()[i]();
				}

				solutionNorm = 0;
				for (int i = 0; i < solvers.Length; i++)
				{
					solutionNorm += solvers[i].LinearSystem.Solution.Norm2();
				}

				error = solutionNorm != 0 ? Math.Abs(solutionNorm - previousSolutionNorm) / solutionNorm : 0;
				Console.WriteLine("Staggered step: {0} - error {1}", staggeredStep, error);
				staggeredStep++;

				if (staggeredStep < maxStaggeredSteps && error > tolerance)
				{
					CreateNewModel(analyzers, solvers);
				}
			}
			while (staggeredStep < maxStaggeredSteps && error > tolerance);
		}

		public abstract void Solve();
	}
}
