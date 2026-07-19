using System;
using System.Collections.Generic;

using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.AlgebraicModel;

namespace MGroup.NumericalAnalyzers
{
	/// <summary>
	/// This class solves the linear system.
	/// </summary>
	public class LinearAnalyzer : IChildAnalyzer
	{
		private readonly IAlgebraicModel algebraicModel;
		private readonly IAnalyzerProvider provider;
		private readonly ISolver solver;
		private IterativeStatistics analysisStatistics = new IterativeStatistics()
		{
			AlgorithmName = "Linear analyzer",
		};

		/// <summary>
		/// This class defines the linear anaylzer.
		/// </summary>
		/// <param name="algebraicModel">Instance of the algebraic model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param> 
		public LinearAnalyzer(IAlgebraicModel algebraicModel, ISolver solver, IAnalyzerProvider provider)
		{
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
		}

		public ILogFactory LogFactory { get; set; }

		public IVector CurrentAnalysisResult { get => solver?.LinearSystem?.Solution; }

		public IAnalysisWorkflowLog[] Logs { get; set; } = new IAnalysisWorkflowLog[0];

		public IParentAnalyzer ParentAnalyzer { get; set; }

		public IVector Responses { get; set; }

		public IVector CurrentAnalysisLinearSystemRhs { get => solver.LinearSystem.RhsVector; }

		public IterativeStatistics AnalysisStatistics => analysisStatistics; 
		
		public void Initialize(bool isFirstAnalysis)
		{
			//if (isFirstAnalysis)
			//	provider.GetProblemDofTypes();
			InitializeLogs();
			AddEquivalentNodalLoadsToRHS();
		}

		public void Solve()
		{
			var start = DateTime.Now;
			solver.Solve();
			
			analysisStatistics.HasConverged = true;
			analysisStatistics.NumIterationsRequired = 1;
			Responses = solver.LinearSystem.Solution.Copy();
			var end = DateTime.Now;

			StoreLogResults(start, end);
		}

		private void AddEquivalentNodalLoadsToRHS()
		{
			//TODO: equivalentNodalLoads will often be 0. Perhaps instead of AddToGlobalVector, we should have AxpyToGlobalVector
			IVector equivalentNodalLoads = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(provider.EnumerateEquivalentNeumannBoundaryConditions, equivalentNodalLoads);
			solver.LinearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads); 
		}

		private void InitializeLogs()
		{
			if (LogFactory != null)
			{
				Logs = LogFactory.CreateLogs();
			}
		}

		private void StoreLogResults(DateTime start, DateTime end)
		{
			foreach (var l in Logs)
			{
				l.StoreResults(start, end, solver.LinearSystem.Solution);
			}
		}

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => CreateState();
			set
			{
			}
		}


		GenericAnalyzerState CreateState() => new GenericAnalyzerState(this, new[]
		{
			(String.Empty, (IVector)null)
		});

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();
	}
}
