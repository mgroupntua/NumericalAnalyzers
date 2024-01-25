using System;
using System.Collections.Generic;

using MGroup.LinearAlgebra.Iterative;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers
{
	/// <summary>
	/// This class solves static problems.
	/// </summary>
	public class StaticAnalyzer : INonLinearParentAnalyzer
	{
		//private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly INonTransientAnalysisProvider provider;
		private readonly IList<IterativeStatistics> analysisStatistics = new[]
		{
			new IterativeStatistics() { AlgorithmName = "Static analyzer" },
		};

		/// <summary>
		/// This class defines the static analyzer.
		/// </summary>
		/// <param name="algebraicModel">Instance of the algebraic model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param> 
		/// <param name="childAnalyzer">Instance of the child analyzer that defines whether the problem is linear or nonlinear</param>
		public StaticAnalyzer(IAlgebraicModel algebraicModel, INonTransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer)
		{
			this.algebraicModel = algebraicModel;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.ChildAnalyzer.ParentAnalyzer = this;
		}

		public IAnalysisWorkflowLog[] Logs { get; set; }

		public IGlobalVector CurrentAnalysisResult { get => ChildAnalyzer?.CurrentAnalysisResult; }

		public IChildAnalyzer ChildAnalyzer { get; }

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => CreateState();
			set
			{
			}
		}

		public IList<IterativeStatistics> AnalysisStatistics => analysisStatistics;

		GenericAnalyzerState CreateState() => new GenericAnalyzerState(this, new[]
		{
			(String.Empty, (IGlobalVector)null)
		});

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();

		/// <summary>
		/// Builds the stiffness matrix of the structure.
		/// </summary>
		public void BuildMatrices()
		{
			var matrix = provider.GetMatrix();

			algebraicModel.LinearSystem.Matrix = matrix;
		}

		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution)
		{
			return algebraicModel.CreateZeroVector();
		}

		/// <summary>
		/// Initializes the values of the system to be solved.
		/// </summary>
		public void Initialize(bool isFirstAnalysis = true)
		{
			if (isFirstAnalysis)
			{
				// Connect data structures of model is called by the algebraic model
				algebraicModel.OrderDofs();
			}

			BuildMatrices();
			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
			}

			ChildAnalyzer.Initialize(isFirstAnalysis);
		}

		/// <summary>
		/// Solves the system and calculates the displacement vector.
		/// </summary>
		public void Solve()
		{
			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
			}

			IGlobalVector rhsVector = provider.GetRhs();
			ChildAnalyzer.CurrentAnalysisLinearSystemRhs.AddIntoThis(rhsVector);
			ChildAnalyzer.Solve();
			AnalysisStatistics[0] = ChildAnalyzer.AnalysisStatistics;
		}
	}
}
