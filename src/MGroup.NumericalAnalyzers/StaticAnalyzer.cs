using System;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers
{
	/// <summary>
	/// This class solves static problems.
	/// </summary>
	public class StaticAnalyzer : INonLinearParentAnalyzer
	{
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly INonTransientAnalysisProvider provider;
		private readonly ISolver solver;

		/// <summary>
		/// This class defines the static analyzer.
		/// </summary>
		/// <param name="model">Instance of the model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param> 
		/// <param name="childAnalyzer">Instance of the child analyzer that defines whether the problem is linear or nonlinear</param>
		public StaticAnalyzer(IModel model, IAlgebraicModel algebraicModel, ISolver solver, INonTransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.ChildAnalyzer.ParentAnalyzer = this;
		}

		public IAnalysisWorkflowLog[] Logs { get; set; }

		public IChildAnalyzer ChildAnalyzer { get; }

		/// <summary>
		/// Builds the stiffness matrix of the structure.
		/// </summary>
		public void BuildMatrices()
		{
			provider.CalculateMatrix();
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
				//provider.GetProblemDofTypes();
				model.ConnectDataStructures();
				algebraicModel.OrderDofs();
			}

			BuildMatrices();

			provider.AssignRhs();//AssignLoads(solver.DistributeNodalLoads);
			//foreach (ILinearSystem linearSystem in solver.LinearSystems.Values)
			//{
			//	linearSystem.RhsVector = linearSystem.Subdomain.Forces;
			//}

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

			ChildAnalyzer.Solve();
		}
	}
}
