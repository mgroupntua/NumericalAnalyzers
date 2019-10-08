namespace MGroup.Analyzers
{
	using System;
	using System.Collections.Generic;

	using MGroup.Analyzers.Interfaces;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.Interfaces;
	using MGroup.MSolve.Logging.Interfaces;
	using MGroup.Solvers;
	using MGroup.Solvers.LinearSystems;

	/// <summary>
	/// This class solves static problems.
	/// </summary>
	public class StaticAnalyzer : INonLinearParentAnalyzer
	{
		private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
		private readonly IModel model;
		private readonly IStaticProvider provider;
		private readonly ISolver solver;

		/// <summary>
		/// This class defines the static analyzer.
		/// </summary>
		/// <param name="model">Instance of the model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param> 
		/// <param name="childAnalyzer">Instance of the child analyzer that defines whether the problem is linear or nonlinear</param>
		public StaticAnalyzer(IModel model, ISolver solver, IStaticProvider provider,
			IChildAnalyzer childAnalyzer)
		{
			this.model = model;
			this.linearSystems = solver.LinearSystems;
			this.solver = solver;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.ChildAnalyzer.ParentAnalyzer = this;
		}

		public Dictionary<int, IAnalyzerLog[]> Logs { get; } = new Dictionary<int, IAnalyzerLog[]>();

		public IChildAnalyzer ChildAnalyzer { get; }

		/// <summary>
		/// Builds the stiffness matrix of the structure.
		/// </summary>
		public void BuildMatrices()
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				linearSystem.Matrix = provider.CalculateMatrix(linearSystem.Subdomain);
			}
		}

		/// <summary>
		/// Calculates other components of the right-hand-side vector
		/// </summary>
		public IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution)
		{
			return linearSystem.CreateZeroVector();
		}

		/// <summary>
		/// Initializes the values of the system to be solved.
		/// </summary>
		public void Initialize(bool isFirstAnalysis = true)
		{
			if (isFirstAnalysis)
			{
				model.ConnectDataStructures();
				solver.OrderDofs(false);
				foreach (ILinearSystem linearSystem in linearSystems.Values)
				{
					linearSystem.Reset();
					linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
				}
			}
			else
			{
				foreach (ILinearSystem linearSystem in linearSystems.Values)
				{
					linearSystem.Reset();
					linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
				}
			}

			BuildMatrices();

			model.AssignLoads(solver.DistributeNodalLoads);
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				linearSystem.RhsVector = linearSystem.Subdomain.Forces;
			}

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
