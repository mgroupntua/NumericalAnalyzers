namespace MGroup.NumericalAnalyzers
{
	using System;
	using System.Collections.Generic;

	using MGroup.MSolve.AnalysisWorkflow;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.FreedomDegrees;
	using MGroup.MSolve.Discretization.Interfaces;
	using MGroup.MSolve.Logging.Interfaces;
	using MGroup.MSolve.Solution;
	using MGroup.MSolve.Solution.LinearSystems;

	/// <summary>
	/// This class solves the linear system.
	/// </summary>
	public class LinearAnalyzer : IChildAnalyzer
	{
		private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
		private readonly IModel model;
		private readonly IAnalyzerProvider provider;
		private readonly ISolver solver;

		/// <summary>
		/// This class defines the linear anaylzer.
		/// </summary>
		/// <param name="model">Instance of the model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param> 
		public LinearAnalyzer(IModel model, ISolver solver, IAnalyzerProvider provider)
		{
			this.model = model;
			this.solver = solver;
			this.linearSystems = solver.LinearSystems;
			this.provider = provider;
		}

		public Dictionary<int, ILogFactory> LogFactories { get; } = new Dictionary<int, ILogFactory>();

		public Dictionary<int, IAnalyzerLog[]> Logs { get; } = new Dictionary<int, IAnalyzerLog[]>();

		public IParentAnalyzer ParentAnalyzer { get; set; }

		public void BuildMatrices()
		{
			if (ParentAnalyzer == null)
			{
				throw new InvalidOperationException("This linear analyzer has no parent.");
			}

			ParentAnalyzer.BuildMatrices();
		}

		public void Initialize(bool isFirstAnalysis)
		{
			InitializeLogs();
		}

		public void Solve()
		{
			DateTime start = DateTime.Now;
			AddEquivalentNodalLoadsToRHS();
			solver.Solve();
			DateTime end = DateTime.Now;
			StoreLogResults(start, end);
		}

		private void AddEquivalentNodalLoadsToRHS()
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				try
				{
					(INode node, IDofType dof, double displacement) = linearSystem.Subdomain.Constraints.Find(du => du != 0.0);
					double scalingFactor = 1;
					IVector initialFreeSolution = linearSystem.CreateZeroVector();
					IVector equivalentNodalLoads = provider.DirichletLoadsAssembler.GetEquivalentNodalLoads(
						linearSystem.Subdomain, initialFreeSolution, scalingFactor);
					linearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads);
				}
				catch (KeyNotFoundException)
				{
				}
			}
		}

		private void InitializeLogs()
		{
			Logs.Clear();
			foreach (int id in LogFactories.Keys)
			{
				Logs.Add(id, LogFactories[id].CreateLogs());
			}
		}

		private void StoreLogResults(DateTime start, DateTime end)
		{
			foreach (int id in Logs.Keys)
			{
				foreach (var l in Logs[id])
				{
					l.StoreResults(start, end, linearSystems[id].Solution);
				}
			}
		}
	}
}
