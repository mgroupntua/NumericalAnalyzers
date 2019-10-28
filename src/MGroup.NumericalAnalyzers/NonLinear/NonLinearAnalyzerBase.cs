namespace MGroup.NumericalAnalyzers.NonLinear
{
	using System;
	using System.Collections.Generic;

	using MGroup.MSolve.AnalysisWorkflow;
	using MGroup.NumericalAnalyzers.Logging;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.Interfaces;
	using MGroup.MSolve.Logging.Interfaces;
	using MGroup.MSolve.Solution;
	using MGroup.MSolve.Solution.LinearSystems;

	/// <summary>
	/// This class represents the base class of all nonlinear anaylsers and contains the basic information necessary for other analyzers
	/// </summary>
	public abstract class NonLinearAnalyzerBase : IChildAnalyzer
	{
		protected readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
		protected readonly int maxIterationsPerIncrement;
		protected readonly IModel model;
		protected readonly int numIncrements;
		protected readonly int numIterationsForMatrixRebuild;
		protected readonly INonLinearProvider provider;
		protected readonly double residualTolerance;
		protected readonly ISolver solver;
		protected readonly IReadOnlyDictionary<int, INonLinearSubdomainUpdater> subdomainUpdaters;
		protected readonly Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
		protected readonly Dictionary<int, IVector> u = new Dictionary<int, IVector>();
		protected readonly Dictionary<int, IVector> du = new Dictionary<int, IVector>();
		protected readonly Dictionary<int, IVector> uPlusdu = new Dictionary<int, IVector>();
		protected Vector globalRhs;
		protected double globalRhsNormInitial;
		protected INonLinearParentAnalyzer parentAnalyzer = null;

		internal NonLinearAnalyzerBase(IModel model, ISolver solver, INonLinearProvider provider,
			IReadOnlyDictionary<int, INonLinearSubdomainUpdater> subdomainUpdaters,
			int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance)
		{
			this.model = model;
			this.solver = solver;
			this.provider = provider;
			this.subdomainUpdaters = subdomainUpdaters;
			this.linearSystems = solver.LinearSystems;
			this.numIncrements = numIncrements;
			this.maxIterationsPerIncrement = maxIterationsPerIncrement;
			this.numIterationsForMatrixRebuild = numIterationsForMatrixRebuild;
			this.residualTolerance = residualTolerance;
		}

		public Dictionary<int, LinearAnalyzerLogFactory> LogFactories { get; } = new Dictionary<int, LinearAnalyzerLogFactory>();

		public Dictionary<int, IAnalyzerLog[]> Logs { get; } = new Dictionary<int, IAnalyzerLog[]>();

		public TotalDisplacementsPerIterationLog TotalDisplacementsPerIterationLog { get; set; }

		public Dictionary<int, TotalLoadsDisplacementsPerIncrementLog> IncrementalLogs { get; }
			= new Dictionary<int, TotalLoadsDisplacementsPerIncrementLog>();

		public IParentAnalyzer ParentAnalyzer
		{
			get => parentAnalyzer;
			set => parentAnalyzer = (INonLinearParentAnalyzer)value;
		}

		/// <summary>
		/// Builds the tangent stiffness matrix of the system.
		/// </summary>
		public void BuildMatrices()
		{
			if (parentAnalyzer == null)
			{
				throw new InvalidOperationException(
				"This Newton-Raphson nonlinear analyzer has no parent.");
			}

			parentAnalyzer.BuildMatrices();
		}

		/// <summary>
		/// Initializes internal vector before the first analysis.
		/// </summary>
		public void Initialize(bool isFirstAnalysis)
		{
			InitializeInternalVectors();
		}

		protected Dictionary<int, IVector> CalculateInternalRhs(int currentIncrement, int iteration)
		{
			var internalRhsVectors = new Dictionary<int, IVector>();
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;

				if (currentIncrement == 0 && iteration == 0)
				{
					du[id].Clear();
					uPlusdu[id].Clear();
					du[id].AddIntoThis(linearSystem.Solution);
					uPlusdu[id].AddIntoThis(linearSystem.Solution);
					du[id].SubtractIntoThis(u[id]);
				}
				else
				{
					du[id].AddIntoThis(linearSystem.Solution);
					uPlusdu[id].Clear();
					uPlusdu[id].AddIntoThis(u[id]);
					uPlusdu[id].AddIntoThis(du[id]);
				}
				IVector internalRhs = subdomainUpdaters[id].GetRhsFromSolution(uPlusdu[id], du[id]);
				provider.ProcessInternalRhs(linearSystem.Subdomain, uPlusdu[id], internalRhs);

				if (parentAnalyzer != null)
				{
					IVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(linearSystem, uPlusdu[id]);
					internalRhs.AddIntoThis(otherRhsComponents);
				}

				internalRhsVectors.Add(id, internalRhs);
			}

			return internalRhsVectors;
		}

		protected double UpdateResidualForcesAndNorm(int currentIncrement, Dictionary<int, IVector> internalRhs)
		{
			globalRhs.Clear();
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;

				linearSystem.RhsVector.Clear();

				for (int j = 0; j <= currentIncrement; j++)
				{
					linearSystem.RhsVector.AddIntoThis(rhs[id]);
				}

				linearSystem.RhsVector.SubtractIntoThis(internalRhs[id]);

				model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
			}
			return provider.CalculateRhsNorm(globalRhs);
		}

		protected void ClearIncrementalSolutionVector()
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				du[linearSystem.Subdomain.ID].Clear();
			}
		}

		protected virtual void InitializeInternalVectors()
		{
			globalRhs = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
			rhs.Clear();
			u.Clear();
			du.Clear();
			uPlusdu.Clear();

			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;

				IVector r = linearSystem.RhsVector.Copy();
				r.ScaleIntoThis(1 / (double)numIncrements);
				rhs.Add(id, r);
				u.Add(id, linearSystem.CreateZeroVector());
				du.Add(id, linearSystem.CreateZeroVector());
				uPlusdu.Add(id, linearSystem.CreateZeroVector());
				model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
			}
			globalRhsNormInitial = provider.CalculateRhsNorm(globalRhs);
		}

		protected void InitializeLogs()
		{
			Logs.Clear();
			foreach (int id in LogFactories.Keys)
			{
				Logs.Add(id, LogFactories[id].CreateLogs());
			}

			foreach (var log in IncrementalLogs.Values)
			{
				log.Initialize();
			}
		}

		protected void SaveMaterialStateAndUpdateSolution()
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				subdomainUpdaters[id].UpdateState();
				u[id].AddIntoThis(du[id]);
			}
		}

		protected void SplitResidualForcesToSubdomains()
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				linearSystem.RhsVector.Clear();
				model.GlobalDofOrdering.ExtractVectorSubdomainFromGlobal(linearSystem.Subdomain, globalRhs,
					linearSystem.RhsVector);
			}
		}

		protected void StoreLogResults(DateTime start, DateTime end)
		{
			foreach (int id in Logs.Keys)
			{
				foreach (var l in Logs[id])
				{
					l.StoreResults(start, end, u[id]);
				}
			}
		}

		protected void UpdateInternalVectors()
		{
			globalRhs.Clear();
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;

				IVector r = linearSystem.RhsVector.Copy();
				r.ScaleIntoThis(1 / (double)numIncrements);
				rhs[id] = r;
				model.GlobalDofOrdering.AddVectorSubdomainToGlobal(linearSystem.Subdomain, linearSystem.RhsVector, globalRhs);
			}
			globalRhsNormInitial = provider.CalculateRhsNorm(globalRhs);
		}

		protected void UpdateRhs(int step)
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				linearSystem.RhsVector.CopyFrom(rhs[linearSystem.Subdomain.ID]);
			}
		}

		/// <summary>
		/// This class solves system and calculates the displacements vector.
		/// </summary>
		public abstract void Solve();
	}
}
