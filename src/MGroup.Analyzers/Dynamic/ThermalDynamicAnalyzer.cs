namespace MGroup.Analyzers.Dynamic
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;

	using MGroup.Analyzers.Interfaces;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.Interfaces;
	using MGroup.MSolve.Logging;
	using MGroup.MSolve.Logging.Interfaces;
	using MGroup.Solvers;
	using MGroup.Solvers.LinearSystems;

	/// <summary>
	/// This class makes the appropriate arrangements for the solution of thermal dynamic equations
	/// according to the Central Diffences method.
	/// Authors: Yannis Kalogeris
	/// </summary>
	public class ThermalDynamicAnalyzer : INonLinearParentAnalyzer
	{
		private readonly double beta, timeStep, totalTime;
		private readonly IModel model;
		private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
		private readonly ISolver solver;
		private readonly IImplicitIntegrationProvider provider;
		private Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
		private Dictionary<int, IVector> rhsPrevious = new Dictionary<int, IVector>();
		private Dictionary<int, IVector> temperature = new Dictionary<int, IVector>();
		private Dictionary<int, IVector> capacityTimesTemperature = new Dictionary<int, IVector>();
		private Dictionary<int, IVector> conductivityTimesTemperature = new Dictionary<int, IVector>();

		/// <summary>
		/// Creates an instance that uses a specific problem type and an appropriate child analyzer for the construction of the system of equations arising from the actual physical problem
		/// </summary>
		/// <param name="model">Instance of the model to be solved</param>
		/// <param name="solver">Instance of the solver that will handle the solution of the system of equations</param>
		/// <param name="provider">Instance of the problem type to be solver</param>
		/// <param name="childAnalyzer">Instance of the child analyzer that will handle the solution of the system of equations</param>
		/// <param name="beta">Instance of parameter "beta" of the method that will be initialized</param>
		/// <param name="timeStep">Instance of the time step of the method that will be initialized</param>
		/// <param name="totalTime">Instance of the total time of the method that will be initialized</param>
		public ThermalDynamicAnalyzer(IModel model, ISolver solver, IImplicitIntegrationProvider provider,
			IChildAnalyzer childAnalyzer, double beta, double timeStep, double totalTime)
		{
			this.model = model;
			this.linearSystems = solver.LinearSystems;
			this.solver = solver;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.beta = beta;
			this.timeStep = timeStep;
			this.totalTime = totalTime;
			this.ChildAnalyzer.ParentAnalyzer = this;
		}

		public Dictionary<int, IAnalyzerLog[]> Logs => null;

		public Dictionary<int, ImplicitIntegrationAnalyzerLog> ResultStorages { get; }
			= new Dictionary<int, ImplicitIntegrationAnalyzerLog>();

		public IChildAnalyzer ChildAnalyzer { get; }

		/// <summary>
		/// Makes the proper solver-specific initializations before the solution of the linear system of equations. This method MUST be called before the actual solution of the aforementioned system
		/// </summary>
		public void BuildMatrices()
		{
			var coeffs = new ImplicitIntegrationCoefficients
			{
				Mass = 1 / timeStep,
				Stiffness = beta,
			};
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				linearSystem.Matrix = provider.LinearCombinationOfMatricesIntoStiffness(coeffs, linearSystem.Subdomain);
			}
		}

		/// <summary>
		/// Calculates inertia forces.
		/// </summary>
		public IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution)
		{
			return provider.MassMatrixVectorProduct(linearSystem.Subdomain, currentSolution);
		}

		/// <summary>
		/// Initializes the models, the solvers, child analyzers, builds the matrices, assigns loads and initializes right-hand-side vectors.
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
				}
			}

			BuildMatrices();

			model.AssignLoads(solver.DistributeNodalLoads);
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				linearSystem.RhsVector = linearSystem.Subdomain.Forces;
			}

			InitializeInternalVectors();

			InitializeRhs();

			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
			}

			ChildAnalyzer.Initialize(isFirstAnalysis);
		}

		/// <summary>
		/// Solves the linear system of equations by calling the corresponding method of the specific solver attached during construction of the current instance
		/// </summary>
		public void Solve()
		{
			int numTimeSteps = (int)(totalTime / timeStep);
			for (int t = 0; t < numTimeSteps; ++t)
			{
				Debug.WriteLine("Newmark step: {0}", t);

				IDictionary<int, IVector> rhsVectors = provider.GetRhsFromHistoryLoad(t);
				foreach (var l in linearSystems.Values)
				{
					l.RhsVector = rhsVectors[l.Subdomain.ID];
				}

				InitializeRhs();
				CalculateRhsImplicit();

				DateTime start = DateTime.Now;
				ChildAnalyzer.Solve();
				DateTime end = DateTime.Now;

				UpdateTemperature(t);
				UpdateResultStorages(start, end);
			}
		}

		/// <summary>
		/// Calculates the right-hand-side of the implicit dyanmic method. This will be used for the solution of the linear dynamic system.
		/// </summary>
		private void CalculateRhsImplicit()
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				linearSystem.RhsVector = CalculateRhsImplicit(linearSystem, true);
			}
		}

		private IVector CalculateRhsImplicit(ILinearSystem linearSystem, bool addRhs)
		{
			int id = linearSystem.Subdomain.ID;
			capacityTimesTemperature[id] = provider.MassMatrixVectorProduct(linearSystem.Subdomain, temperature[id]);
			conductivityTimesTemperature[id] = provider.DampingMatrixVectorProduct(linearSystem.Subdomain, temperature[id]);

			IVector rhsResult = rhsPrevious[id].LinearCombination(1 - beta, rhs[id], beta);
			rhsResult.AxpyIntoThis(capacityTimesTemperature[id], 1 / timeStep);
			rhsResult.AxpyIntoThis(conductivityTimesTemperature[id], -(1 - beta));

			rhsPrevious[id] = rhs[id];
			return rhsResult;
		}

		private void InitializeInternalVectors()
		{
			temperature.Clear();
			rhs.Clear();
			rhsPrevious.Clear();
			capacityTimesTemperature.Clear();
			conductivityTimesTemperature.Clear();

			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				capacityTimesTemperature.Add(id, linearSystem.CreateZeroVector());
				conductivityTimesTemperature.Add(id, linearSystem.CreateZeroVector());
				rhs.Add(id, linearSystem.CreateZeroVector());
				rhsPrevious.Add(id, linearSystem.CreateZeroVector());

				if (linearSystem.Solution != null)
				{
					temperature.Add(id, linearSystem.Solution.Copy());
				}
				else
				{
					temperature.Add(id, linearSystem.CreateZeroVector());
				}
			}
		}

		private void InitializeRhs()
		{
			ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
			{
				Mass = 0,
				Stiffness = 0,
			};
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				provider.ProcessRhs(coeffs, linearSystem.Subdomain, linearSystem.RhsVector);
				rhs[linearSystem.Subdomain.ID] = linearSystem.RhsVector.Copy();
			}
		}

		private void UpdateResultStorages(DateTime start, DateTime end)
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				if (ResultStorages.ContainsKey(id))
				{
					if (ResultStorages[id] != null)
					{
						foreach (var l in ChildAnalyzer.Logs[id])
						{
							ResultStorages[id].StoreResults(start, end, l);
						}
					}
				}
			}
		}

		private void UpdateTemperature(int timeStep)
		{
			foreach (ILinearSystem linearSystem in linearSystems.Values)
			{
				int id = linearSystem.Subdomain.ID;
				temperature[id].CopyFrom(linearSystem.Solution);
			}
		}
	}
}
