using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Constitutive;

namespace MGroup.NumericalAnalyzers.Dynamic
{
	public class PseudoTransientAnalyzer : INonLinearParentAnalyzer, IStepwiseAnalyzer
	{
		private const string TIME = TransientLiterals.TIME;
		private const string CURRENTTIMESTEP = "Current timestep";
		private const string CURRENTSOLUTION = "Current solution";

		private readonly TransientAnalysisCoefficients transientCoeffs = new ()
		{
			SecondOrderDerivativeCoefficient = 0,
			FirstOrderDerivativeCoefficient = 0,
			ZeroOrderDerivativeCoefficient = 1,
		};

		/// <summary>
		/// This class implements a Pseudo-Transient Analyzer
		/// Authors: George Stavroulakis, Theofilos Christodoulou.
		/// </summary>
		private readonly double timeStep;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Newmark method
		/// Authors: George Stavroulakis, Theofilos Christodoulou.
		/// </summary>
		private readonly double totalTime;
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ISolver solver;
		private readonly ITransientAnalysisProvider provider;
		private IGlobalVector rhs;
		private IGlobalVector solution;
		private int currentStep;
		private DateTime start, end;
		private GenericAnalyzerState currentState;

		/// <summary>
		/// Creates an instance that uses a specific problem type and an appropriate child analyzer for the construction of the system of equations arising from the actual physical problem.
		/// </summary>
		/// <param name="model">Instance of the model to be solved.</param>
		/// <param name="solver">Instance of the solver that will handle the solution of the system of equations.</param>
		/// <param name="provider">Instance of the problem type to be solver.</param>
		/// <param name="childAnalyzer">Instance of the child analyzer that will handle the solution of the system of equations.</param>
		/// <param name="timeStep">Instance of the time step of the method that will be initialized.</param>
		/// <param name="totalTime">Instance of the total time of the method that will be initialized.</param>
		private PseudoTransientAnalyzer(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer, double timeStep, double totalTime, int currentStep)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.timeStep = timeStep;
			this.totalTime = totalTime;
			this.currentStep = currentStep;
			this.ChildAnalyzer.ParentAnalyzer = this;
		}

		public IAnalysisWorkflowLog[] Logs => null;

		public ImplicitIntegrationAnalyzerLog ResultStorage { get; set; }

		public IChildAnalyzer ChildAnalyzer { get; }

		public int CurrentStep { get => currentStep; }

		public int Steps { get => (int)(totalTime / timeStep); }

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				currentStep = (int)currentState.StateValues[CURRENTTIMESTEP];
				currentStep = (int)currentState.StateValues[CURRENTTIMESTEP];

				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;
				solution.CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
			}
		}

		/// <summary>
		/// Makes the proper solver-specific initializations before the solution of the linear system of equations. This method MUST be called before the actual solution of the aforementioned system
		/// </summary>
		public void BuildMatrices()
		{
			provider.LinearCombinationOfMatricesIntoEffectiveMatrix(transientCoeffs);
		}

		/// <summary>
		/// Calculates inertia forces and damping forces.
		/// </summary>
		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution) => algebraicModel.CreateZeroVector();

		/// <summary>
		/// Initializes the models, the solvers, child analyzers, builds the matrices, assigns loads and initializes right-hand-side vectors.
		/// </summary>
		public void Initialize(bool isFirstAnalysis = true)
		{
			if (isFirstAnalysis)
			{
				model.ConnectDataStructures();
				algebraicModel.OrderDofs();
			}

			BuildMatrices();

			provider.AssignRhs();

			InitializeInternalVectors();

			InitializeRhs();

			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Pseudo-Transient analyzer must contain an embedded analyzer.");
			}

			ChildAnalyzer.Initialize(isFirstAnalysis);
		}

		private void SolveCurrentTimestep()
		{
			Debug.WriteLine("Pseudo-Transient Analyzer step: {0}", currentStep);

			IGlobalVector rhsVector = provider.GetRhs(currentStep * timeStep);
			solver.LinearSystem.RhsVector = rhsVector;

			InitializeRhs();
			CalculateRhsImplicit();

			start = DateTime.Now;
			ChildAnalyzer.Initialize(false);
			ChildAnalyzer.Solve();
			end = DateTime.Now;
			Debug.WriteLine("Pseudo-Transient Analyzer elapsed time: {0}", end - start);
		}

		/// <summary>
		/// Solves the linear system of equations by calling the corresponding method of the specific solver attached during construction of the current instance
		/// </summary>
		public void Solve()
		{
			for (int i = 0; i < Steps; ++i)
			{
				SolveCurrentTimestep();
				AdvanceStep();
			}
		}

		/// <summary>
		/// Calculates the right-hand-side of the implicit dyanmic method. This will be used for the solution of the linear dynamic system.
		/// </summary>
		private void CalculateRhsImplicit()
		{
			solver.LinearSystem.RhsVector = rhs;
		}

		private void InitializeInternalVectors()
		{
			rhs = algebraicModel.CreateZeroVector();

			if (solver.LinearSystem.Solution != null)
			{
				solution = solver.LinearSystem.Solution.Copy();
			}
			else
			{
				solution = algebraicModel.CreateZeroVector();
			}
		}

		private void InitializeRhs()
		{
			provider.ProcessRhs(transientCoeffs, solver.LinearSystem.RhsVector);
			rhs.CopyFrom(solver.LinearSystem.RhsVector);
		}

		private void UpdateResultStorages(DateTime start, DateTime end)
		{
			if (ResultStorage != null)
			{
				foreach (var l in ChildAnalyzer.Logs)
				{
					ResultStorage.StoreResults(start, end, l);
				}
			}
		}

		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this,
				new[]
				{
					(CURRENTSOLUTION, solution),
				},
				new[]
				{
					(TIME, (double)currentStep * (double) timeStep),
					(CURRENTTIMESTEP, (double)currentStep),
				});

			return currentState;
		}

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();

		/// <summary>
		/// Solves the linear system of equations of the current timestep
		/// </summary>
		void IStepwiseAnalyzer.Solve()
		{
			SolveCurrentTimestep();
		}

		public void AdvanceStep()
		{
			Debug.WriteLine("Advancing step");
			UpdateResultStorages(start, end);
			currentStep++;
		}

		public class Builder
		{
			private readonly double timeStep;
			private readonly double totalTime;
			private readonly IChildAnalyzer childAnalyzer;
			private readonly IModel model;
			private readonly IAlgebraicModel algebraicModel;
			private readonly ISolver solver;
			private readonly ITransientAnalysisProvider provider;
			private int currentStep = 0;

			public Builder(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
				IChildAnalyzer childAnalyzer, double timeStep, double totalTime, int currentStep = 0)
			{
				this.model = model;
				this.algebraicModel = algebraicModel;
				this.solver = solver;
				this.provider = provider;
				this.childAnalyzer = childAnalyzer;
				this.currentStep = currentStep;

				this.timeStep = timeStep;
				this.totalTime = totalTime;
			}

			public PseudoTransientAnalyzer Build()
				=> new PseudoTransientAnalyzer(model, algebraicModel, solver, provider, childAnalyzer, timeStep, totalTime, currentStep);
		}
	}
}
