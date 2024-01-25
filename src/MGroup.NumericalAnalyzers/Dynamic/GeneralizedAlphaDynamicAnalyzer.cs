using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.NumericalAnalyzers.Logging;
using System.Linq;
using MGroup.LinearAlgebra.Iterative;
using System.Collections.Generic;

namespace MGroup.NumericalAnalyzers.Dynamic
{
	public class GeneralizedAlphaDynamicAnalyzer : INonLinearParentAnalyzer, IStepwiseAnalyzer
	{
		private const string TIME = TransientLiterals.TIME;
		private const string CURRENTTIMESTEP = "Current timestep";
		private const string CURRENTSOLUTION = "Current solution";
		private const string PREVIOUSSOLUTION = "Previous solution";
		private const string FIRSTORDERSOLUTION = "First order derivative of solution";
		private const string SECONDORDERSOLUTION = "Second order derivative of solution";

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Generalized alpha method
		/// Authors: George Stavroulakis
		/// </summary>
		private readonly double beta;
		private readonly bool calculateInitialDerivativeVectors = true;
		private TransientAnalysisPhase analysisPhase = TransientAnalysisPhase.InitialConditionEvaluation;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Generalized alpha method
		/// Authors: George Stavroulakis
		/// </summary>
		private readonly double gamma;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Generalized alpha method
		/// Authors: George Stavroulakis
		/// </summary>
		private readonly double timeStep;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Generalized alpha method
		/// Authors: George Stavroulakis
		/// </summary>
		private int currentStep;
		private DateTime start, end;
		private readonly double totalTime;
		private readonly double am, af;
		private readonly double a0, a0N;
		private readonly double a1;
		private readonly double a2, a2N;
		private readonly double a3, a3N;
		private readonly double a4;
		private readonly double a5;
		private readonly double a6N;
		private readonly double a7N;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ITransientAnalysisProvider provider;
		private IGlobalVector[] solutions;
		private IGlobalVector rhs;
		private IGlobalVector zeroOrderDerivativeSolutionOfPreviousStep;
		private IGlobalVector firstOrderDerivativeOfSolutionForRhs;
		private IGlobalVector firstOrderDerivativeComponentOfRhs;
		private IGlobalVector secondOrderDerivativeOfSolutionForRhs;
		private IGlobalVector secondOrderDerivativeComponentOfRhs;
		private GenericAnalyzerState currentState;
		private IList<IterativeStatistics> analysisStatistics;

		/// <summary>
		/// Creates an instance that uses a specific problem type and an appropriate child analyzer for the construction of the system of equations arising from the actual physical problem
		/// </summary>
		/// <param name="model">Instance of the model to be solved</param>
		/// <param name="provider">Instance of the problem type to be solver</param>
		/// <param name="childAnalyzer">Instance of the child analyzer that will handle the solution of the system of equations</param>
		/// <param name="timeStep">Instance of the time step of the method that will be initialized</param>
		/// <param name="totalTime">Instance of the total time of the method that will be initialized</param>
		/// <param name="alpha">Instance of parameter "alpha" of the method that will be initialized</param>
		/// <param name="delta">Instance of parameter "delta" of the method that will be initialized</param>
		/// <param name="currentStep">Starts the analysis from step equal to this parameter</param>
		/// <param name="calculateInitialDerivativeVectors">Set to false to skip initial condition calculation based on initial values (default is true)</param>
		private GeneralizedAlphaDynamicAnalyzer(IAlgebraicModel algebraicModel, ITransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer, double timeStep, double totalTime, double alpha, double delta, double am, double af, int currentStep, bool calculateInitialDerivativeVectors)
		{
			this.calculateInitialDerivativeVectors = calculateInitialDerivativeVectors;
			this.algebraicModel = algebraicModel;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.beta = alpha;
			this.gamma = delta;
			this.am = am;
			this.af = af;
			this.timeStep = timeStep;
			this.totalTime = totalTime;
			this.currentStep = currentStep;
			this.ChildAnalyzer.ParentAnalyzer = this;

			/// <summary>
			/// Initialize coefficients. It would make sense for them to be initialized in a different function, if they could
			/// change during the analysis
			/// </summary>
			a0N = 1 / (alpha * timeStep * timeStep);
			a0 = (1 - am) / (alpha * timeStep * timeStep);
			a1 = delta * (1 - af) / (alpha * timeStep);
			a2N = 1 / (alpha * timeStep);
			a2 = (1 - am) / (alpha * timeStep);
			a3N = (1 / (2 * alpha)) - 1;
			a3 = ((1 - am) / (2 * alpha)) - 1;
			a4 = ((delta - delta * af) / alpha) - 1;
			a5 = (1 - af) * timeStep * 0.5 * ((delta / alpha) - 2);
			a6N = timeStep * (1 - delta);
			a7N = delta * timeStep;

			if (provider.ProblemOrder > DifferentiationOrder.Second)
			{
				throw new ArgumentException($"Wrong problem order. Must be zero, first or second order and it is {provider.ProblemOrder}");
			}

			this.calculateInitialDerivativeVectors = calculateInitialDerivativeVectors;
			this.analysisStatistics = Enumerable.Range(0, Steps).Select(x => new IterativeStatistics() { AlgorithmName = "Generalized alpha dynamic analyzer" }).ToArray();
		}

		public IAnalysisWorkflowLog[] Logs => null;

		public IGlobalVector CurrentAnalysisResult { get => solutions[0]; }

		public ImplicitIntegrationAnalyzerLog ResultStorage { get; set; }

		public IChildAnalyzer ChildAnalyzer { get; }

		public int CurrentStep { get => currentStep; }

		public int Steps { get => (int)(totalTime / timeStep); }

		public IList<IterativeStatistics> AnalysisStatistics => analysisStatistics;

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				currentStep = (int)currentState.StateValues[CURRENTTIMESTEP];
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[PREVIOUSSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[SECONDORDERSOLUTION].CheckForCompatibility = false;

				solutions[0].CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);
				zeroOrderDerivativeSolutionOfPreviousStep.CopyFrom(currentState.StateVectors[PREVIOUSSOLUTION]);
				solutions[(int)DifferentiationOrder.First].CopyFrom(currentState.StateVectors[FIRSTORDERSOLUTION]);
				solutions[(int)DifferentiationOrder.Second].CopyFrom(currentState.StateVectors[SECONDORDERSOLUTION]);

				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[PREVIOUSSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[SECONDORDERSOLUTION].CheckForCompatibility = true;
			}
		}

		private void SetAnalysisPhase(TransientAnalysisPhase phase)
		{
			analysisPhase = phase;
			provider.SetTransientAnalysisPhase(phase);
		}

		/// <summary>
		/// Makes the proper solver-specific initializations before the solution of the linear system of equations. This method MUST be called before the actual solution of the aforementioned system
		/// </summary>
		public void BuildMatrices()
		{
			var matrix = provider.GetMatrix(DifferentiationOrder.Zero).Copy();
			matrix.LinearCombinationIntoThis(1d, provider.GetMatrix(DifferentiationOrder.Second), a0);
			matrix.AxpyIntoThis(provider.GetMatrix(DifferentiationOrder.First), a1);

			algebraicModel.LinearSystem.Matrix = matrix;
		}

		/// <summary>
		/// Calculates equivalent right-hand side for first- and second-order time derivatives for use in non-linear solvers.
		/// Returns zero vector if transient analysis phase is TransientAnalysisPhase.InitialConditionEvaluation.
		/// </summary>
		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution)
		{
			var result = algebraicModel.CreateZeroVector();
			if (analysisPhase == TransientAnalysisPhase.InitialConditionEvaluation)
			{
				return result;
			}

			provider.GetMatrix(DifferentiationOrder.Second).MultiplyVector(currentSolution, result);
			var temp = algebraicModel.CreateZeroVector();
			provider.GetMatrix(DifferentiationOrder.First).MultiplyVector(currentSolution, temp);
			result.LinearCombinationIntoThis(a0, temp, a1);

			return result;
		}

		private void SolveForInitialConditions()
		{
			if (calculateInitialDerivativeVectors == false || currentStep != 0 || provider.ProblemOrder == DifferentiationOrder.Zero)
			{
				return;
			}

			SetAnalysisPhase(TransientAnalysisPhase.InitialConditionEvaluation);
			var rhsFromDerivatives = algebraicModel.CreateZeroVector();
			var temp = algebraicModel.CreateZeroVector();
			for (int i = 0; i < (int)provider.ProblemOrder - 1; i++)
			{
				var d = (DifferentiationOrder)i;
				var lhs = provider.GetVectorFromModelConditions(d, 0);
				if (lhs.Norm2() != 0)
				{
					provider.GetMatrix(d).MultiplyVector(lhs, temp);
					rhsFromDerivatives.SubtractIntoThis(temp);
				}
			}

			algebraicModel.LinearSystem.Matrix = provider.GetMatrix(provider.ProblemOrder).Copy();
			var rhsVector = provider.GetRhs(currentStep * timeStep);
			ChildAnalyzer.CurrentAnalysisLinearSystemRhs.CopyFrom(rhsVector);
			ChildAnalyzer.CurrentAnalysisLinearSystemRhs.AddIntoThis(rhsFromDerivatives);

			ChildAnalyzer.Initialize(false);
			ChildAnalyzer.Solve();

			solutions[(int)provider.ProblemOrder].CopyFrom(algebraicModel.LinearSystem.Solution);
		}

		/// <summary>
		/// Initializes the models, the solvers, child analyzers, builds the matrices, solves for initial values and initializes right-hand-side vectors.
		/// </summary>
		public void Initialize(bool isFirstAnalysis = true)
		{
			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Generalized alpha analyzer must contain an embedded analyzer.");
			}

			if (isFirstAnalysis)
			{
				// Connect data structures of model is called by the algebraic model
				algebraicModel.OrderDofs();
			}

			ChildAnalyzer.Initialize(isFirstAnalysis);
			InitializeInternalVectors();
			SolveForInitialConditions();

			BuildMatrices();
			InitializeRhs();
		}

		private void SolveCurrentTimestep()
		{
			SetAnalysisPhase(TransientAnalysisPhase.Solution);
			Debug.WriteLine("Generalized alpha step: {0}", currentStep);

			AddHigherOrderContributions(currentStep * timeStep);
			IGlobalVector rhsVector = provider.GetRhs(currentStep * timeStep);
			ChildAnalyzer.CurrentAnalysisLinearSystemRhs.CopyFrom(rhsVector);

			InitializeRhs();
			CalculateRhsImplicit();

			start = DateTime.Now;
			ChildAnalyzer.Initialize(false);
			ChildAnalyzer.Solve();
			analysisStatistics[currentStep] = ChildAnalyzer.AnalysisStatistics;
			end = DateTime.Now;
			Debug.WriteLine("Generalized alpha elapsed time: {0}", end - start);
		}

		/// <summary>
		/// Perform the transient analysis by employing the assigned child analyzer for every timestep.
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
			if (provider.ProblemOrder == DifferentiationOrder.Second)
			{
				secondOrderDerivativeOfSolutionForRhs = solutions[0].LinearCombination(a0, solutions[(int)DifferentiationOrder.First], a2);
				secondOrderDerivativeOfSolutionForRhs.AxpyIntoThis(solutions[(int)DifferentiationOrder.Second], a3);
				secondOrderDerivativeComponentOfRhs.Clear();
				provider.GetMatrix(DifferentiationOrder.Second).MultiplyVector(secondOrderDerivativeOfSolutionForRhs, secondOrderDerivativeComponentOfRhs);

				firstOrderDerivativeOfSolutionForRhs = solutions[0].LinearCombination(a1, solutions[(int)DifferentiationOrder.First], a4);
				firstOrderDerivativeOfSolutionForRhs.AxpyIntoThis(solutions[(int)DifferentiationOrder.Second], a5);
				firstOrderDerivativeComponentOfRhs.Clear();
				provider.GetMatrix(DifferentiationOrder.First).MultiplyVector(firstOrderDerivativeOfSolutionForRhs, firstOrderDerivativeComponentOfRhs);
			}

			if (provider.ProblemOrder == DifferentiationOrder.First)
			{
				firstOrderDerivativeOfSolutionForRhs = solutions[0].LinearCombination(a1, solutions[(int)DifferentiationOrder.First], a4);
				firstOrderDerivativeOfSolutionForRhs.AxpyIntoThis(solutions[(int)DifferentiationOrder.Second], a5);
				provider.GetMatrix(DifferentiationOrder.First).MultiplyVector(firstOrderDerivativeOfSolutionForRhs, firstOrderDerivativeComponentOfRhs);
			}

			IGlobalVector rhsResult = secondOrderDerivativeComponentOfRhs.Add(firstOrderDerivativeComponentOfRhs);
			bool addRhs = true;
			if (addRhs)
			{
				rhsResult.AddIntoThis(rhs);
			}

			ChildAnalyzer.CurrentAnalysisLinearSystemRhs.CopyFrom(rhsResult);
		}

		private void InitializeInternalVectors()
		{
			solutions = new IGlobalVector[3];
			for (int i = 0; i < 3; i++)
			{
				solutions[i] = algebraicModel.CreateZeroVector();
			}

			secondOrderDerivativeOfSolutionForRhs = algebraicModel.CreateZeroVector();
			secondOrderDerivativeComponentOfRhs = algebraicModel.CreateZeroVector();
			firstOrderDerivativeOfSolutionForRhs = algebraicModel.CreateZeroVector();
			firstOrderDerivativeComponentOfRhs = algebraicModel.CreateZeroVector();
			zeroOrderDerivativeSolutionOfPreviousStep = algebraicModel.CreateZeroVector();
			rhs = algebraicModel.CreateZeroVector();

			//Code in comments in case we need replicate previous behavior
			//if (ChildAnalyzer?.CurrentAnalysisResult != null)
			//{
			//	solutions[0] = ChildAnalyzer.CurrentAnalysisResult.Copy();
			//}
			//else
			//{
			//	solutions[0] = provider.GetVectorFromModelConditions(DifferentiationOrder.Zero, 0);
			//	//TODO na upologisthoun solutions[1] = provider.
			//}
			solutions[0] = provider.GetVectorFromModelConditions(DifferentiationOrder.Zero, 0);
		}

		private void InitializeRhs()
		{
			rhs.CopyFrom(ChildAnalyzer.CurrentAnalysisLinearSystemRhs);
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

		private void AddHigherOrderContributions(double time)
		{
			for (int i = 1; i <= (int)provider.ProblemOrder; i++)
			{
				solutions[i].AddIntoThis(provider.GetVectorFromModelConditions((DifferentiationOrder)i, time));
			}
		}

		private void UpdateHigherOrderDerivatives()
		{
			zeroOrderDerivativeSolutionOfPreviousStep.CopyFrom(solutions[0]);
			solutions[0].CopyFrom(ChildAnalyzer.CurrentAnalysisResult);

			var secondOrderDerivativeOfSolutionOfPreviousStep = solutions[(int)DifferentiationOrder.Second].Copy();

			solutions[(int)DifferentiationOrder.Second] = solutions[0].Subtract(zeroOrderDerivativeSolutionOfPreviousStep);
			solutions[(int)DifferentiationOrder.Second].LinearCombinationIntoThis(a0N, solutions[(int)DifferentiationOrder.First], -a2N);
			solutions[(int)DifferentiationOrder.Second].AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, -a3N);

			solutions[(int)DifferentiationOrder.First].AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, a6N);
			solutions[(int)DifferentiationOrder.First].AxpyIntoThis(solutions[(int)DifferentiationOrder.Second], a7N);
		}

		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this,
				new[]
				{
					(CURRENTSOLUTION, solutions[0]),
					(PREVIOUSSOLUTION, zeroOrderDerivativeSolutionOfPreviousStep),
					(FIRSTORDERSOLUTION, solutions[(int)DifferentiationOrder.First]),
					(SECONDORDERSOLUTION, solutions[(int)DifferentiationOrder.Second]),
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

			UpdateHigherOrderDerivatives();
			UpdateResultStorages(start, end);

			currentStep++;
		}

		public class Builder
		{
			private readonly bool calculateInitialDerivativeVectors = true;
			private readonly double timeStep;
			private readonly double totalTime;
			private readonly IChildAnalyzer childAnalyzer;
			private readonly IAlgebraicModel algebraicModel;
			private readonly ITransientAnalysisProvider provider;
			private double am = 0;
			private double af = 0;
			private double beta = 0.25;
			private double gamma = 0.5;
			private int currentStep = 0;

			public Builder(IAlgebraicModel algebraicModel, ITransientAnalysisProvider provider,
				IChildAnalyzer childAnalyzer, double timeStep, double totalTime, bool calculateInitialDerivativeVectors = true, int currentStep = 0)
			{
				this.calculateInitialDerivativeVectors = calculateInitialDerivativeVectors;
				this.algebraicModel = algebraicModel;
				this.provider = provider;
				this.childAnalyzer = childAnalyzer;
				this.currentStep = currentStep;

				this.timeStep = timeStep;
				this.totalTime = totalTime;
			}

			public void SetSpectralRadius(double spectralRadius)
			{
				this.am = (2 * spectralRadius - 1) / (spectralRadius + 1);
				this.af = spectralRadius / (spectralRadius + 1);
				this.beta = 0.25 * Math.Pow(1 - am + af, 2);
				this.gamma = 0.5 - am + af;
			}

			/// <summary>
			///
			/// </summary>
			/// <param name="beta">
			/// Used in the intepolation between the accelerations of the previous and current time step, in order to obtain the
			/// current displacements. Also called alpha by Bathe.
			/// </param>
			/// <param name="gamma">
			/// Used in the intepolation between the accelerations of the previous and current time step, in order to obtain the
			/// current velocities. Also called delta by Bathe.
			/// </param>
			/// <param name="allowConditionallyStable">
			/// If set to true, the user must make sure that the time step chosen is lower than the critical step size
			/// corresponding to these particular <paramref name="beta"/>, <paramref name="gamma"/> parameters.
			/// </param>
			public void SetNewmarkParameters(double beta, double gamma, bool allowConditionallyStable = false)
			{
				if (!allowConditionallyStable)
				{
					if (gamma < 0.5)
					{
						throw new ArgumentException(
						"Newmark delta has to be bigger than 0.5 to ensure unconditional stability.");
					}

					if (beta < 0.25)
					{
						throw new ArgumentException(
						"Newmark alpha has to be bigger than 0.25 to ensure unconditional stability.");
					}
				}
				if (gamma < 0.5)
				{
					throw new ArgumentException("Newmark delta has to be bigger than 0.5.");
				}

				double aLimit = 0.25 * Math.Pow(0.5 + gamma, 2);
				if (beta < aLimit)
				{
					throw new ArgumentException($"Newmark alpha has to be bigger than {aLimit}.");
				}

				this.gamma = gamma;
				this.beta = beta;
			}

			/// <summary>
			/// Central diffences: gamma = 1/2, beta = 0. Newmark results in central diffences, a conditionally stable explicit
			/// method. To ensure stability, the time step must be &lt;= the critical step size = 2 / w,  where w is the maximum
			/// natural radian frequency. It would be more efficient to use an explicit dynamic analyzer.
			/// </summary>
			public void SetNewmarkParametersForCentralDifferences()
			{
				gamma = 0.5;
				beta = 0.0;
			}

			/// <summary>
			/// Constant acceleration (also called average acceleration or trapezoid rule): gamma = 1/2, beta = 1/4.
			/// This is the most common scheme and is unconditionally stable. In this analyzer, it is used as the default.
			/// </summary>
			public void SetNewmarkParametersForConstantAcceleration()
			{
				gamma = 0.5;
				beta = 0.25;
			}

			/// <summary>
			/// Linear acceleration: gamma = 1/2, beta = 1/6. This is more accurate than the default constant acceleration,
			/// but it conditionally stable. To ensure stability, the time step must be &lt;= the critical step size = 3.464 / w
			/// = 0.551 * T, where w is the maximum natural radian frequency and T is the minimum natural period.
			/// </summary>
			public void SetNewmarkParametersForLinearAcceleration()
			{
				gamma = 0.5;
				beta = 1.0 / 6.0;
			}

			public GeneralizedAlphaDynamicAnalyzer Build()
				=> new GeneralizedAlphaDynamicAnalyzer(algebraicModel, provider, childAnalyzer, timeStep, totalTime, beta, gamma, am, af, currentStep, calculateInitialDerivativeVectors);
		}
	}
}
