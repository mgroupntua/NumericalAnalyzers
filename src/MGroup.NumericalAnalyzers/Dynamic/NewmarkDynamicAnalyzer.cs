using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Constitutive;

namespace MGroup.NumericalAnalyzers.Dynamic
{
	public class NewmarkDynamicAnalyzer : INonLinearParentAnalyzer, IStepwiseAnalyzer
	{
		private const string TIME = TransientLiterals.TIME;
		private const string CURRENTTIMESTEP = "Current timestep";
		private const string CURRENTSOLUTION = "Current solution";
		private const string PREVIOUSSOLUTION = "Previous solution";
		private const string FIRSTORDERSOLUTION = "First order derivative of solution";
		private const string SECONDORDERSOLUTION = "Second order derivative of solution";

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Newmark method
		/// Authors: George Stavroulakis, George Soimiris
		/// </summary>
		private readonly double beta;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Newmark method
		/// Authors: George Stavroulakis, George Soimiris
		/// </summary>
		private readonly double gamma;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Newmark method
		/// Authors: George Stavroulakis, George Soimiris
		/// </summary>
		private readonly double timeStep;

		/// <summary>
		/// This class makes the appropriate arrangements for the solution of linear dynamic equations
		/// according to implicit Newmark method
		/// Authors: George Stavroulakis, George Soimiris
		/// </summary>
		private readonly double totalTime;
		private readonly double a0;
		private readonly double a1;
		private readonly double a2;
		private readonly double a3;
		private readonly double a4;
		private readonly double a5;
		private readonly double a6;
		private readonly double a7;
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ISolver solver;
		private readonly ITransientAnalysisProvider provider;
		private IGlobalVector rhs;
		private IGlobalVector solution;
		private IGlobalVector solutionOfPreviousStep;
		private IGlobalVector firstOrderDerivativeOfSolution;
		private IGlobalVector firstOrderDerivativeOfSolutionForRhs;
		private IGlobalVector firstOrderDerivativeComponentOfRhs;
		private IGlobalVector secondOrderDerivativeOfSolution;
		private IGlobalVector secondOrderDerivativeOfSolutionForRhs;
		private IGlobalVector secondOrderDerivativeComponentOfRhs;
		private int currentStep;
		private DateTime start, end;
		private GenericAnalyzerState currentState;

		/// <summary>
		/// Creates an instance that uses a specific problem type and an appropriate child analyzer for the construction of the system of equations arising from the actual physical problem
		/// </summary>
		/// <param name="model">Instance of the model to be solved</param>
		/// <param name="solver">Instance of the solver that will handle the solution of the system of equations</param>
		/// <param name="provider">Instance of the problem type to be solver</param>
		/// <param name="childAnalyzer">Instance of the child analyzer that will handle the solution of the system of equations</param>
		/// <param name="timeStep">Instance of the time step of the method that will be initialized</param>
		/// <param name="totalTime">Instance of the total time of the method that will be initialized</param>
		/// <param name="alpha">Instance of parameter "alpha" of the method that will be initialized</param>
		/// <param name="delta">Instance of parameter "delta" of the method that will be initialized</param>
		private NewmarkDynamicAnalyzer(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer, double timeStep, double totalTime, double alpha, double delta, int currentStep)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.beta = alpha;
			this.gamma = delta;
			this.timeStep = timeStep;
			this.totalTime = totalTime;
			this.currentStep = currentStep;
			this.ChildAnalyzer.ParentAnalyzer = this;

			/// <summary>
			/// Initialize coefficients. It would make sense for them to be initialized in a different function, if they could
			/// change during the analysis
			/// </summary>
			a0 = 1 / (alpha * timeStep * timeStep);
			a1 = delta / (alpha * timeStep);
			a2 = 1 / (alpha * timeStep);
			a3 = (1 / (2 * alpha)) - 1;
			a4 = (delta / alpha) - 1;
			a5 = timeStep * 0.5 * ((delta / alpha) - 2);
			a6 = timeStep * (1 - delta);
			a7 = delta * timeStep;
		}

		public IAnalysisWorkflowLog[] Logs => null;

		public IGlobalVector CurrentAnalysisResult { get => solution; }

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
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[PREVIOUSSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[SECONDORDERSOLUTION].CheckForCompatibility = false;

				solution.CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);
				solutionOfPreviousStep.CopyFrom(currentState.StateVectors[PREVIOUSSOLUTION]);
				firstOrderDerivativeOfSolution.CopyFrom(currentState.StateVectors[FIRSTORDERSOLUTION]);
				secondOrderDerivativeOfSolution.CopyFrom(currentState.StateVectors[SECONDORDERSOLUTION]);

				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[PREVIOUSSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[FIRSTORDERSOLUTION].CheckForCompatibility = true;
				currentState.StateVectors[SECONDORDERSOLUTION].CheckForCompatibility = true;
			}
		}

		/// <summary>
		/// Makes the proper solver-specific initializations before the solution of the linear system of equations. This method MUST be called before the actual solution of the aforementioned system
		/// </summary>
		public void BuildMatrices()
		{
			var coeffs = new TransientAnalysisCoefficients
			{
				SecondOrderDerivativeCoefficient = a0,
				FirstOrderDerivativeCoefficient = a1,
				ZeroOrderDerivativeCoefficient = 1,
			};
			provider.LinearCombinationOfMatricesIntoEffectiveMatrix(coeffs);
		}

		/// <summary>
		/// Calculates inertia forces and damping forces.
		/// </summary>
		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution)
		{
			IGlobalVector result = provider.SecondOrderDerivativeMatrixVectorProduct(currentSolution);
			IGlobalVector temp = provider.FirstOrderDerivativeMatrixVectorProduct(currentSolution);
			result.LinearCombinationIntoThis(a0, temp, a1);
			return result;
		}

		/// <summary>
		/// Initializes the models, the solvers, child analyzers, builds the matrices, assigns loads and initializes right-hand-side vectors.
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

			provider.AssignRhs();

			InitializeInternalVectors();

			InitializeRhs();

			if (ChildAnalyzer == null)
			{
				throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
			}

			ChildAnalyzer.Initialize(isFirstAnalysis);
		}

		private void SolveCurrentTimestep()
		{
			Debug.WriteLine("Newmark step: {0}", currentStep);

			AddExternalVelocitiesAndAccelerations(currentStep * timeStep);
			IGlobalVector rhsVector = provider.GetRhs(currentStep * timeStep);
			solver.LinearSystem.RhsVector = rhsVector;

			InitializeRhs();
			CalculateRhsImplicit();

			start = DateTime.Now;
			ChildAnalyzer.Initialize(false);
			ChildAnalyzer.Solve();
			end = DateTime.Now;
			Debug.WriteLine("Newmark elapsed time: {0}", end - start);
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
			secondOrderDerivativeOfSolutionForRhs = solution.LinearCombination(a0, firstOrderDerivativeOfSolution, a2);
			secondOrderDerivativeOfSolutionForRhs.AxpyIntoThis(secondOrderDerivativeOfSolution, a3);

			firstOrderDerivativeOfSolutionForRhs = solution.LinearCombination(a1, firstOrderDerivativeOfSolution, a4);
			firstOrderDerivativeOfSolutionForRhs.AxpyIntoThis(secondOrderDerivativeOfSolution, a5);

			secondOrderDerivativeComponentOfRhs = provider.SecondOrderDerivativeMatrixVectorProduct(secondOrderDerivativeOfSolutionForRhs);
			firstOrderDerivativeComponentOfRhs = provider.FirstOrderDerivativeMatrixVectorProduct(firstOrderDerivativeOfSolutionForRhs);

			IGlobalVector rhsResult = secondOrderDerivativeComponentOfRhs.Add(firstOrderDerivativeComponentOfRhs);
			bool addRhs = true;
			if (addRhs)
			{
				rhsResult.AddIntoThis(rhs);
			}

			solver.LinearSystem.RhsVector = rhsResult;
		}

		private void InitializeInternalVectors()
		{
			secondOrderDerivativeOfSolutionForRhs = algebraicModel.CreateZeroVector();
			secondOrderDerivativeComponentOfRhs = algebraicModel.CreateZeroVector();
			firstOrderDerivativeOfSolutionForRhs = algebraicModel.CreateZeroVector();
			firstOrderDerivativeComponentOfRhs = algebraicModel.CreateZeroVector();
			solutionOfPreviousStep = algebraicModel.CreateZeroVector();
			firstOrderDerivativeOfSolution = algebraicModel.CreateZeroVector();
			secondOrderDerivativeOfSolution = algebraicModel.CreateZeroVector();
			rhs = algebraicModel.CreateZeroVector();

			if (ChildAnalyzer?.CurrentAnalysisResult != null)
			{
				solution = ChildAnalyzer.CurrentAnalysisResult.Copy();
			}
			else
			{
				solution = algebraicModel.CreateZeroVector();
			}
			//if (solver.LinearSystem.Solution != null)
			//{
			//	solution = solver.LinearSystem.Solution.Copy();
			//}
			//else
			//{
			//	solution = algebraicModel.CreateZeroVector();
			//}
		}

		private void InitializeRhs()
		{
			TransientAnalysisCoefficients coeffs = new TransientAnalysisCoefficients
			{
				SecondOrderDerivativeCoefficient = a0,
				FirstOrderDerivativeCoefficient = a1,
				ZeroOrderDerivativeCoefficient = 1,
			};
			provider.ProcessRhs(coeffs, solver.LinearSystem.RhsVector);
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

		private void AddExternalVelocitiesAndAccelerations(double time)
		{
			IGlobalVector externalAccelerations = provider.GetSecondOrderDerivativeVectorFromBoundaryConditions(time);
			secondOrderDerivativeOfSolution.AddIntoThis(externalAccelerations);

			IGlobalVector externalVelocities = provider.GetFirstOrderDerivativeVectorFromBoundaryConditions(time);
			firstOrderDerivativeOfSolution.AddIntoThis(externalVelocities);
		}

		private void UpdateVelocityAndAcceleration()
		{
			solutionOfPreviousStep.CopyFrom(solution);
			solution.CopyFrom(ChildAnalyzer.CurrentAnalysisResult);

			var secondOrderDerivativeOfSolutionOfPreviousStep = secondOrderDerivativeOfSolution.Copy();

			secondOrderDerivativeOfSolution = solution.Subtract(solutionOfPreviousStep);
			secondOrderDerivativeOfSolution.LinearCombinationIntoThis(a0, firstOrderDerivativeOfSolution, -a2);
			secondOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, -a3);

			firstOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, a6);
			firstOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolution, a7);
		}
		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this,
				new[]
				{
					(CURRENTSOLUTION, solution),
					(PREVIOUSSOLUTION, solutionOfPreviousStep),
					(FIRSTORDERSOLUTION, firstOrderDerivativeOfSolution),
					(SECONDORDERSOLUTION, secondOrderDerivativeOfSolution),
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

			UpdateVelocityAndAcceleration();
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
			private double beta = 0.25;
			private double gamma = 0.5;
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

			public NewmarkDynamicAnalyzer Build()
				=> new NewmarkDynamicAnalyzer(model, algebraicModel, solver, provider, childAnalyzer, timeStep, totalTime, beta, gamma, currentStep);
		}
	}
}
