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

namespace MGroup.NumericalAnalyzers.Dynamic
{
	public class CentralDifferences : INonLinearParentAnalyzer
	{
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
		private IGlobalVector zeroOrderDerivativeOfSolutionForRhs;
		private IGlobalVector zeroOrderDerivativeComponentOfRhs;
		private IGlobalVector firstOrderDerivativeOfSolution;
		private IGlobalVector firstOrderDerivativeOfSolutionForRhs;
		private IGlobalVector firstOrderDerivativeComponentOfRhs;
		private IGlobalVector secondOrderDerivativeOfSolution;
		private IGlobalVector secondOrderDerivativeOfSolutionForRhs;
		private IGlobalVector secondOrderDerivativeComponentOfRhs;

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
		private CentralDifferences(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
			IChildAnalyzer childAnalyzer, double timeStep, double totalTime, double alpha, double delta)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.ChildAnalyzer = childAnalyzer;
			this.timeStep = timeStep;
			this.totalTime = totalTime;
			this.ChildAnalyzer.ParentAnalyzer = this;

			/// <summary>
			/// Initialize coefficients. It would make sense for them to be initialized in a different function, if they could
			/// change during the analysis
			/// </summary>
			a0 = 1.0 / (timeStep * timeStep);
			a1 = 1.0 / (2.0 * timeStep);
			a2 = 2.0 * a0;
			a3 = 1.0 / a2;
		}

		public IAnalysisWorkflowLog[] Logs => null;

		public ImplicitIntegrationAnalyzerLog ResultStorage { get; set; }

		public IChildAnalyzer ChildAnalyzer { get; }

		/// <summary>
		/// Makes the proper solver-specific initializations before the solution of the linear system of equations. This method MUST be called before the actual solution of the aforementioned system
		/// </summary>
		public void BuildMatrices()
		{
			var coeffs = new TransientAnalysisCoefficients
			{
				SecondOrderDerivativeCoefficient = a0,
				FirstOrderDerivativeCoefficient = a1,
				ZeroOrderDerivativeCoefficient = 0,
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

		/// <summary>
		/// Solves the linear system of equations by calling the corresponding method of the specific solver attached during construction of the current instance
		/// </summary>
		public void Solve()
		{
			int numTimeSteps = (int)(totalTime / timeStep);
			for (int i = 0; i < numTimeSteps; ++i)
			{
				Debug.WriteLine("Newmark step: {0}", i);

				AddExternalVelocitiesAndAccelerations(i * timeStep);
				IGlobalVector rhsVector = provider.GetRhs(i * timeStep);
				solver.LinearSystem.RhsVector = rhsVector; //TODOGoat: Perhaps the provider should set the rhs vector, like it does for the matrix. Either way the provider does this as a side effect

				InitializeRhs();
				CalculateRhsImplicit();

				DateTime start = DateTime.Now;
				ChildAnalyzer.Solve();
				DateTime end = DateTime.Now;

				UpdateVelocityAndAcceleration();
				UpdateResultStorages(start, end);
			}
		}

		/// <summary>
		/// Calculates the right-hand-side of the implicit dyanmic method. This will be used for the solution of the linear dynamic system.
		/// </summary>
		private void CalculateRhsImplicit()
		{
			secondOrderDerivativeOfSolutionForRhs = solution.LinearCombination(a2, solutionOfPreviousStep, -a0);

			firstOrderDerivativeOfSolutionForRhs = solutionOfPreviousStep.Scale(a1);

			zeroOrderDerivativeOfSolutionForRhs = solution.Scale(-1.0);

			secondOrderDerivativeComponentOfRhs = provider.SecondOrderDerivativeMatrixVectorProduct(secondOrderDerivativeOfSolutionForRhs);
			firstOrderDerivativeComponentOfRhs = provider.FirstOrderDerivativeMatrixVectorProduct(firstOrderDerivativeOfSolutionForRhs);
			zeroOrderDerivativeComponentOfRhs = provider.ZeroOrderDerivativeMatrixVectorProduct(zeroOrderDerivativeOfSolutionForRhs);

			IGlobalVector rhsResult = secondOrderDerivativeComponentOfRhs.Add(firstOrderDerivativeComponentOfRhs);
			rhsResult.AddIntoThis(zeroOrderDerivativeComponentOfRhs);
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
			zeroOrderDerivativeOfSolutionForRhs = algebraicModel.CreateZeroVector();
			zeroOrderDerivativeComponentOfRhs = algebraicModel.CreateZeroVector();
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
			IGlobalVector solutionBeforePreviousStep = solutionBeforePreviousStep.CopyFrom(solutionOfPreviousStep);
			solutionOfPreviousStep.CopyFrom(solution);
			solution.CopyFrom(solver.LinearSystem.Solution);

			

			secondOrderDerivativeOfSolution = solution.Subtract(solutionOfPreviousStep);
			secondOrderDerivativeOfSolution.LinearCombinationIntoThis(a0, firstOrderDerivativeOfSolution, -a2);
			secondOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, -a3);

			firstOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, a6);
			firstOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolution, a7);
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

			public Builder(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
				IChildAnalyzer childAnalyzer, double timeStep, double totalTime)
			{
				this.model = model;
				this.algebraicModel = algebraicModel;
				this.solver = solver;
				this.provider = provider;
				this.childAnalyzer = childAnalyzer;

				this.timeStep = timeStep;
				this.totalTime = totalTime;
			}

			
		}
	}
}
