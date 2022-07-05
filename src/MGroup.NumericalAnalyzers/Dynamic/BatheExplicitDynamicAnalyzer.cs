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
	public class BatheExplicitDynamicAnalyzer : INonLinearParentAnalyzer
	{
		public IChildAnalyzer ChildAnalyzer { get; }

		public IAnalysisWorkflowLog[] Logs => null;
		
		private readonly double timeStep;

		private readonly double totalTime;
		
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

		private BatheExplicitDynamicAnalyzer(IModel model, IAlgebraicModel algebraicModel, ISolver solver, ITransientAnalysisProvider provider,
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
		}
		private List<double> Calculate_qValues(double p)
		{
			double q1 = (1.0 - 2.0 * p) / (2.0 * p * (1.0 - p));
			double q2 = 0.5 - p * q1;
			double q0 = -q1 - q2 + 0.5;
			List<double> q = new List<double>();
			q.Add(q0);
			q.Add(q1);
			q.Add(q2);
			return q;
		}

		private List<double> Calculate_aValues(double p, List<double> q, double tStep)
		{
			double q0 = q[0];
			double q1 = q[1];
			double q2 = q[2];
			double deltat = tStep;
			double a0 = p * deltat;
			double a1 = 0.5 * Math.Pow(p * deltat, 2);
			double a2 = a0 / 2.0;
			double a3 = (1.0 - p) * deltat;
			double a4 = 0.5 * Math.Pow((1.0 - p) * deltat, 2);
			double a5 = q0 * a3;
			double a6 = (0.5 + q1) * a3;
			double a7 = q2 * a3;
			List<double> a = new List<double>();
			a.Add(a0);
			a.Add(a1);
			a.Add(a2);
			a.Add(a3);
			a.Add(a4);
			a.Add(a5);
			a.Add(a6);
			a.Add(a7);
			return a;
		}
		public void BuildMatrices()
		{
			throw new NotImplementedException();
		}

		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution)
		{
			throw new NotImplementedException();
		}

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


		private void UpdateVelocityAndAcceleration()
		{
			solutionOfPreviousStep.CopyFrom(solution);
			solution.CopyFrom(solver.LinearSystem.Solution);

			var secondOrderDerivativeOfSolutionOfPreviousStep = secondOrderDerivativeOfSolution.Copy();

			secondOrderDerivativeOfSolution = solution.Subtract(solutionOfPreviousStep);
			secondOrderDerivativeOfSolution.LinearCombinationIntoThis(a0, firstOrderDerivativeOfSolution, -a2);
			secondOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, -a3);

			firstOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolutionOfPreviousStep, a6);
			firstOrderDerivativeOfSolution.AxpyIntoThis(secondOrderDerivativeOfSolution, a7);
		}

		private IGlobalVector U_middle(IGlobalVector u_previous, IGlobalVector du_previous, IGlobalVector ddu_previous, double a0, double a1)
		{
			IGlobalVector u_middle = algebraicModel.CreateZeroVector();
			u_middle = du_previous.LinearCombination(a0, ddu_previous, a1);
			u_middle.AddIntoThis(u_previous);
			return u_middle;
		}

		private IGlobalVector R_hat_middle(IGlobalVector r_previous, IGlobalVector r_current, double p)
		{

			IGlobalVector r_hat_middle = algebraicModel.CreateZeroVector();
			r_hat_middle = r_previous.LinearCombination(1.0 - p, r_current, p);
			return r_hat_middle;
		}

		private IGlobalVector R_roundhat_middle(IGlobalVector r_hat_middle, IGlobalVector u_middle, IGlobalVector du_previous, IGlobalVector ddu_previous,
		   IGlobalMatrix kMatrix, IGlobalMatrix cMatrix, double a0)
		{
			IGlobalVector r_roundhat_middle = algebraicModel.CreateZeroVector();
			r_roundhat_middle = kMatrix.MultiplyVector()
			double[] stifPart = VectorOperations.MatrixVectorProduct(kMatrix, u_middle);
			double[] duplusa0ddu = VectorOperations.VectorVectorAddition(du_previous,
				VectorOperations.VectorScalarProductNew(ddu_previous, a0));
			double[] dampPart = VectorOperations.MatrixVectorProduct(cMatrix, duplusa0ddu);
			double[] r_roundhat1 = VectorOperations.VectorVectorSubtraction(r_hat_middle, stifPart);
			double[] r_roundhat2 = VectorOperations.VectorVectorSubtraction(r_roundhat1, dampPart);
			return r_roundhat2;
		}

		public void Solve()
		{
			throw new NotImplementedException();
		}
	}
}
