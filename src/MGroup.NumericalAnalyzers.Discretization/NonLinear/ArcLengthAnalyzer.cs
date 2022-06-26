using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.NumericalAnalyzers.NonLinear;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.NumericalAnalyzers.Discretization.NonLinear
{
	/// <summary>
	/// This class solves the nonlinear system of equations using the arc length method.
	/// </summary>
	public class ArcLengthAnalyzer : NonLinearAnalyzerBase
	{
		private int sign;
		private int previousIterations;
		private int increment;
		private int iteration;
		private IGlobalVector incrSolution;
		private IGlobalVector resSolution;
		private IGlobalVector curSolution;
		private IGlobalVector prevIncrSolution;
		private IGlobalVector polynomialSolution;
		private IGlobalVector externalRhs;
		private IGlobalVector internalRhs;
		private IGlobalVector rhsResidual;
		private IGlobalVector globalRhs; //TODO: This was originally readonly 
		private double globalRhsNormInitial; //TODO: This can probably be a local variable.
		private double incrementalRhsNormInitial;
		private double initialLoadFactor;
		private double lambda;
		private double dlambda;
		private double dlambda1;
		private double dlambda2;
		private double firstError;
		bool append = false;

		//polynomial coefficients
		private double a;
		private double b;
		private double c;

		//parameters
		public double shape = 0; // 0 for cylindrical, 1 for spherical, intermediate values (0-1) for elliptic
		public double initialdlambda; //load increment at first step
		public double deltaS; //constraint radius
		public double Id = 4; //desired number of iterations (3-5)
		public bool constS = true;

		/// <summary>
		/// This class solves the linearized geoemtrically nonlinear system of equations according to Newton-Raphson's load control incremental-iterative method.
		/// </summary>
		/// <param name="model">Instance of the model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param>
		/// <param name="numIncrements">Number of total load increments</param>
		/// <param name="maxIterationsPerIncrement">Number of maximum iterations within a load increment</param>
		/// <param name="numIterationsForMatrixRebuild">Number of iterations for the rebuild of the siffness matrix within a load increment</param>
		/// <param name="residualTolerance">Tolerance for the convergence criterion of the residual forces</param>
		private ArcLengthAnalyzer(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider,
			INonLinearModelUpdater modelUpdater,
			int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance)
			: base(algebraicModel, solver, provider, modelUpdater, numIncrements, maxIterationsPerIncrement,
				numIterationsForMatrixRebuild, residualTolerance)
		{ }


		/// <summary>
		/// Solves the nonlinear equations and calculates the displacements vector.
		/// </summary>
		public override void Solve()
		{
			InitializeLogs();

			DateTime start = DateTime.Now;
			UpdateInternalVectors();//TODOMaria this divides the externally applied load by the number of increments and scatters it to all subdomains and stores it in the class subdomain dictionary and total external load vector
			initialdlambda = 1d / numIncrements;
			dlambda = initialdlambda;
			lambda += dlambda;
			for (increment = 0; increment < numIncrements; increment++)
			{
				double errorNorm = 0;
				ClearIncrementalSolutionVector();//TODOMaria this sets du to 0

				for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
				{
					UpdateRhs(increment);//TODOMaria this copies the residuals stored in the class dictionary to the subdomains
					solver.Solve();
					UpdateSolutionsDictionary(incrSolution);

					UpdateRhsResidual(lambda);

					UpdateRhs();
					solver.Solve();
					UpdateSolutionsDictionary(resSolution);

					if (iteration == 0 && increment == 0)
					{
						deltaS = dlambda * Math.Sqrt((Math.Pow(shape, 2) * externalRhs.DotProduct(externalRhs)) + incrSolution.DotProduct(incrSolution));
						curSolution = incrSolution.Scale(dlambda);
						UpdateSolution(increment, iteration, curSolution);
					}
					else if (iteration == 0 && increment != 0)
					{
						if (constS == false)
						{
							deltaS = (double)deltaS * Id / previousIterations;
						}

						if (prevIncrSolution.DotProduct(incrSolution) > 0)
						{
							sign = +1;
						}
						else
						{
							sign = -1;
						}

						dlambda = sign * deltaS / Math.Sqrt((Math.Pow(shape, 2) * externalRhs.DotProduct(externalRhs)) + incrSolution.DotProduct(incrSolution));
						lambda += dlambda;
						curSolution = resSolution.Axpy(incrSolution, dlambda);
						UpdateSolution(increment, iteration, curSolution);
					}
					else if (iteration != 0)
					{
						polynomialSolution = du.Add(resSolution);
						a = incrSolution.DotProduct(incrSolution);
						b = 2 * incrSolution.DotProduct(polynomialSolution);
						c = polynomialSolution.DotProduct(polynomialSolution) - Math.Pow(deltaS, 2);
						var det = (b * b) - (4 * a * c);
						if (det > 0) //two real roots
						{
							var dlambda1 = (-b - Math.Sqrt(det)) / (2 * a);
							var curSolution1 = resSolution.Axpy(incrSolution, dlambda1);
							var cosine1 = du.DotProduct(uPlusdu.Add(curSolution1).Subtract(u)) / (du.Norm2() * uPlusdu.Add(curSolution1).Subtract(u).Norm2());

							var dlambda2 = (-b + Math.Sqrt(det)) / (2 * a);
							var curSolution2 = resSolution.Axpy(incrSolution, dlambda2);
							var cosine2 = du.DotProduct(uPlusdu.Add(curSolution2).Subtract(u)) / (du.Norm2() * uPlusdu.Add(curSolution2).Subtract(u).Norm2());

							if (cosine1 > cosine2)
							{
								dlambda = dlambda1;
								lambda += dlambda;
								curSolution = curSolution1.Copy();
							}
							else
							{
								dlambda = dlambda2;
								lambda += dlambda;
								curSolution = curSolution2.Copy();
							}
						}
						else //global optimum
						{
							dlambda = -b / (2 * a);
							lambda += dlambda;
							curSolution = resSolution.Axpy(incrSolution, dlambda);
						}
					}

					UpdateSolution(increment, iteration, curSolution);
					UpdateInternalRhs();
					double residualNormCurrent = UpdateResidualForcesAndNorm(lambda);
					errorNorm = incrementalRhsNormInitial != 0 ? residualNormCurrent / incrementalRhsNormInitial : 0;// (rhsNorm*increment/increments) : 0;//TODOMaria this calculates the internal force vector and subtracts it from the external one (calculates the residual)

					IncrementalDisplacementsLog?.StoreDisplacements(uPlusdu);
					TotalDisplacementsPerIterationLog?.StoreDisplacements(uPlusdu);

					if (iteration == 0)
					{
						firstError = errorNorm;
					}

					if (errorNorm < residualTolerance)
					{
						if (IncrementalLog != null)
						{
							IncrementalLog.LogTotalDataForIncrement(increment, iteration, errorNorm, uPlusdu, internalRhs);
						}
						break;
					}

					if ((iteration + 1) % numIterationsForMatrixRebuild == 0)
					{
						provider.Reset();
						BuildMatrices();
					}
				}

				Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
				previousIterations = iteration + 1;

				prevIncrSolution = du.Copy();
				SaveMaterialStateAndUpdateSolution();
			}

			DateTime end = DateTime.Now;
			StoreLogResults(start, end);
		}

		protected override void InitializeInternalVectors()
		{
			base.InitializeInternalVectors();
		}

		private void UpdateRhs()
		{
			solver.LinearSystem.RhsVector.CopyFrom(rhsResidual);
		}

		private void UpdateInternalRhs()
		{
			internalRhs.Clear();
			var internalRHS = modelUpdater.CalculateResponseIntegralVector(uPlusdu);
			provider.ProcessInternalRhs(uPlusdu, internalRHS);

			if (parentAnalyzer != null)
			{
				var otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(uPlusdu);
				internalRHS.AddIntoThis(otherRhsComponents);
			}

			internalRhs.Add(internalRHS);
		}

		private void UpdateSolution(int currentIncrement, int iteration, IGlobalVector solution)
		{
			if (currentIncrement == 0 && iteration == 0)
			{
				du.Clear();
				uPlusdu.Clear();
				du.AddIntoThis(solution);
				uPlusdu.AddIntoThis(solution);
				du.SubtractIntoThis(u);
			}
			else
			{
				du.AddIntoThis(solution);
				uPlusdu.Clear();
				uPlusdu.AddIntoThis(u);
				uPlusdu.AddIntoThis(du);
			}
		}

		private double UpdateResidualForcesAndNorm(double lambda)
		{
			globalRhs.Clear();
			solver.LinearSystem.RhsVector.Clear();
			solver.LinearSystem.RhsVector.AddIntoThis(externalRhs.Scale(lambda));
			solver.LinearSystem.RhsVector.SubtractIntoThis(internalRhs);
			rhsResidual = solver.LinearSystem.RhsVector.Copy();

			return provider.CalculateRhsNorm(globalRhs);
		}

		private void UpdateSolutionsDictionary(IGlobalVector solutionsDictionary)
		{
			solutionsDictionary.Clear();
			solutionsDictionary.Add(solver.LinearSystem.Solution.Copy());
		}

		private void UpdateRhsResidual(double lambda)
		{
			rhsResidual = externalRhs.Scale(lambda).Subtract(internalRhs);
		}

		public class Builder : NonLinearAnalyzerBuilderBase
		{
			public Builder(IModel model, IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider, int numIncrements)
				: base(model, algebraicModel, solver, provider, numIncrements)
			{
				MaxIterationsPerIncrement = 1000;
				NumIterationsForMatrixRebuild = 1;
				ResidualTolerance = 1E-3;
			}

			public ArcLengthAnalyzer Build() => new ArcLengthAnalyzer(algebraicModel, solver, provider, ModelUpdater, 
					numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance);
		}
	}
}
