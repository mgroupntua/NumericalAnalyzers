using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.NumericalAnalyzers.NonLinear;
using MGroup.NumericalAnalyzers.Logging;

namespace MGroup.NumericalAnalyzers.Discretization.NonLinear
{
	/// <summary>
	/// This class solves the nonlinear system of equations using the arc length method.
	/// </summary>
	public class ArcLengthAnalyzer : NonLinearAnalyzerBase
	{
		private int sign;
		private int previousIterations;
		private IGlobalVector incrSolution;
		private IGlobalVector resSolution;
		private IGlobalVector curSolution;
		private IGlobalVector prevIncrSolution;
		private IGlobalVector polynomialSolution;
		private IGlobalVector rhsResidual;
		private IGlobalVector globalRhs;
		private double incrementalRhsNormInitial;
		private double lambda;
		private double dlambda;
		private double firstError;

		// polynomial coefficients
		private double a;
		private double b;
		private double c;

		// parameters
		private double shape = 0;
		private double initialdlambda;
		private double deltaS;
		private double numOfIterations = 4;
		private bool constConstraint = true;
		private bool stopIfNotConverged = false;

		/// <summary>
		/// This class solves the linearized geometrically nonlinear system of equations according to the Arc Length incremental-iterative method.
		/// </summary>
		/// <param name="algebraicModel">Instance of the algebraic model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param>
		/// <param name="numIncrements">Number of total load increments</param>
		/// <param name="maxIterationsPerIncrement">Number of maximum iterations within a load increment</param>
		/// <param name="numIterationsForMatrixRebuild">Number of iterations for the rebuild of the siffness matrix within a load increment</param>
		/// <param name="residualTolerance">Tolerance for the convergence criterion of the residual forces</param>
		/// <param name="shape">Option for the shape of the constraint - 0 for cylindrical, 1 for spherical, intermediate values for elliptic (default : shape = 0)</param>
		/// <param name="constConstraint">Option for constant radius of the constraint (default : constConstraint = 'true')</param>
		/// <param name="numOfIterations">(only usefull for constConstraint = false) Number of expected iterations within a load increment (default : numOfIterations = 4)</param>
		private ArcLengthAnalyzer(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider, int numIncrements, int maxIterationsPerIncrement, 
			int numIterationsForMatrixRebuild, double residualTolerance, double shape, bool constConstraint, int numOfIterations, bool stopIfNotConverged)
			: base(algebraicModel, solver, provider, numIncrements, maxIterationsPerIncrement,
				numIterationsForMatrixRebuild, residualTolerance)
		{
			this.shape = shape;
			this.constConstraint = constConstraint;
			this.numOfIterations = numOfIterations;
			this.stopIfNotConverged = stopIfNotConverged;
			this.analysisStatistics.AlgorithmName = "Arc length analyzer";
		}

		/// <summary>
		/// Solves the nonlinear equations and calculates the displacements vector.
		/// </summary>
		public override void Solve()
		{
			bool notConverged = false;

			analysisStatistics.NumIterationsRequired = 0;
			InitializeLogs();
			InitializeVectors();

			DateTime start = DateTime.Now;
			UpdateInternalVectors();
			initialdlambda = 1d;
			dlambda = initialdlambda;
			lambda += dlambda;
			for (int increment = 0; increment < numIncrements; increment++)
			{
				double errorNorm = 0;
				ClearIncrementalSolutionVector();

				int iteration = 0;
				for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
				{
					analysisStatistics.NumIterationsRequired++;

					if (iteration == maxIterationsPerIncrement - 1)
					{
						notConverged = true;
						analysisStatistics.ResidualNormRatioEstimation = errorNorm;
						if (stopIfNotConverged)
						{
							return;
						}
						else
						{
							break;
						}
					}

					if (double.IsNaN(errorNorm))
					{
						notConverged = true;
						analysisStatistics.ResidualNormRatioEstimation = errorNorm;
						if (stopIfNotConverged)
						{
							return;
						}
						else
						{
							break;
						}
					}

					UpdateRhs(increment);
					solver.Solve();
					UpdateSolution(incrSolution);

					UpdateRhsWithResidual();
					solver.Solve();
					UpdateSolution(resSolution);

					if (iteration == 0 && increment == 0)
					{
						deltaS = dlambda * Math.Sqrt((Math.Pow(shape, 2) * Math.Pow(lambda, 2) * rhsIncrement.DotProduct(rhsIncrement)) + incrSolution.DotProduct(incrSolution));
						curSolution = incrSolution.Scale(dlambda);
						UpdateSolution(increment, iteration, curSolution);
					}
					else if (iteration == 0 && increment != 0)
					{
						if (constConstraint == false)
						{
							deltaS = (double)deltaS * numOfIterations / previousIterations;
						}

						if (prevIncrSolution.DotProduct(incrSolution) > 0)
						{
							sign = +1;
						}
						else
						{
							sign = -1;
						}

						dlambda = sign * deltaS / Math.Sqrt((Math.Pow(shape, 2) * Math.Pow(lambda, 2) * rhsIncrement.DotProduct(rhsIncrement)) + incrSolution.DotProduct(incrSolution));
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
						UpdateSolution(increment, iteration, curSolution);
					}

					IGlobalVector internalRhsVector = CalculateInternalRhs(increment, iteration);
					double residualNormCurrent = UpdateResidualForcesAndNorm(internalRhsVector);
					errorNorm = incrementalRhsNormInitial != 0 ? residualNormCurrent / incrementalRhsNormInitial : 0;// (rhsNorm*increment/increments) : 0;//TODOMaria this calculates the internal force vector and subtracts it from the external one (calculates the residual)

					if (iteration == 0)
					{
						firstError = errorNorm;
					}

					if (IncrementalDisplacementsLog != null)
					{
						IncrementalDisplacementsLog.StoreDisplacements(uPlusdu);
					}

					if (errorNorm < residualTolerance)
					{
						if (analysisStatistics.ResidualNormRatioEstimation < errorNorm)
						{
							analysisStatistics.ResidualNormRatioEstimation = errorNorm;
						}

						if (IncrementalLog != null)
						{
							IncrementalLog.LogTotalDataForIncrement(increment, iteration, errorNorm, uPlusdu, internalRhsVector); //internalRhs
						}

						break;
					}

					if ((iteration + 1) % numIterationsForMatrixRebuild == 0)
					{
						provider.Reset();
						parentAnalyzer.BuildMatrices();
					}
				}

				if (TotalDisplacementsPerIterationLog != null)
				{
					TotalDisplacementsPerIterationLog.StoreDisplacements(uPlusdu);
				}

				Debug.WriteLine("AL {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
				previousIterations = iteration + 1;

				prevIncrSolution = du.Copy();
				SaveMaterialStateAndUpdateSolution();
			}

			analysisStatistics.HasConverged = !notConverged;
			DateTime end = DateTime.Now;
			StoreLogResults(start, end);
		}

		private void UpdateRhsWithResidual()
		{
			solver.LinearSystem.RhsVector.CopyFrom(rhsResidual);
		}

		new private IGlobalVector CalculateInternalRhs(int increment, int iteration)
		{
			IGlobalVector internalRhs = provider.CalculateResponseIntegralVector(uPlusdu);
			provider.ProcessInternalRhs(uPlusdu, internalRhs);

			if (parentAnalyzer != null)
			{
				IGlobalVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(uPlusdu);
				internalRhs.AddIntoThis(otherRhsComponents);
			}

			return internalRhs;
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

		private double UpdateResidualForcesAndNorm(IGlobalVector internalRhsVector)
		{
			globalRhs.Clear();
			solver.LinearSystem.RhsVector.Clear();
			solver.LinearSystem.RhsVector.AddIntoThis(rhsIncrement.Scale(lambda));
			solver.LinearSystem.RhsVector.SubtractIntoThis(internalRhsVector);
			rhsResidual = solver.LinearSystem.RhsVector.Copy();

			return provider.CalculateRhsNorm(solver.LinearSystem.RhsVector);
		}

		private void UpdateSolution(IGlobalVector solution)
		{
			solution.Clear();
			solution.CopyFrom(solver.LinearSystem.Solution);
		}

		private void InitializeVectors()
		{
			incrSolution = algebraicModel.CreateZeroVector();
			resSolution = algebraicModel.CreateZeroVector();
			curSolution = algebraicModel.CreateZeroVector();
			prevIncrSolution = algebraicModel.CreateZeroVector();
			polynomialSolution = algebraicModel.CreateZeroVector();
			rhsResidual = algebraicModel.CreateZeroVector();
			rhsResidual.CopyFrom(rhsIncrement);
			globalRhs = algebraicModel.CreateZeroVector();
			incrementalRhsNormInitial = provider.CalculateRhsNorm(rhsIncrement);
		}

		public class Builder : NonLinearAnalyzerBuilderBase
		{
			private double shape;
			private bool constConstraint;
			private int numOfIterations;
			private bool stopIfNotConverged = true;

			public Builder(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider, int numIncrements, double shape = 0, int numOfIterations = 4, 
				bool constConstraint = true, bool stopIfNotConverged = true)
				: base(algebraicModel, solver, provider, numIncrements)
			{
				MaxIterationsPerIncrement = 1000;
				NumIterationsForMatrixRebuild = 1;
				ResidualTolerance = 1E-3;
				this.shape = shape;
				this.constConstraint = constConstraint;
				this.numOfIterations = numOfIterations;
				this.stopIfNotConverged = stopIfNotConverged;
			}

			public ArcLengthAnalyzer Build() => new ArcLengthAnalyzer(algebraicModel, solver, provider, 
					numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance, shape, constConstraint, numOfIterations, stopIfNotConverged);
		}
	}
}
