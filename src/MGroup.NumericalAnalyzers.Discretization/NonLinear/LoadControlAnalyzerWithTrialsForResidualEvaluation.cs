using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.NonLinear;
using MGroup.LinearAlgebra;

namespace MGroup.NumericalAnalyzers.Discretization.NonLinear
{
	public class LoadControlAnalyzerWithTrialsForResidualEvaluation : NonLinearAnalyzerBase
	{
		/// <summary>
		/// This class solves the linearized geoemtrically nonlinear system of equations according to Newton-Raphson's load control incremental-iterative method.
		/// </summary>
		/// <param name="model">Instance of the model that will be solved</param>
		/// <param name="solver">Instance of the solver that will solve the linear system of equations</param>
		/// <param name="provider">Instance of the problem type to be solved</param>
		/// <param name="subdomainUpdaters">Instance that updates constraints, right-hand-side vector, updates and resets state</param>
		/// <param name="numIncrements">Number of total load increments</param>
		/// <param name="maxIterationsPerIncrement">Number of maximum iterations within a load increment</param>
		/// <param name="numIterationsForMatrixRebuild">Number of iterations for the rebuild of the siffness matrix within a load increment</param>
		/// <param name="residualTolerance">Tolerance for the convergence criterion of the residual forces</param>
		private LoadControlAnalyzerWithTrialsForResidualEvaluation(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider,
			int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance)
			: base(algebraicModel, solver, provider, numIncrements, maxIterationsPerIncrement,
				numIterationsForMatrixRebuild, residualTolerance)
		{ }


		/// <summary>
		/// Solves the nonlinear equations and calculates the displacements vector.
		/// </summary>
		public override void Solve()
		{
			InitializeLogs();

			DateTime start = DateTime.Now;
			UpdateInternalVectors();
			for (int increment = 0; increment < numIncrements; increment++)
			{
				AnalysisState.newmarkIncrementNumber = increment;
				double errorNorm = 0;
				ClearIncrementalSolutionVector();
				UpdateRhs(increment);

				double firstError = 0;
				int iteration = 0;
				for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
				{
					AnalysisState.loadControlIteration = iteration;
					if (iteration == maxIterationsPerIncrement - 1)
					{
						return;
					}

					if (double.IsNaN(errorNorm))
					{
						return;
					}

					solver.Solve();
					StoreSystem();
					int nTrials = 4; int ratio = 4; double[] trialErrors = new double[nTrials]; double trialSize = (double)1 / ratio;
					for (int i1 = 0; i1 < nTrials; i1++)
					{
						var trialInternalRhsVector = CalculateInternalTrialRHS(increment, iteration, (i1 + 1) * trialSize);
						trialErrors[i1] = globalRhsNormInitial != 0 ? UpdateResidualForcesAndNorm(increment, trialInternalRhsVector) / globalRhsNormInitial : 0;
						RestoreSystem();
					}
					int minErrorTrialNo = 0;
					double minErrorTrialValue = trialErrors[0];
					for (int i1 = 1; i1 < nTrials; i1++)
					{
						if (trialErrors[i1] < minErrorTrialValue)
						{
							minErrorTrialNo = i1;
							minErrorTrialValue = trialErrors[i1];
						}
					}

					Console.WriteLine($"chosen trial no is {minErrorTrialNo + 1}"); Debug.WriteLine($"chosen trial no is {minErrorTrialNo + 1}");

					var internalRhsVector = CalculateInternalTrialRHS(increment, iteration, (minErrorTrialNo + 1) * trialSize);
					double residualNormCurrent = UpdateResidualForcesAndNorm(increment, internalRhsVector);
					errorNorm = globalRhsNormInitial != 0 ? residualNormCurrent / globalRhsNormInitial : 0;

					if (iteration == 0)
					{
						firstError = errorNorm;
					}

					if (IncrementalDisplacementsLog != null)
					{
						IncrementalDisplacementsLog.StoreDisplacements(uPlusdu);
					}

					if (TotalDisplacementsPerIterationLog != null)
					{
						TotalDisplacementsPerIterationLog.StoreDisplacements(uPlusdu);
					}

					if (errorNorm < residualTolerance)
					{
						if (IncrementalLog != null)
						{
							IncrementalLog.LogTotalDataForIncrement(increment, iteration, errorNorm, uPlusdu, internalRhsVector);
						}
						break;
					}

					if ((iteration + 1) % numIterationsForMatrixRebuild == 0)
					{
						provider.Reset();
						parentAnalyzer.BuildMatrices();
					}
				}
				Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
				SaveMaterialStateAndUpdateSolution();
			}
			DateTime end = DateTime.Now;
			StoreLogResults(start, end);
		}

		public class Builder : NonLinearAnalyzerBuilderBase
		{
			public Builder(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider, int numIncrements)
				: base(algebraicModel, solver, provider, numIncrements)
			{
				MaxIterationsPerIncrement = 1000;
				NumIterationsForMatrixRebuild = 1;
				ResidualTolerance = 1E-3;
			}

			public LoadControlAnalyzerWithTrialsForResidualEvaluation Build() => new LoadControlAnalyzerWithTrialsForResidualEvaluation(algebraicModel, solver, provider,
				numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance);
		}
	}
}
