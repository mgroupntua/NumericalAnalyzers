using System;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.AlgebraicModel;

namespace MGroup.NumericalAnalyzers.Discretization.NonLinear
{
	public abstract class NonLinearAnalyzerBuilderBase
	{
		protected int maxIterationsPerIncrement = 1000;
		//protected readonly IModel model;
		protected readonly IAlgebraicModel algebraicModel;
		protected readonly int numIncrements;
		protected int numIterationsForMatrixRebuild = 1;
		protected readonly INonLinearProvider provider;
		protected double residualTolerance = 1e-8;
		protected readonly ISolver solver;

		protected NonLinearAnalyzerBuilderBase(IAlgebraicModel algebraicModel, ISolver solver, 
			INonLinearProvider provider, int numIncrements)
		{
			//TODO: this should belong to all (child) analyzer builders
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.numIncrements = numIncrements;
		}

		public int MaxIterationsPerIncrement
		{
			get => maxIterationsPerIncrement;
			set
			{
				if (value < 1)
				{
					throw new ArgumentException("Max iterations per increment must be >= 1");
				}

				maxIterationsPerIncrement = value;
			}
		}

		public int NumIterationsForMatrixRebuild
		{
			get => numIterationsForMatrixRebuild;
			set
			{
				if (value < 1)
				{
					throw new ArgumentException("Iterations number for matrix rebuild must be >= 1");
				}

				numIterationsForMatrixRebuild = value;
			}
		}

		public double ResidualTolerance
		{
			get => residualTolerance;
			set
			{
				if (value <= 0.0)
				{
					throw new ArgumentException("Residual tolerance must be positive");
				}

				residualTolerance = value;
			}
		}
	}
}
