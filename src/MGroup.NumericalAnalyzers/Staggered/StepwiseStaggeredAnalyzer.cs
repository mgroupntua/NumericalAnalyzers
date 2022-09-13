using System;
using System.Diagnostics;
using System.Linq;

using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers.Staggered
{
	public class StepwiseStaggeredAnalyzer : StaggeredAnalyzerBase
	{
		public StepwiseStaggeredAnalyzer(IParentAnalyzer[] analyzers, ISolver[] solvers, CreateNewModelDelegate createNewModel, int maxStaggeredSteps, double tolerance)
		{
			this.maxStaggeredSteps = maxStaggeredSteps;
			this.tolerance = tolerance;
			this.analyzers = analyzers;
			this.solvers = solvers;
			this.CreateNewModel = createNewModel;

			currentSolutions = new IGlobalVector[analyzers.Length];
		}

		public void SolveCurrentStep() => base.Solve(analyzers.Select(x => x is IStepwiseAnalyzer ? (Action)((IStepwiseAnalyzer)x).Solve : (Action)x.Solve).ToArray);

		public override void Solve()
		{
			int steps = analyzers.Max(x => x is IStepwiseAnalyzer ? ((IStepwiseAnalyzer)x).Steps : 0);
			for (int i = 0; i < steps; i++)
			{
				SolveCurrentStep();

				for (int j = 0; j < analyzers.Length; j++)
				{
					if (analyzers[j] is IStepwiseAnalyzer)
					{
						var analyzer = (IStepwiseAnalyzer)analyzers[j];
						analyzer.AdvanceStep();
					}
				}
			}
		}
	}
}
