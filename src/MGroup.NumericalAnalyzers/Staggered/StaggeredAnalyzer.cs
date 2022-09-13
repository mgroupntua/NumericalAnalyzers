using System;
using System.Diagnostics;
using System.Linq;

using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers.Staggered
{
	public class StaggeredAnalyzer : StaggeredAnalyzerBase
	{
		public StaggeredAnalyzer(IParentAnalyzer[] analyzers, ISolver[] solvers, CreateNewModelDelegate createNewModel, int maxStaggeredSteps, double tolerance)
		{
			this.maxStaggeredSteps = maxStaggeredSteps;
			this.tolerance = tolerance;
			this.analyzers = analyzers;
			this.solvers = solvers;
			this.CreateNewModel = createNewModel;

			currentSolutions = new IGlobalVector[analyzers.Length];
		}

		public override void Solve() => Solve(analyzers.Select(x => (Action)x.Solve).ToArray);
	}
}
