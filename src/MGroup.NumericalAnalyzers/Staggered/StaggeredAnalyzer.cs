using System;
using System.Collections.Generic;
using System.Linq;

using MGroup.LinearAlgebra.Iterative;
using MGroup.MSolve.AnalysisWorkflow;
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

			if (analyzers == null)
			{
				throw new ArgumentException("Analyzers is null");
			}

			if (solvers == null)
			{
				throw new ArgumentException("Solvers is null");
			}

			if (solvers.Length != analyzers.Length) 
			{
				throw new ArgumentException($"Solvers and analyzer length mismatch (analyzers: {analyzers.Length}, solvers: {solvers.Length})");
			}

			this.currentSolutions = new IGlobalVector[analyzers.Length];
			this.CreateNewModel = createNewModel;
		}

		public override void Solve() => Solve(analyzers.Select(x => (Action)x.Solve).ToArray);
	}
}
