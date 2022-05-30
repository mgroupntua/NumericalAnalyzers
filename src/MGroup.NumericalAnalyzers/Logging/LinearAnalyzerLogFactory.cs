using System.Collections.Generic;

using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution.AlgebraicModel;

namespace MGroup.NumericalAnalyzers.Logging
{
	public class LinearAnalyzerLogFactory : ILogFactory
    {
		private readonly IList<(INode, IDofType)> dofs;
		private readonly IElementType[] stressElements, forceElements;
		private readonly IVectorValueExtractor resultsExtractor;

		public LinearAnalyzerLogFactory(IList<(INode, IDofType)> dofs, IElementType[] stressElements, IElementType[] dofElements,
			IVectorValueExtractor resultsExtractor)
        {
            this.dofs = dofs;
			this.stressElements = stressElements;
            this.forceElements = dofElements;
			this.resultsExtractor = resultsExtractor;
        }

        public LinearAnalyzerLogFactory(IList<(INode, IDofType)> dofs, IVectorValueExtractor resultsExtractor) 
			: this(dofs, new IElementType[0], new IElementType[0], resultsExtractor)
        {
        }

        public IAnalysisWorkflowLog[] CreateLogs()
        {
            var l = new List<IAnalysisWorkflowLog>();
            l.Add(new DOFSLog(dofs, resultsExtractor));
            if (stressElements.Length > 0)
			{
				l.Add(new StressesLog(stressElements, resultsExtractor));
			}
			if (forceElements.Length > 0)
			{
				l.Add(new ForcesLog(forceElements, resultsExtractor));
			}
			return l.ToArray();
        }
    }
}
