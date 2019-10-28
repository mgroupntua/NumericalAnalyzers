using System.Collections.Generic;
using MGroup.MSolve.Discretization.Interfaces;
using MGroup.MSolve.Logging.Interfaces;

namespace MGroup.NumericalAnalyzers.Logging
{
    public class LinearAnalyzerLogFactory : ILogFactory
    {
        private readonly int[] dofs;
        private readonly IElement[] stressElements, forceElements;

        public LinearAnalyzerLogFactory(int[] dofs, IElement[] stressElements, IElement[] dofElements)
        {
            this.dofs = dofs;
            this.stressElements = stressElements;
            this.forceElements = dofElements;
        }

        public LinearAnalyzerLogFactory(int[] dofs) : this(dofs, new IElement[0], new IElement[0])
        {
        }

        public IAnalyzerLog[] CreateLogs()
        {
            var l = new List<IAnalyzerLog>();
            l.Add(new DOFSLog(dofs));
            if (stressElements.Length > 0)
                l.Add(new StressesLog(stressElements));
            if (forceElements.Length > 0)
                l.Add(new ForcesLog(forceElements));
            return l.ToArray();
        }
    }
}
