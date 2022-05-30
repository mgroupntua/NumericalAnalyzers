using System;
using System.Collections.Generic;

using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers.Logging
{
	public class TotalDisplacementsPerIterationLog : IAnalysisWorkflowLog
    {
		private readonly IVectorValueExtractor resultsExtractor;
		private readonly List<Table<INode, IDofType, double>> dofDisplacementsPerIter;

        /// <summary>
        /// Initializes a new instance of <see cref="TotalDisplacementsPerIterationLog"/>.
        /// </summary>
        /// <param name="watchDofs">Which freedom degrees to track for each subdomain.</param>
        public TotalDisplacementsPerIterationLog(IList<(INode, IDofType)> watchDofs, IVectorValueExtractor resultsExtractor)
        {
            this.WatchDofs = watchDofs;
			this.resultsExtractor = resultsExtractor;
			this.dofDisplacementsPerIter = new List<Table<INode, IDofType, double>>();
        }

		public IList<(INode, IDofType)> WatchDofs { get; }

		/// <summary>
		/// Stores the total displacements = u_converged + du, for a new iteration.
		/// </summary>
		/// <param name="totalDisplacements">The total displacements for each subdomain.</param>
		public void StoreDisplacements(IGlobalVector totalDisplacements)
        {
            var currentIterDisplacements = new Table<INode, IDofType, double>();
			foreach ((INode node, IDofType dof) in WatchDofs)
			{
				currentIterDisplacements[node, dof] = resultsExtractor.ExtractSingleValue(totalDisplacements, node, dof);
			}
			dofDisplacementsPerIter.Add(currentIterDisplacements);
        }

        public double GetTotalDisplacement(int iteration, INode node, IDofType dof) 
            => dofDisplacementsPerIter[iteration][node, dof];

        public void StoreResults(DateTime startTime, DateTime endTime, IGlobalVector solution)
        {
            throw new NotImplementedException();
        }
    }
}
