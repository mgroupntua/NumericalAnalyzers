using System;
using System.Collections.Generic;
using System.Text;

using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers.Logging
{
	public class DOFSLog : IAnalysisWorkflowLog
    {
		private readonly IList<(INode, IDofType)> dofs;
		private readonly IVectorValueExtractor resultsExtractor;

        public DOFSLog(IList<(INode, IDofType)> dofs, IVectorValueExtractor resultsExtractor)
        {
            this.dofs = dofs;
			this.resultsExtractor = resultsExtractor;
		}

		public Table<INode, IDofType, double> DOFValues { get; } = new Table<INode, IDofType, double>();
        public DateTime StartTime { get; set; }
        public DateTime EndTime { get; set; }

        public override string ToString()
        {
            StringBuilder s = new StringBuilder();
            foreach ((INode node, IDofType dof, double val) in DOFValues)
			{
				s.Append($"(node {node.ID}, dof {dof}): {val}; ");
			}
			return s.ToString();
        }

        #region IResultStorage Members

        public void StoreResults(DateTime startTime, DateTime endTime, IGlobalVector solution)
        {
            StartTime = startTime;
            EndTime = endTime;
            foreach ((INode node, IDofType dof) in dofs)
            {
				DOFValues[node, dof] = resultsExtractor.ExtractSingleValue(solution, node, dof);
            }
        }

        #endregion
    }
}
