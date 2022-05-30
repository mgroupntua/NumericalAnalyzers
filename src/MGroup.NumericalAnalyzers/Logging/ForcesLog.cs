using System;
using System.Collections.Generic;
using System.Text;

using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers.Logging
{
	public class ForcesLog : IAnalysisWorkflowLog
    {
        private readonly IElementType[] elements;
		private readonly IVectorValueExtractor resultsExtractor;
		private readonly Dictionary<int, double[]> forces = new Dictionary<int, double[]>();

        public ForcesLog(IElementType[] elements, IVectorValueExtractor resultsExtractor)
        {
            this.elements = elements;
			this.resultsExtractor = resultsExtractor;
        }

        public Dictionary<int, double[]> Forces { get { return forces; } }
        public DateTime StartTime { get; set; }
        public DateTime EndTime { get; set; }

        public override string ToString()
        {
            StringBuilder s = new StringBuilder();
            foreach (int id in forces.Keys)
            {
                s.Append(String.Format("({0}): ", id));
                for (int i = 0; i < forces[id].Length; i++)
                    s.Append(String.Format("{0:0.00000}/", forces[id][i]));
                s.Append("; ");
            }
            return s.ToString();
        }

        #region IResultStorage Members

        public void StoreResults(DateTime startTime, DateTime endTime, IGlobalVector solutionVector)
        {
			StartTime = startTime;
			EndTime = endTime;
			foreach (IElementType e in elements)
			{
				double[] localVector = resultsExtractor.ExtractElementVector(solutionVector, e);
				forces[e.ID] = e.CalculateResponseIntegralForLogging(localVector);

				//for (int i = 0; i < stresses[e.ID].Length; i++)
				//    Debug.Write(stresses[e.ID][i]);
				//Debug.WriteLine("");
			}
		}

        #endregion
    }
}
