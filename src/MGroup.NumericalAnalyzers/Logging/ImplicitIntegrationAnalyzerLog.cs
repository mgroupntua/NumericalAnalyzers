using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.AnalysisWorkflow.Logging;

namespace MGroup.NumericalAnalyzers.Logging
{
    public class ImplicitIntegrationAnalyzerLog
    {
        private List<IAnalysisWorkflowLog> logs = new List<IAnalysisWorkflowLog>();
        private readonly HashSet<Type> logTypes = new HashSet<Type>();

        public List<IAnalysisWorkflowLog> Logs => logs;

        public void StoreResults(DateTime startTime, DateTime endTime, IAnalysisWorkflowLog log)
        {
            if (logTypes.Contains(log.GetType()) == false)
                logTypes.Add(log.GetType());
            logs.Add(log);
        }

        public void ClearResults()
        {
            logs.Clear();
        }

        public void GroupLogsByKind()
        {
            int interleaving = logTypes.Count;
            int items = logs.Count / interleaving;
            var deInterlacedLogs = new IAnalysisWorkflowLog[logs.Count];
            for (int i = 0; i < logs.Count; i++)
            {
                int item = i / interleaving;
                int group = i % interleaving;
                deInterlacedLogs[group * items + item] = logs[i];
            }
            logs = deInterlacedLogs.ToList<IAnalysisWorkflowLog>();
        }
    }
}
