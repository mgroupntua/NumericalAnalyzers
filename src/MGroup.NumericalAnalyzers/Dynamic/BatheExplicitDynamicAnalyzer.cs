using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution.AlgebraicModel;

namespace MGroup.NumericalAnalyzers.Dynamic
{
	public class BatheExplicitDynamicAnalyzer : INonLinearParentAnalyzer
	{
		public IChildAnalyzer ChildAnalyzer { get; }

		public IAnalysisWorkflowLog[] Logs => null;

		private List<double> Calculate_qValues(double p)
		{
			double q1 = (1.0 - 2.0 * p) / (2.0 * p * (1.0 - p));
			double q2 = 0.5 - p * q1;
			double q0 = -q1 - q2 + 0.5;
			List<double> q = new List<double>();
			q.Add(q0);
			q.Add(q1);
			q.Add(q2);
			return q;
		}

		private List<double> Calculate_aValues(double p, List<double> q, double tStep)
		{
			double q0 = q[0];
			double q1 = q[1];
			double q2 = q[2];
			double deltat = tStep;
			double a0 = p * deltat;
			double a1 = 0.5 * Math.Pow(p * deltat, 2);
			double a2 = a0 / 2.0;
			double a3 = (1.0 - p) * deltat;
			double a4 = 0.5 * Math.Pow((1.0 - p) * deltat, 2);
			double a5 = q0 * a3;
			double a6 = (0.5 + q1) * a3;
			double a7 = q2 * a3;
			List<double> a = new List<double>();
			a.Add(a0);
			a.Add(a1);
			a.Add(a2);
			a.Add(a3);
			a.Add(a4);
			a.Add(a5);
			a.Add(a6);
			a.Add(a7);
			return a;
		}
		public void BuildMatrices()
		{
			throw new NotImplementedException();
		}

		public IGlobalVector GetOtherRhsComponents(IGlobalVector currentSolution)
		{
			throw new NotImplementedException();
		}

		public void Initialize(bool isFirstAnalysis)
		{
			throw new NotImplementedException();
		}

		public void Solve()
		{
			throw new NotImplementedException();
		}
	}
}
