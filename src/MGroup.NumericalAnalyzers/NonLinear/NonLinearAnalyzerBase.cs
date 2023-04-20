using System;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.DataStructures;
using MGroup.LinearAlgebra.Iterative;

namespace MGroup.NumericalAnalyzers.NonLinear
{

	/// <summary>
	/// This class represents the base class of all nonlinear anaylsers and contains the basic information necessary for other analyzers
	/// </summary>
	public abstract class NonLinearAnalyzerBase : IChildAnalyzer
	{
		private const string CURRENTSOLUTION = "Current solution";
		private const string LASTRHS = "Last RHS";

		protected readonly int maxIterationsPerIncrement;
		protected readonly IAlgebraicModel algebraicModel;
		protected readonly int numIncrements;
		protected readonly int numIterationsForMatrixRebuild;
		protected readonly INonLinearProvider provider;
		protected readonly double residualTolerance;
		protected readonly ISolver solver;
		protected IGlobalVector rhsIncrement;
		protected IGlobalVector lastRhs;
		protected IGlobalVector u;
		protected IGlobalVector du;
		protected IGlobalVector uPlusdu;
		protected double globalRhsNormInitial;
		protected INonLinearParentAnalyzer parentAnalyzer = null;
		private GenericAnalyzerState currentState;
		protected IterativeStatistics analysisStatistics = new IterativeStatistics()
		{
			AlgorithmName = "Non-linear analyzer",
		};

		public NonLinearAnalyzerBase(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider,
			int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance)
		{
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.numIncrements = numIncrements;
			this.maxIterationsPerIncrement = maxIterationsPerIncrement;
			this.numIterationsForMatrixRebuild = numIterationsForMatrixRebuild;
			this.residualTolerance = residualTolerance;
		}

		public LinearAnalyzerLogFactory LogFactory { get; set; }

		public IGlobalVector CurrentAnalysisResult { get => u; }

		public IAnalysisWorkflowLog[] Logs { get; set; } = new IAnalysisWorkflowLog[0];

		public TotalDisplacementsPerIterationLog TotalDisplacementsPerIterationLog { get; set; }

		public IncrementalDisplacementsLog IncrementalDisplacementsLog { get; set; }

		public TotalLoadsDisplacementsPerIncrementLog IncrementalLog { get; set; }

		public IParentAnalyzer ParentAnalyzer
		{
			get => parentAnalyzer;
			set => parentAnalyzer = (INonLinearParentAnalyzer)value;
		}

		public IGlobalVector Responses { get; set; }

		public IGlobalVector CurrentAnalysisLinearSystemRhs { get => solver.LinearSystem.RhsVector; }

		public IterativeStatistics AnalysisStatistics => analysisStatistics;
		
		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;
				currentState.StateVectors[LASTRHS].CheckForCompatibility = false;

				u.CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);
				lastRhs.CopyFrom(currentState.StateVectors[LASTRHS]);

				currentState.StateVectors[LASTRHS].CheckForCompatibility = true;
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
			}
		}

		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this, new[]
			{
				(CURRENTSOLUTION, u),
				(LASTRHS, lastRhs),
			});

			return currentState;
		}

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();

		/// <summary>
		/// Initializes internal vector before the first analysis.
		/// </summary>
		public void Initialize(bool isFirstAnalysis)
		{
			InitializeInternalVectors(isFirstAnalysis);
		}

		protected IGlobalVector CalculateInternalRhs(int currentIncrement, int iteration)
		{
			if (currentIncrement == 0 && iteration == 0)
			{
				du.Clear();
				uPlusdu.Clear();
				du.AddIntoThis(solver.LinearSystem.Solution);
				uPlusdu.AddIntoThis(solver.LinearSystem.Solution);
				du.SubtractIntoThis(u);
			}
			else
			{
				du.AddIntoThis(solver.LinearSystem.Solution);
				uPlusdu.Clear();
				uPlusdu.AddIntoThis(u);
				uPlusdu.AddIntoThis(du);
			}

			IGlobalVector internalRhs = provider.CalculateResponseIntegralVector(uPlusdu);
			provider.ProcessInternalRhs(uPlusdu, internalRhs);

			if (parentAnalyzer != null)
			{
				IGlobalVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(uPlusdu);
				internalRhs.AddIntoThis(otherRhsComponents);
			}

			return internalRhs;
		}

		protected double UpdateResidualForcesAndNorm(int currentIncrement, int iteration, IGlobalVector internalRhs)
		{
			if (iteration == 0)
			{
				lastRhs.AddIntoThis(solver.LinearSystem.RhsVector);
			}

			solver.LinearSystem.RhsVector.CopyFrom(lastRhs);
			solver.LinearSystem.RhsVector.SubtractIntoThis(internalRhs);
			return provider.CalculateRhsNorm(solver.LinearSystem.RhsVector);
		}

		protected void ClearIncrementalSolutionVector()
		{
			du.Clear();
		}

		protected virtual void InitializeInternalVectors(bool isFirstAnalysis)
		{
			if (lastRhs == null)
			{
				lastRhs = algebraicModel.CreateZeroVector();
			}
			else
			{
				if (isFirstAnalysis)
				{
					lastRhs.Clear();
				}
			}

			if (rhsIncrement == null)
			{
				rhsIncrement = algebraicModel.CreateZeroVector();
			}
			else
			{
				if (isFirstAnalysis)
				{
					rhsIncrement.Clear();
				}
			}

			if (u == null)
			{
				u = algebraicModel.CreateZeroVector();
			}
			else
			{
				if (isFirstAnalysis)
				{
					u.Clear();
				}
			}

			if (du == null)
			{
				du = algebraicModel.CreateZeroVector();
			}
			else
			{
				du.Clear();
			}

			if (uPlusdu == null)
			{
				uPlusdu = algebraicModel.CreateZeroVector();
			}
			else
			{
				uPlusdu.Clear();
			}
		}

		protected void InitializeLogs()
		{
			if (LogFactory != null)
			{
				Logs = LogFactory.CreateLogs();
			}

			if (IncrementalLog != null)
			{
				IncrementalLog.Initialize();
			}
		}

		protected void SaveMaterialStateAndUpdateSolution()
		{
			ParentAnalyzer.CreateState();
			provider.UpdateState(ParentAnalyzer.CurrentState);
			u.AddIntoThis(du);
		}

		protected void StoreLogResults(DateTime start, DateTime end)
		{
			foreach (var l in Logs)
			{
				l.StoreResults(start, end, u);
			}
		}

		protected void UpdateInternalVectors()
		{
			rhsIncrement = solver.LinearSystem.RhsVector.Subtract(lastRhs);
			rhsIncrement.ScaleIntoThis(1 / (double)numIncrements);
			globalRhsNormInitial = provider.CalculateRhsNorm(solver.LinearSystem.RhsVector);
		}

		protected virtual void UpdateRhs(int step)
		{
			solver.LinearSystem.RhsVector.CopyFrom(rhsIncrement);
		}

		/// <summary>
		/// This class solves system and calculates the displacements vector.
		/// </summary>
		public abstract void Solve();
	}
}
