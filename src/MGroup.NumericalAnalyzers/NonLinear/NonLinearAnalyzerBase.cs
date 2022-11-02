using System;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.DataStructures;

namespace MGroup.NumericalAnalyzers.NonLinear
{

	/// <summary>
	/// This class represents the base class of all nonlinear anaylsers and contains the basic information necessary for other analyzers
	/// </summary>
	public abstract class NonLinearAnalyzerBase : IChildAnalyzer
	{
		private const string CURRENTSOLUTION = "Current solution";

		protected readonly int maxIterationsPerIncrement;
		protected readonly IAlgebraicModel algebraicModel;
		protected readonly int numIncrements;
		protected readonly int numIterationsForMatrixRebuild;
		protected readonly INonLinearProvider provider;
		protected readonly double residualTolerance;
		protected readonly ISolver solver;
		protected readonly INonLinearModelUpdater modelUpdater;
		protected IGlobalVector rhs;
		protected IGlobalVector u;
		protected IGlobalVector du;
		protected IGlobalVector uPlusdu;
		protected double globalRhsNormInitial;
		protected INonLinearParentAnalyzer parentAnalyzer = null;
		private GenericAnalyzerState currentState;

		public NonLinearAnalyzerBase(IAlgebraicModel algebraicModel, ISolver solver, INonLinearProvider provider,
			INonLinearModelUpdater modelUpdater,
			int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance)
		{
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.provider = provider;
			this.modelUpdater = modelUpdater;
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

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;

				u.CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);

				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
			}
		}

		//public IGlobalVector CurrentAnalysisLinearSystemRhs => throw new NotImplementedException();

		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this, new[]
			{
				(CURRENTSOLUTION, u),
			});

			return currentState;
		}

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();

		/// <summary>
		/// Builds the tangent stiffness matrix of the system.
		/// </summary>
		public void BuildMatrices()
		{
			if (parentAnalyzer == null)
			{
				throw new InvalidOperationException(
				"This Newton-Raphson nonlinear analyzer has no parent.");
			}

			parentAnalyzer.BuildMatrices();
		}

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

			IGlobalVector internalRhs = modelUpdater.CalculateResponseIntegralVector(uPlusdu);
			provider.ProcessInternalRhs(uPlusdu, internalRhs);

			if (parentAnalyzer != null)
			{
				IGlobalVector otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(uPlusdu);
				internalRhs.AddIntoThis(otherRhsComponents);
			}

			return internalRhs;
		}

		protected double UpdateResidualForcesAndNorm(int currentIncrement, IGlobalVector internalRhs)
		{
			solver.LinearSystem.RhsVector.Clear();
			for (int j = 0; j <= currentIncrement; j++)
			{
				solver.LinearSystem.RhsVector.AddIntoThis(rhs);
			}

			solver.LinearSystem.RhsVector.SubtractIntoThis(internalRhs);
			return provider.CalculateRhsNorm(solver.LinearSystem.RhsVector);
		}

		protected void ClearIncrementalSolutionVector()
		{
			du.Clear();
		}

		protected virtual void InitializeInternalVectors(bool isFirstAnalysis)
		{
			rhs = solver.LinearSystem.RhsVector.Copy();
			rhs.ScaleIntoThis(1 / (double)numIncrements);
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

			//u = algebraicModel.CreateZeroVector();
			//du = algebraicModel.CreateZeroVector();
			//uPlusdu = algebraicModel.CreateZeroVector();
			globalRhsNormInitial = provider.CalculateRhsNorm(solver.LinearSystem.RhsVector);
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
			modelUpdater.UpdateState(ParentAnalyzer.CurrentState);
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
			rhs = solver.LinearSystem.RhsVector.Copy();
			rhs.ScaleIntoThis(1 / (double)numIncrements);
			globalRhsNormInitial = provider.CalculateRhsNorm(solver.LinearSystem.RhsVector);
		}

		protected void UpdateRhs(int step)
		{
			solver.LinearSystem.RhsVector.CopyFrom(rhs);
		}

		/// <summary>
		/// This class solves system and calculates the displacements vector.
		/// </summary>
		public abstract void Solve();
	}
}
