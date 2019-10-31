using System;
using System.Collections.Generic;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization;

namespace MGroup.NumericalAnalyzers.NonLinear
{
	/// <summary>
	/// Subdomain state update class that accounts for non zero initial conditions (displacements).
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public class NonLinearSubdomainUpdaterWithInitialConditions : INonLinearSubdomainUpdater
	{
		private readonly ISubdomain subdomain;

		public NonLinearSubdomainUpdaterWithInitialConditions(ISubdomain subdomain)
		{
			this.subdomain = subdomain;
		}

		public IVector GetRHSFromSolutionWithInitialDisplacemntsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{
			return this.subdomain.GetRHSFromSolutionWithInitialDisplacementsEffect(solution, dSolution, boundaryNodes,
				initialConvergedBoundaryDisplacements, totalBoundaryDisplacements,
				nIncrement, totalIncrements);
		}

		public void ResetState()
		{
			this.subdomain.ClearMaterialStresses();
		}

		public void UpdateState()
		{
			this.subdomain.SaveMaterialState();
		}

		public void ScaleConstraints(double scalingFactor)
		{
			throw new NotSupportedException();
		}

		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
		{
			throw new NotSupportedException();
			return this.subdomain.GetRhsFromSolution(solution, dSolution);
		}
	}
}
