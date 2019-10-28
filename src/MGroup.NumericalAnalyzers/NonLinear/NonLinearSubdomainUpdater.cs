namespace MGroup.NumericalAnalyzers.NonLinear
{
	using MGroup.MSolve.AnalysisWorkflow;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.Interfaces;

	public class NonLinearSubdomainUpdater : INonLinearSubdomainUpdater
	{
		private readonly ISubdomain subdomain;

		public NonLinearSubdomainUpdater(ISubdomain subdomain)
		{
			this.subdomain = subdomain;
		}

		public void ScaleConstraints(double scalingFactor)
		{
			this.subdomain.ScaleConstraints(scalingFactor);
		}

		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
		{
			return subdomain.GetRhsFromSolution(solution, dSolution);
		}

		public void ResetState()
		{
			this.subdomain.ClearMaterialStresses();
		}

		public void UpdateState()
		{
			this.subdomain.SaveMaterialState();
		}
	}
}

