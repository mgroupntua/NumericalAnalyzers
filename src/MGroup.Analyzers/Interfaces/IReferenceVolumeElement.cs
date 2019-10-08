namespace MGroup.Analyzers.Interfaces
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.Discretization.Interfaces;

	public interface IReferenceVolumeElement
	{
		void ApplyBoundaryConditions();

		IMatrixView CalculateKinematicRelationsMatrix(ISubdomain subdomain);

		double CalculateRveVolume();
	}
}
