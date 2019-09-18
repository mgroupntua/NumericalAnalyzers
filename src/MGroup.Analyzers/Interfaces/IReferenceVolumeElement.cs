using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.Analyzers.Interfaces
{
    public interface IReferenceVolumeElement
    {
        void ApplyBoundaryConditions();
        IMatrixView CalculateKinematicRelationsMatrix(ISubdomain subdomain);
        double CalculateRveVolume();
    }
}
