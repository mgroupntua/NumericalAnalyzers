//TODO:  delete the original one (IEquivalentLoadsAssembler) in FEM.Interfaces

using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.Analyzers.Interfaces
{
    public interface IDirichletEquivalentLoadsAssembler
    {
        IVector GetEquivalentNodalLoads(ISubdomain subdomain, IVectorView solution, double constraintScalingFactor);
    }
}
