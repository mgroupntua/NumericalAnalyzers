//TODO:  delete the original one (IEquivalentLoadsAssembler) in FEM.Interfaces

namespace MGroup.Analyzers.Interfaces
{
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.MSolve.Discretization.Interfaces;

	public interface IDirichletEquivalentLoadsAssembler
	{
		IVector GetEquivalentNodalLoads(ISubdomain subdomain, IVectorView solution, double constraintScalingFactor);
	}
}
