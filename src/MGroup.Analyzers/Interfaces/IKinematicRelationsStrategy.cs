namespace MGroup.Analyzers.Interfaces
{
	using MGroup.MSolve.Discretization.Interfaces;

	public interface IKinematicRelationsStrategy
	{
		double[,] GetNodalKinematicRelationsMatrix(INode boundaryNode);
	}
}
