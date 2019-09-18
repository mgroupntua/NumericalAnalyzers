using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.Analyzers.Interfaces
{
    public interface IKinematicRelationsStrategy
    {
        double[,] GetNodalKinematicRelationsMatrix(INode boundaryNode);
    }
}
