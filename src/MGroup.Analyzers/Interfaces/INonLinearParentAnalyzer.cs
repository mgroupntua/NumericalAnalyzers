using MGroup.LinearAlgebra.Vectors;
using MGroup.Solvers.LinearSystems;

namespace MGroup.Analyzers.Interfaces
{
    //TODO: Confusing name. The child analyzer of this is a nonlinear analyzer.
    public interface INonLinearParentAnalyzer : IParentAnalyzer
    {
        IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution);
    }
}
