using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;

//TODO: finding the contributing elements and the corresponding local dof indices can be done only once in the constructor.
namespace MGroup.NumericalAnalyzers.Logging
{
	/// <summary>
	/// This does not work if the requested node belongs to an element that contains embedded elements.
	/// </summary>
	internal class ConstrainedDofForcesCalculator
    {
		private readonly IVectorValueExtractor resultsExtractor;
		private readonly IEnumerable<INodalBoundaryCondition> boundaryConditions;

		internal ConstrainedDofForcesCalculator(IEnumerable<INodalBoundaryCondition> boundaryConditions, IVectorValueExtractor resultsExtractor)
        {
			this.boundaryConditions = boundaryConditions;
			this.resultsExtractor = resultsExtractor;
		}

        internal double CalculateForceAt(INode node, IDofType dofType, IGlobalVector totalDisplacements)
        {
            double totalForce = 0.0;

            foreach (IElementType element in node.ElementsDictionary.Values)
            {
                // It is possible that one of the elements at this node does not engage this dof type, in which case -1 will be returned.
                // We will not have any contribution from them. If none of the elements engage this dof type, the total force will always be 0.
                int monitorDofIdx = FindLocalDofIndex(element, node, dofType);
                if (monitorDofIdx == -1) continue;

                //TODO: if an element has embedded elements, then we must also take into account their forces.
                double[] totalElementDisplacements = CalculateElementDisplacements(element, boundaryConditions, totalDisplacements);
                double[] elementForces = element.CalculateResponseIntegralForLogging(totalElementDisplacements);

                totalForce += elementForces[monitorDofIdx];
            }

            return totalForce;
        }

		private double[] CalculateElementDisplacements(IElementType element, IEnumerable<INodalBoundaryCondition> boundaryConditions, IGlobalVector displacements)
		{
			double[] elementNodalDisplacements = resultsExtractor.ExtractElementVector(displacements, element);
			//DirichletElementLoad.ApplyBoundaryConditions(element, elementNodalDisplacements);
			element.MapNodalBoundaryConditionsToElementVector(boundaryConditions, elementNodalDisplacements);
			return elementNodalDisplacements;
		}

		/// <summary>
		/// Returns -1 if the element does not engage the requested <see cref="IDofType"/>
		/// </summary>
		private int FindLocalDofIndex(IElementType element, INode node, IDofType dofType)
        {
            int localNodeIdx = element.Nodes.Select(x=>x).ToList().IndexOf(node);
            Debug.Assert(localNodeIdx != -1, "The element does not contain this node.");
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            int localDofIdx = elementDofs[localNodeIdx].FindFirstIndex(dofType);
            int multNum = elementDofs[localNodeIdx].Count;
            int dofIdx = multNum * (localNodeIdx + 1) - (localDofIdx + 1);
            return dofIdx;
        }
    }
}
