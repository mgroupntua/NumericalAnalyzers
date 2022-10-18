using System.Linq;

using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Providers;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.NumericalAnalyzers.Discretization.NonLinear
{
	public class NonLinearModelUpdater : INonLinearModelUpdater
	{
		private readonly IAlgebraicModel algebraicModel;
		private readonly ElementInternalRhsProvider rhsProvider = new ElementInternalRhsProvider();


		public NonLinearModelUpdater(IAlgebraicModel algebraicModel)
		{
			this.algebraicModel = algebraicModel;
		}

		//TODO: I suggest splitting this into 2 methods. One for updating the elements/materials and one for calculating the internal rhs
		public IGlobalVector CalculateResponseIntegralVector(IGlobalVector solution)
		{
			var dirichletBoundaryConditions = algebraicModel.BoundaryConditionsInterpreter.GetDirichletBoundaryConditionsWithNumbering()
				.Select(x => new NodalBoundaryCondition(x.Value.Node, x.Key.DOF, x.Value.Amount));
			// First update the state of the elements
			algebraicModel.DoPerElement<IElementType>(element =>
			{
				double[] elementDisplacements = algebraicModel.ExtractElementVector(solution, element);
				element.MapNodalBoundaryConditionsToElementVector(dirichletBoundaryConditions, elementDisplacements);
				element.CalculateResponse(elementDisplacements);
			});

			// Then calculate the internal rhs vector
			IGlobalVector internalRhs = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(internalRhs, rhsProvider);
			return internalRhs;
		}

		public void UpdateState(IHaveState externalState)
		{
			algebraicModel.DoPerElement<IElementType>(element =>
			{
				element.SaveConstitutiveLawState(externalState);
			});
		}
	}
}
