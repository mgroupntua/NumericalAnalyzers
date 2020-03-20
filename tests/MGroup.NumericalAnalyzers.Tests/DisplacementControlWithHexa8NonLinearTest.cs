namespace MGroup.NumericalAnalyzers.Tests
{
	using System.Collections.Generic;

	using MGroup.Constitutive.Structural;
	using MGroup.Constitutive.Structural.ContinuumElements;
	using MGroup.FEM.Entities;
	using MGroup.FEM.Structural.Elements;
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.Integration.Quadratures;
	using MGroup.NumericalAnalyzers;
	using MGroup.NumericalAnalyzers.Logging;
	using MGroup.NumericalAnalyzers.NonLinear;
	using MGroup.Solvers.Direct;
	using Xunit;

	public class DisplacementControlWithHexa8NonLinearTest
	{
		private const int subdomainID = 0;

		[Fact]
		public void DisplacementControlWithHexa8NonLinear()
		{
			const int increments = 10;
			const double nodalDisplacement = -5.0;
			const double youngModulus = 4.0;
			const double poissonRatio = 0.4;

			// Create Model
			Model model = new Model();

			// Create Subdomain
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			// Create Elastic Material
			var solidMaterial = new ElasticMaterial3D()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			// Node creation
			Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			Node node2 = new Node(id: 2, x: 10.0, y: 0.0, z: 0.0);
			Node node3 = new Node(id: 3, x: 0.0, y: 10.0, z: 0.0);
			Node node4 = new Node(id: 4, x: 10.0, y: 10.0, z: 0.0);
			Node node5 = new Node(id: 5, x: 0.0, y: 0.0, z: 10.0);
			Node node6 = new Node(id: 6, x: 10.0, y: 0.0, z: 10.0);
			Node node7 = new Node(id: 7, x: 0.0, y: 10.0, z: 10.0);
			Node node8 = new Node(id: 8, x: 10.0, y: 10.0, z: 10.0);

			// Create List of nodes
			IList<Node> nodes = new List<Node>();
			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);
			nodes.Add(node4);
			nodes.Add(node5);
			nodes.Add(node6);
			nodes.Add(node7);
			nodes.Add(node8);

			// Add nodes to the nodes dictonary of the model
			for (int i = 0; i < nodes.Count; ++i)
			{
				model.NodesDictionary.Add(i + 1, nodes[i]);
			}

			// Hexa8NonLinear element definition
			var hexa8NLelement = new Element()
			{
				ID = 1,
				ElementType = new Hexa8NonLinear(solidMaterial, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))
			};

			// Add nodes to the created element
			hexa8NLelement.AddNode(model.NodesDictionary[node8.ID]);
			hexa8NLelement.AddNode(model.NodesDictionary[node7.ID]);
			hexa8NLelement.AddNode(model.NodesDictionary[node5.ID]);
			hexa8NLelement.AddNode(model.NodesDictionary[node6.ID]);
			hexa8NLelement.AddNode(model.NodesDictionary[node4.ID]);
			hexa8NLelement.AddNode(model.NodesDictionary[node3.ID]);
			hexa8NLelement.AddNode(model.NodesDictionary[node1.ID]);
			hexa8NLelement.AddNode(model.NodesDictionary[node2.ID]);

			// Add Hexa element to the element and subdomains dictionary of the model
			model.ElementsDictionary.Add(hexa8NLelement.ID, hexa8NLelement);
			model.SubdomainsDictionary[subdomainID].Elements.Add(hexa8NLelement);

			// Boundary Condtitions
			for (int iNode = 1; iNode <= 4; iNode++)
			{
				model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
				model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
				model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}

			// Boundary Condtitions - Prescribed DOFs
			for (int iNode = 5; iNode <= 8; iNode++)
			{
				model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ, Amount = nodalDisplacement });
			}

			// Choose linear equation system solver
			var solverBuilder = new SkylineSolver.Builder();
			SkylineSolver solver = solverBuilder.BuildSolver(model);

			// Choose the provider of the problem -> here a structural problem
			var provider = new ProblemStructural(model, solver);

			// Choose child analyzer -> Child: DisplacementControlAnalyzer
			var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
			var childAnalyzerBuilder = new DisplacementControlAnalyzer.Builder(model, solver, provider, increments)
			{
				MaxIterationsPerIncrement = 50,
				NumIterationsForMatrixRebuild = 1,
				ResidualTolerance = 1E-03,
			};
			var childAnalyzer = childAnalyzerBuilder.Build();

			// Choose parent analyzer -> Parent: Static
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 0 });

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
			Assert.Equal(-1.019828463478385, log.DOFValues[0], 8);
		}
	}
}
