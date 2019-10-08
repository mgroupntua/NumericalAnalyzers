namespace MGroup.Analyzers.Tests
{
	using System.Collections.Generic;

	using MGroup.FEM.Elements;
	using MGroup.FEM.Entities;
	using MGroup.Materials;
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.FreedomDegrees;
	using MGroup.MSolve.Discretization.Loads;
	using MGroup.Problems;
	using MGroup.Solvers;
	using MGroup.Solvers.Direct;
	using Xunit;

	public class EulerBeam2DLinearTest
	{
		/// <summary>
		/// Example of Linear Beam2D Structure.
		/// Authors: George Soimiris
		/// </summary>
		[Fact]
		public void TestEulerBeam2DLinearBendingExample()
		{
			/// <summary>
			/// Define mechanical properties
			/// </summary>
			double youngModulus = 21000.0;
			double poissonRatio = 0.3;
			double nodalLoad = 2000.0;
			int nElems = 2;
			int monitorNode = 3;

			/// <summary>
			/// Create new material
			/// </summary>
			var material = new ElasticMaterial()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			/// <summary>
			/// Node creation. Define geometry.
			/// </summary>
			IList<Node> nodes = new List<Node>();
			Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			Node node2 = new Node(id: 2, x: 100.0, y: 0.0, z: 0.0);
			Node node3 = new Node(id: 3, x: 200.0, y: 0.0, z: 0.0);
			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);

			/// <summary>
			/// Model creation.
			/// </summary>
			Model model = new Model();

			/// <summary>
			/// Add a single subdomain to the model.
			/// </summary>
			model.SubdomainsDictionary.Add(1, new Subdomain(1));

			/// <summary>
			/// Add nodes to the nodes dictonary of the model.
			/// </summary>
			for (int i = 0; i < nodes.Count; ++i)
			{
				model.NodesDictionary.Add(i + 1, nodes[i]);
			}

			/// <summary>
			/// Constrain bottom nodes of the model. Define Boundary Conditions.
			/// </summary>
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
			model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });

			/// <summary>
			/// Generate elements of the structure.
			/// </summary>
			int iNode = 1;
			for (int iElem = 0; iElem < nElems; iElem++)
			{
				/// <summary>
				/// Create new Beam2D section and element.
				/// </summary>
				var beam = new EulerBeam2D(youngModulus)
				{
					Density = 7.85,
					SectionArea = 91.04,
					MomentOfInertia = 8091.00,
				};

				/// <summary>
				/// Create elements.
				/// </summary>
				var element = new Element()
				{
					ID = iElem + 1,
					ElementType = beam
				};

				/// <summary>
				/// Add nodes to the created element.
				/// </summary>
				element.AddNode(model.NodesDictionary[iNode]);
				element.AddNode(model.NodesDictionary[iNode + 1]);

				var a = beam.StiffnessMatrix(element);

				/// <summary>
				/// Adds Beam2D element to the element and subdomains dictionary of the model.
				/// </summary>
				model.ElementsDictionary.Add(element.ID, element);
				model.SubdomainsDictionary[1].Elements.Add(element);
				iNode++;
			}

			/// <summary>
			/// Add nodal load values at the top nodes of the model.
			/// </summary>
			model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = StructuralDof.TranslationY });

			/// <summary>
			/// Defines Skyline Solver.
			/// </summary>
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			/// <summary>
			/// Defines Problem type as Structural.
			/// </summary>
			var provider = new ProblemStructural(model, solver);

			/// <summary>
			/// Defines Analyzers.
			/// Chlid Analyzer: Linear
			/// Parent Analyzer: Static
			/// </summary>
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			/// <summary>
			/// Run the anlaysis.
			/// </summary>
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			/// <summary>
			/// Solution and right-hand-side norms.
			/// </summary>
			double solutionNorm = solver.LinearSystems[1].Solution.Norm2();
			double rhsNorm = solver.LinearSystems[1].RhsVector.Norm2();

			/// <summary>
			/// Check solution.
			/// </summary>
			Assert.Equal(31.388982074929341, solver.LinearSystems[1].Solution[4], 12);
		}
	}
}
