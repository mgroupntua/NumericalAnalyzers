## Example 1: Euler Beam2D Linear Test

The following example demostrates the linear, static beahavior of a cantelever, 2-dimensional beam structure under bending.
The examples splits in three parts. The first part refers to the model creation, the second to the loading and boundary conditions
and the third and final defines the analysis type. 

At first we define the mechanical properties of the beam cross-section and the material of the structure,
according to the following piece of code:

```csharp
// Define mechanical properties
double youngModulus = 21000.0;
double poissonRatio = 0.3;
double nodalLoad = 2000.0;
int nElems = 2;
int monitorNode = 3;

// Create new material
var material = new ElasticMaterial()
{
	YoungModulus = youngModulus,
	PoissonRatio = poissonRatio,
};
```

The geometry of the structure consists of two (2) elements along the x-axis, with a total of three (3) nodes and
six (6) degrees of freedom (DOFs). The geometry of the nodes is give by:

```csharp
// Define geometry: Node creation. 
IList<Node> nodes = new List<Node>();
Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
Node node2 = new Node(id: 2, x: 100.0, y: 0.0, z: 0.0);
Node node3 = new Node(id: 3, x: 200.0, y: 0.0, z: 0.0);
nodes.Add(node1);
nodes.Add(node2);
nodes.Add(node3);
```

A model and a sub-domain are then generated as:

```csharp
// Model creation.
Model model = new Model();

// Add a single sub-domain to the model.
model.SubdomainsDictionary.Add(1, new Subdomain(1));
```

and the nodes are inserted to the model by a dictionary according to the following:

```csharp
// Add nodes to the nodes dictionary of the model.
for (int i = 0; i < nodes.Count; ++i)
{
  model.NodesDictionary.Add(i + 1, nodes[i]);
}
```

The boundary condtions are applied on each node using its key and on the degree of freedom 
that is constrained (translational along x,y,z-axes or rotational about x,y,z-axes), as:

```csharp
// Define Boundary Conditions: Constrain bottom nodes of the model.
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
```

In order to define a beam element, a beam-section has to be defined first and contains the mechanical properties.
Then each element has to be defined by two (2) nodes, since the EulerBeam2D is a two-noded element. Finally, the elements
are added in an elements dictionary and passed into the model and the subdomain that we wish.

```csharp
// Generate elements of the structure.
int iNode = 1;
for (int iElem = 0; iElem < nElems; iElem++)
{
	// Create new Beam2D section and element.
	var beam = new EulerBeam2D(youngModulus)
	{
		Density = 7.85,
		SectionArea = 91.04,
		MomentOfInertia = 8091.00,
	};

	// Create elements.
	var element = new Element()
	{
		ID = iElem + 1,
		ElementType = beam
	};

	// Add nodes to the created element.
	element.AddNode(model.NodesDictionary[iNode]);
	element.AddNode(model.NodesDictionary[iNode + 1]);
				
	// Show element stiffness matrix.
	var a = beam.StiffnessMatrix(element);

	// Adds Beam2D element to the element and sub-domains dictionary of the model.
	model.ElementsDictionary.Add(element.ID, element);
	model.SubdomainsDictionary[1].Elements.Add(element);
	iNode++;
}
```

The loading conditions are defined by the keay of the node and the direction of the load, as:

```csharp
// Add nodal load values at the top nodes of the model.
model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = StructuralDof.TranslationY });
```

In final part of the examples we have to define the problem type to be solved and the type of the solver that we 
want to use. In this case we want to use a Skyline solver, while our problem type is structural and the analyzers 
must combine a linear and a static one, which can be done by the definition of child and a parent analyzer, as can be seen
in the following code lines:

```csharp
// Defines Skyline Solver.
var solverBuilder = new SkylineSolver.Builder();
ISolver solver = solverBuilder.BuildSolver(model);

// Defines Problem type as Structural.
var provider = new ProblemStructural(model, solver);

// Defines Analyzers: Child Analyzer => Linear, Parent Analyzer => Static
var childAnalyzer = new LinearAnalyzer(model, solver, provider);
var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
```

Before the analysis starts, it is mendatory to initialize the parent analyzer:

```csharp
// Run the analysis.
parentAnalyzer.Initialize();
parentAnalyzer.Solve();
```

In the following lines, we just demonstrate the solution norm and the residual (right-hand-side vector) norm, as:

```csharp
// Solution and right-hand-side norms.
double solutionNorm = solver.LinearSystems[1].Solution.Norm2();
double rhsNorm = solver.LinearSystems[1].RhsVector.Norm2();
```

While the final check of the test is to compare the y-axis dispacement of the free node of the structure with 
the real solution:

```csharp
// Check solution.
Assert.Equal(31.388982074929341, solver.LinearSystems[1].Solution[4], 12);
```

Since this critetion is satisfied, the test is valid and this means that the solution we get is equal or close enough to the 
rright solution.



