## Example 3: Beam3D Quaternion Nonlinear Test

The present test demonstrates the geometrically nonlinear beahvior of 3-D corotational
beam portal frame structure and can be found in **Beam3DQuaternionNonlinearTest.PlaneFrameTest()**.
The plane frame structure consists of three (3) elements and four (4) nodes.

Define the mechanical properties, nodal loads, number of nodes, number of elements and 
monitor node:
```csharp
double youngModulus = 21000.0;
double poissonRatio = 0.3;
double nodalLoad = 500000.0;
double area = 91.04;
double inertiaY = 2843.0;
double inertiaZ = 8091.0;
double torsionalInertia = 76.57;
double effectiveAreaY = 91.04;
double effectiveAreaZ = 91.04;
int nNodes = 4;
int nElems = 3;
int monitorNode = 2;
```

Generate nodal geometry:
```csharp
// Node creation
IList<Node> nodes = new List<Node>();
Node node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
Node node2 = new Node(id: 2, x: 0.0, y: 100.0, z: 0.0);
Node node3 = new Node(id: 3, x: 100.0, y: 100.0, z: 0.0);
Node node4 = new Node(id: 4, x: 100.0, y: 0.0, z: 0.0);

nodes.Add(node1);
nodes.Add(node2);
nodes.Add(node3);
nodes.Add(node4);
```

Define the model and the subdomain:
```csharp
// Model creation
Model model = new Model();

// Add a single subdomain to the model
model.SubdomainsDictionary.Add(1, new Subdomain(1));

// Add nodes to the nodes dictionary of the model
for (int i = 0; i < nodes.Count; ++i)
{
	model.NodesDictionary.Add(i + 1, nodes[i]);
}
```

Constrain the first and the last nodes:
```csharp
// Constrain first and last nodes of the model
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
```

Generate the elements of the structure:
```csharp
// Generate elements of the structure
int iNode = 1;
for (int iElem = 0; iElem < nElems; iElem++)
{
	// element nodes
	IList<Node> elementNodes = new List<Node>();
	elementNodes.Add(model.NodesDictionary[iNode]);
	elementNodes.Add(model.NodesDictionary[iNode + 1]);

	// Create new Beam3D section and element
	var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);
	var beam = new Beam3DCorotationalQuaternion(elementNodes, material, 7.85, beamSection);

	// Create elements
	var element = new Element()
	{
		ID = iElem + 1,
		ElementType = beam
	};

	var a = beam.StiffnessMatrix(element);

	// Add nodes to the created element
	element.AddNode(model.NodesDictionary[iNode]);
	element.AddNode(model.NodesDictionary[iNode + 1]);

	// Add beam element to the element and subdomains dictionary of the model
	model.ElementsDictionary.Add(element.ID, element);
	model.SubdomainsDictionary[1].Elements.Add(element);
	iNode++;
}
```

Loading conditions:
```csharp
// Add nodal load values at the top nodes of the model
model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode], DOF = StructuralDof.TranslationX });
```

Define Skyline solver, Structural provider, Load Control (i.e Newton-Raphson) chlid analyzer and
Static parent analyzer:
```csharp
// Solver
var solverBuilder = new SkylineSolver.Builder();
ISolver solver = solverBuilder.BuildSolver(model);

// Problem type
var provider = new ProblemStructural(model, solver);

// Analyzers
int increments = 10;
var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
childAnalyzerBuilder.ResidualTolerance = 1E-3;
//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
```

Output request and run the analysis:
```csharp
// Request output
childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 0 });

parentAnalyzer.Initialize();
parentAnalyzer.Solve();
```

Check output results:
```csharp
// Check output
DOFSLog log = (DOFSLog)childAnalyzer.Logs[1][0];
Assert.Equal(120.1108698752, log.DOFValues[0], 2);
```



