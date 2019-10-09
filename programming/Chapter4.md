## Example 4: Newmark Dynamic Analysis Example

The present example demonstrates the behavior of the Newmark dynamic analyzer. To this purpose
Example 9.4 from [1] is solved. 

The Mocking technique will be used, since no elements are used in this example. This method 
will help us define the mass, stiffness and damping matrices along with the initial conditions
without generating any particular element or structure.

The initial model subdomain is generated according to the following:
```csharp
var model = new Model();
int subdomainID = 0;
model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
```

Create an empty node and an element without nodes:
```csharp
var n = new Node(id: 0, x: double.NaN);
var e = new Element() { ID = 0 };
e.NodesDictionary.Add(0, n);
```

Create a "mock" finite element with the following matrices as preperties:
```csharp
var m = new Mock<IFiniteElement>();
m.Setup(x => x.StiffnessMatrix(e)).Returns(Matrix.CreateFromArray(new double[,] { { 6, -2, }, { -2, 4 } }));
m.Setup(x => x.MassMatrix(e)).Returns(Matrix.CreateFromArray(new double[,] { { 2, 0, }, { 0, 1 } }));
m.Setup(x => x.DampingMatrix(e)).Returns(Matrix.CreateFromArray(new double[,] { { 0, 0, }, { 0, 0 } }));
m.Setup(x => x.GetElementDofTypes(e)).Returns(new[] { new[] { StructuralDof.TranslationX, StructuralDof.TranslationY } });
m.SetupGet(x => x.DofEnumerator).Returns(new GenericDofEnumerator());
```

And pass these properties to the element:
```csharp
e.ElementType = m.Object;
```

Define Nodes and Elements Disctionaries to be inserted to the model and the subdomain:
```csharp
model.NodesDictionary.Add(0, n);
model.ElementsDictionary.Add(0, e);
model.SubdomainsDictionary[subdomainID].Elements.Add(e);
model.Loads.Add(new Load() { Amount = 10, Node = n, DOF = StructuralDof.TranslationY });
```

Calculate mass accelaration history loads along x and y-axis:
```csharp
var lX = new Mock<IMassAccelerationHistoryLoad>();
lX.SetupGet(x => x.DOF).Returns(StructuralDof.TranslationX);
lX.SetupGet(x => x[It.IsAny<int>()]).Returns(0);
var lY = new Mock<IMassAccelerationHistoryLoad>();
lY.SetupGet(x => x.DOF).Returns(StructuralDof.TranslationY);
lY.SetupGet(x => x[0]).Returns(10);
lY.SetupGet(x => x[It.IsInRange(1, 100, Range.Inclusive)]).Returns(0);
model.MassAccelerationHistoryLoads.Add(lX.Object);
model.MassAccelerationHistoryLoads.Add(lY.Object);
m.Setup(x => x.CalculateAccelerationForces(It.IsAny<Element>(), It.IsAny<IList<MassAccelerationLoad>>()))
	.Returns<Element, IList<MassAccelerationLoad>>((element, loads) =>
	{
		double[] accelerations = { loads[0].Amount, loads[1].Amount };
		var massMatrix = Matrix.CreateFromArray(new double[,] { { 2, 0, }, { 0, 1 } });
		return massMatrix.Multiply(accelerations);
	}
);
```

Define Skyline solver, Structural provider, Linear child analyzer and Newmark Dynamic
parent anayzer:
```csharp
// Solver
var solverBuilder = new SkylineSolver.Builder();
ISolver solver = solverBuilder.BuildSolver(model);

// Problem type
var provider = new ProblemStructural(model, solver);

// Analyzers
var childAnalyzer = new LinearAnalyzer(model, solver, provider);
var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.28, 3.36);
parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();
```

Request output and run the analysis:
```csharp
// Request output
childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 0, 1 });

// Run the analysis
parentAnalyzer.Initialize();
parentAnalyzer.Solve();
```

## References
[1] Finite Element Procedures, K-J. Bathe, Prentice Hall, 2nd edition, 2014.

