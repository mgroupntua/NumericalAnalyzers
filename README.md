![alt text](http://mgroup.ntua.gr/wp-content/uploads/2018/05/MGroup52.png "MGroup")


# NumericalAnalyzers

## Features

The present documentation provides information about the analyzers that MSolve
contains in order to solve static or dynamic and geometrically linear or nonlinear 
problems.

The response of linear static systems can be calculated using the following
procedures:

- Direct Solution:
    - Cholesky factorization  
- Iterative Solution:
    - Preconditioned Conjugate Gradient (PCG) method

More information about the above mentioned methods can be found in the section of
**Solvers**.

The solution of nonlinear static equations can be obtained by:

- The Newton-Raphson (or Load Control) Method and
- The Displacement Control Mathod

The solution of equilibrium equations in dynamic analysis can be found using:

- Direct Methods:
    - Central differences
    - Newmark
- Nonlinear Dynamic Equations:
    - Implicit (i.e. Newmark)
- Nonstructural Problems:
    - Heat Tranfer (i.e. Central Differences)

## Installation instructions
You can choose either to clone the solution or downloads it as a zip file.

### Clone solution
1. Under the repository name, click **Clone or Download** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. In the popup appearing choose the **Use HTTPS** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/2.png "2")

3. Use the ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/3.png "3") to copy the link provided.

4. Open Visual Studio. In Team Explorer window appearing in your screen under Local Git Repositories click the **Clone** option. If Team Explorer window is not visible you can enable in View -> Team Explorer

  ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/4.png "4")
  
5. In the text box appearing paste the link.

 ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/5.png "5")

6. Click clone and Visual Studio will automatically download and import **MSolve.Analyzers**


### Download as ZIP
1. Under the repository name, click **Clone or Download** option

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. Click **Download ZIP** option. **MGroup.NumericalAnalyzers** will be downloaded as a ZIP file.

3. Extract the ZIP file to the folder of choice.

4. Double click on **MGroup.NumericalAnalyzers.sln** file to open the code with Visual Studio





