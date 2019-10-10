## 1. Nonlinear Analyzers

### 1.1 Load Control Method

The following notes describe the load control method as described by the Newton-Raphson 
scheme. The nonlinear equation of equilibrium can be expressed in terms of the out-of-balance
force vector, as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{r}&space;=&space;\mathbf{0}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{r}&space;=&space;\mathbf{0}" title="\mathbf{r} = \mathbf{0}" /></a>

or in terms of external and internal force vectors:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}^{j}&space;=&space;\mathbf{f}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}^{j}&space;=&space;\mathbf{f}^{j}" title="\mathbf{R}^{j} = \mathbf{f}^{j}" /></a>

where the superscript j indicates that the external nodal forces are applied 
incrementally. The external loading vector **R** is independent from the displacements, 
while the internal force vector **f** is nonlinearly dependent from the displacements 
vector **u**. Using the Taylor series expansion for eq.1, we get:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{r}(\mathbf{u}_{i}^{j})&space;=&space;\mathbf{r}(\mathbf{u}_{i-1}^{j})&space;&plus;&space;\left[\frac{\partial&space;\mathbf{r}}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;\delta&space;\mathbf{u}_{i}^{j}&space;&plus;&space;..." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{r}(\mathbf{u}_{i}^{j})&space;=&space;\mathbf{r}(\mathbf{u}_{i-1}^{j})&space;&plus;&space;\left[\frac{\partial&space;\mathbf{r}}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;\delta&space;\mathbf{u}_{i}^{j}&space;&plus;&space;..." title="\mathbf{r}(\mathbf{u}_{i}^{j}) = \mathbf{r}(\mathbf{u}_{i-1}^{j}) + \left[\frac{\partial \mathbf{r}}{\partial \mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}} \delta \mathbf{u}_{i}^{j} + ..." /></a>

where:

<a href="https://www.codecogs.com/eqnedit.php?latex=\left[\frac{\partial&space;\mathbf{r}}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;=&space;\left[\frac{\partial&space;\mathbf{r}(\mathbf{u}_{i-1}^{j})}{\partial&space;\mathbf{u}_{i-1}^{j}}\right]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left[\frac{\partial&space;\mathbf{r}}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;=&space;\left[\frac{\partial&space;\mathbf{r}(\mathbf{u}_{i-1}^{j})}{\partial&space;\mathbf{u}_{i-1}^{j}}\right]" title="\left[\frac{\partial \mathbf{r}}{\partial \mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}} = \left[\frac{\partial \mathbf{r}(\mathbf{u}_{i-1}^{j})}{\partial \mathbf{u}_{i-1}^{j}}\right]" /></a>

and the i-th iterative solution at the j-th load increment is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta&space;\mathbf{u}_{i}^{j}&space;=&space;\mathbf{u}_{i}^{j}&space;-&space;\mathbf{u}_{i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta&space;\mathbf{u}_{i}^{j}&space;=&space;\mathbf{u}_{i}^{j}&space;-&space;\mathbf{u}_{i-1}^{j}" title="\delta \mathbf{u}_{i}^{j} = \mathbf{u}_{i}^{j} - \mathbf{u}_{i-1}^{j}" /></a>

According to the above, we obtain that for the (i-1)th iteration, 
the out-of-balance force vector becomes:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{r}_{i-1}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{r}_{i-1}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}" title="\mathbf{r}_{i-1}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}" /></a>

Combining the above, we get:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}&space;&plus;&space;\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;\delta&space;\mathbf{u}_{i}^{j}&space;=&space;\mathbf{0}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}&space;&plus;&space;\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;\delta&space;\mathbf{u}_{i}^{j}&space;=&space;\mathbf{0}" title="\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j} + \left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial \mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}} \delta \mathbf{u}_{i}^{j} = \mathbf{0}" /></a>

that can also be written as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;\delta&space;\mathbf{u}_{i}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}&space;\delta&space;\mathbf{u}_{i}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}" title="\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial \mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}} \delta \mathbf{u}_{i}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}" /></a>

and takes the final incremental-iterative form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{K}_{i-1}^{j}\delta&space;\mathbf{u}_{i}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}&space;\Rightarrow&space;\mathbf{K}_{i-1}^{j}\delta\mathbf{u}_{i}^{j}&space;=&space;\mathbf{r}_{i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{K}_{i-1}^{j}\delta&space;\mathbf{u}_{i}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j}&space;\Rightarrow&space;\mathbf{K}_{i-1}^{j}\delta\mathbf{u}_{i}^{j}&space;=&space;\mathbf{r}_{i-1}^{j}" title="\mathbf{K}_{i-1}^{j}\delta \mathbf{u}_{i}^{j}=\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j} \Rightarrow \mathbf{K}_{i-1}^{j}\delta\mathbf{u}_{i}^{j} = \mathbf{r}_{i-1}^{j}" /></a>

where:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{K}_{i-1}^{j}=\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{K}_{i-1}^{j}=\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial&space;\mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}" title="\mathbf{K}_{i-1}^{j}=\left[\frac{\partial(\mathbf{R}^{j}-\mathbf{f}_{i-1}^{j})}{\partial \mathbf{u}}\right]_{\mathbf{u}_{i-1}^{j}}" /></a>

is the tangent stiffness matrix.

**Incremental-iterative implementation**

The iterative solution is calculated as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{u}_{i}^{j}&space;=&space;\left[&space;\mathbf{K}_{i-1}^{j}&space;\right]^{-1}\mathbf{r}_{i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{u}_{i}^{j}&space;=&space;\left[&space;\mathbf{K}_{i-1}^{j}&space;\right]^{-1}\mathbf{r}_{i-1}^{j}" title="\delta\mathbf{u}_{i}^{j} = \left[ \mathbf{K}_{i-1}^{j} \right]^{-1}\mathbf{r}_{i-1}^{j}" /></a>

and the total solution up to this point becomes:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{u}_{i}^{j}&space;=&space;\mathbf{u}_{i-1}^{j}&space;&plus;&space;\delta&space;\mathbf{u}_{i}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{u}_{i}^{j}&space;=&space;\mathbf{u}_{i-1}^{j}&space;&plus;&space;\delta&space;\mathbf{u}_{i}^{j}" title="\mathbf{u}_{i}^{j} = \mathbf{u}_{i-1}^{j} + \delta \mathbf{u}_{i}^{j}" /></a>

Since the total displacements vector **u** is known, 
we can estimate the total internal force **f**, 
at the i-th iteration. According to the above, 
the updated residual force vector at this iteration can be calculated from the relation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{r}_{i}^{j}&space;=&space;\mathbf{R}^{j}&space;-&space;\mathbf{f}_{i}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{r}_{i}^{j}&space;=&space;\mathbf{R}^{j}&space;-&space;\mathbf{f}_{i}^{j}" title="\mathbf{r}_{i}^{j} = \mathbf{R}^{j} - \mathbf{f}_{i}^{j}" /></a>

Convergence has been achieved when the norm of the residual vector is small enough 
compared to the norm of the actual internal forces. 
This can be described by the following equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\|\mathbf{r}_{i}^{j}\|}{\|\mathbf{R}_{i}^{j}\|}&space;<&space;tolerance" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\|\mathbf{r}_{i}^{j}\|}{\|\mathbf{R}_{i}^{j}\|}&space;<&space;tolerance" title="\frac{\|\mathbf{r}_{i}^{j}\|}{\|\mathbf{R}_{i}^{j}\|} < tolerance" /></a>

where the tolerance indicates the convergence of our method and in practice it takes 
values between 10^6^ and 10^3^.

### 1.2 Displacement Control Method

The following notes describe the displacement control method for multiple prescribed DOFs as 
referred to BeBorst book. 

Equilibrium equation is written as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{K}_{i-1}^{j}\delta\mathbf{u}_{i}^{j}&space;=&space;\mathbf{r}_{i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{K}_{i-1}^{j}\delta\mathbf{u}_{i}^{j}&space;=&space;\mathbf{r}_{i-1}^{j}" title="\mathbf{K}_{i-1}^{j}\delta\mathbf{u}_{i}^{j} = \mathbf{r}_{i-1}^{j}" /></a>

where:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{r}_{i-1}^{j}=\mathbf{R}_{i}^{j}-\mathbf{f}_{i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{r}_{i-1}^{j}=\mathbf{R}_{i}^{j}-\mathbf{f}_{i-1}^{j}" title="\mathbf{r}_{i-1}^{j}=\mathbf{R}_{i}^{j}-\mathbf{f}_{i-1}^{j}" /></a>

is the residual force vector.
The displacements vector is partitioned according to "free" and "prescribed" DOFs, as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{u}_{i}^{j}=[\delta\mathbf{u}_{f,i}^{j},&space;\delta\mathbf{u}_{p,i}^{j}]^{T}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{u}_{i}^{j}=[\delta\mathbf{u}_{f,i}^{j},&space;\delta\mathbf{u}_{p,i}^{j}]^{T}" title="\delta\mathbf{u}_{i}^{j}=[\delta\mathbf{u}_{f,i}^{j}, \delta\mathbf{u}_{p,i}^{j}]^{T}" /></a>

The partitioned stiffness matrix, external loading vector and internal forces vector can be written as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{K}_{i-1}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{K}_{ff,i-1}^{j}&space;&&space;\mathbf{K}_{fp,i-1}^{j}&space;\\&space;\mathbf{K}_{pf,i-1}^{j}&space;&&space;\mathbf{K}_{pp,i-1}^{j}&space;\\&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{K}_{i-1}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{K}_{ff,i-1}^{j}&space;&&space;\mathbf{K}_{fp,i-1}^{j}&space;\\&space;\mathbf{K}_{pf,i-1}^{j}&space;&&space;\mathbf{K}_{pp,i-1}^{j}&space;\\&space;\end{bmatrix}" title="\mathbf{K}_{i-1}^{j} = \begin{bmatrix} \mathbf{K}_{ff,i-1}^{j} & \mathbf{K}_{fp,i-1}^{j} \\ \mathbf{K}_{pf,i-1}^{j} & \mathbf{K}_{pp,i-1}^{j} \\ \end{bmatrix}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}_{i}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{R}_{f,i}^{j}&space;\\&space;\mathbf{R}_{p,i}^{j}&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}_{i}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{R}_{f,i}^{j}&space;\\&space;\mathbf{R}_{p,i}^{j}&space;\end{bmatrix}" title="\mathbf{R}_{i}^{j} = \begin{bmatrix} \mathbf{R}_{f,i}^{j} \\ \mathbf{R}_{p,i}^{j} \end{bmatrix}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{f}_{i-1}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{f}_{f,i-1}^{j}&space;\\&space;\mathbf{f}_{p,i-1}^{j}&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{f}_{i-1}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{f}_{f,i-1}^{j}&space;\\&space;\mathbf{f}_{p,i-1}^{j}&space;\end{bmatrix}" title="\mathbf{f}_{i-1}^{j} = \begin{bmatrix} \mathbf{f}_{f,i-1}^{j} \\ \mathbf{f}_{p,i-1}^{j} \end{bmatrix}" /></a>

According to the above, the system of equations to be solved becomes:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{K}_{ff,i-1}^{j}\delta\mathbf{u}_{f,i}^{j}&space;&plus;&space;\mathbf{K}_{fp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}&space;=&space;\mathbf{R}_{f,i}^{j}&space;-&space;\mathbf{f}_{f,i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{K}_{ff,i-1}^{j}\delta\mathbf{u}_{f,i}^{j}&space;&plus;&space;\mathbf{K}_{fp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}&space;=&space;\mathbf{R}_{f,i}^{j}&space;-&space;\mathbf{f}_{f,i-1}^{j}" title="\mathbf{K}_{ff,i-1}^{j}\delta\mathbf{u}_{f,i}^{j} + \mathbf{K}_{fp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j} = \mathbf{R}_{f,i}^{j} - \mathbf{f}_{f,i-1}^{j}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{K}_{pf,i-1}^{j}\delta\mathbf{u}_{f,i}^{j}&space;&plus;&space;\mathbf{K}_{pp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}&space;=&space;\mathbf{R}_{p,i}^{j}&space;-&space;\mathbf{f}_{p,i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{K}_{pf,i-1}^{j}\delta\mathbf{u}_{f,i}^{j}&space;&plus;&space;\mathbf{K}_{pp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}&space;=&space;\mathbf{R}_{p,i}^{j}&space;-&space;\mathbf{f}_{p,i-1}^{j}" title="\mathbf{K}_{pf,i-1}^{j}\delta\mathbf{u}_{f,i}^{j} + \mathbf{K}_{pp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j} = \mathbf{R}_{p,i}^{j} - \mathbf{f}_{p,i-1}^{j}" /></a>

The incremental displacements vector of the "free" DOFs can be calculated from the
following relatioship as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{u}_{f,i}^{j}=(\mathbf{K}_{ff,i-1}^{j})^{-1}[\mathbf{R}_{f,i}^{j}-\mathbf{f}_{f,i-1}^{j}-\mathbf{K}_{fp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{u}_{f,i}^{j}=(\mathbf{K}_{ff,i-1}^{j})^{-1}[\mathbf{R}_{f,i}^{j}-\mathbf{f}_{f,i-1}^{j}-\mathbf{K}_{fp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}]" title="\delta\mathbf{u}_{f,i}^{j}=(\mathbf{K}_{ff,i-1}^{j})^{-1}[\mathbf{R}_{f,i}^{j}-\mathbf{f}_{f,i-1}^{j}-\mathbf{K}_{fp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}]" /></a>

The equivalent external loading of the prescribed DOFs can be now calculated as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}_{p,i}^{j}&space;=&space;\mathbf{K}_{pf,i-1}^{j}\delta\mathbf{u}_{f,i}^{j}&plus;\mathbf{K}_{pp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}&plus;\mathbf{f}_{p,i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}_{p,i}^{j}&space;=&space;\mathbf{K}_{pf,i-1}^{j}\delta\mathbf{u}_{f,i}^{j}&plus;\mathbf{K}_{pp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}&plus;\mathbf{f}_{p,i-1}^{j}" title="\mathbf{R}_{p,i}^{j} = \mathbf{K}_{pf,i-1}^{j}\delta\mathbf{u}_{f,i}^{j}+\mathbf{K}_{pp,i-1}^{j}\delta\mathbf{u}_{p,i}^{j}+\mathbf{f}_{p,i-1}^{j}" /></a>

Since the iterative displacements vector of the free DOFs is known, we have:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}_{p,i}^{j}&space;=&space;[\mathbf{K}_{pf,i-1}(\mathbf{K}_{ff,i-1}^{j})^{-1}\mathbf{K}_{fp,i-1}^{j}&plus;\mathbf{K}_{pp,i-1}^{j}]\delta\mathbf{u}_{p,i}^{j}&space;&plus;&space;\mathbf{K}_{pf,i-1}(\mathbf{K}_{ff,i-1}^{j})^{-1}(\mathbf{R}_{f,i}^{j}&space;-&space;\mathbf{f}_{f,i-1}^{j})&space;&plus;&space;\mathbf{f}_{p,i-1}^{j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}_{p,i}^{j}&space;=&space;[\mathbf{K}_{pf,i-1}(\mathbf{K}_{ff,i-1}^{j})^{-1}\mathbf{K}_{fp,i-1}^{j}&plus;\mathbf{K}_{pp,i-1}^{j}]\delta\mathbf{u}_{p,i}^{j}&space;&plus;&space;\mathbf{K}_{pf,i-1}(\mathbf{K}_{ff,i-1}^{j})^{-1}(\mathbf{R}_{f,i}^{j}&space;-&space;\mathbf{f}_{f,i-1}^{j})&space;&plus;&space;\mathbf{f}_{p,i-1}^{j}" title="\mathbf{R}_{p,i}^{j} = [\mathbf{K}_{pf,i-1}(\mathbf{K}_{ff,i-1}^{j})^{-1}\mathbf{K}_{fp,i-1}^{j}+\mathbf{K}_{pp,i-1}^{j}]\delta\mathbf{u}_{p,i}^{j} + \mathbf{K}_{pf,i-1}(\mathbf{K}_{ff,i-1}^{j})^{-1}(\mathbf{R}_{f,i}^{j} - \mathbf{f}_{f,i-1}^{j}) + \mathbf{f}_{p,i-1}^{j}" /></a>

The internal forces vector of the current iteration is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{f}_{i}^{j}=[\delta\mathbf{f}_{f,i}^{j},&space;\delta\mathbf{f}_{p,i}^{j}]^{T}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{f}_{i}^{j}=[\delta\mathbf{f}_{f,i}^{j},&space;\delta\mathbf{f}_{p,i}^{j}]^{T}" title="\delta\mathbf{f}_{i}^{j}=[\delta\mathbf{f}_{f,i}^{j}, \delta\mathbf{f}_{p,i}^{j}]^{T}" /></a>

While, the residual becomes:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{r}_{i}^{j}=\mathbf{R}_{i}^{j}-\mathbf{f}_{i}^{j}=\begin{bmatrix}&space;\mathbf{R}_{f,i}^{j}&space;\\&space;\mathbf{R}_{p,i}^{j}&space;\\&space;\end{bmatrix}&space;-&space;\begin{bmatrix}&space;\mathbf{f}_{f,i}^{j}&space;\\&space;\mathbf{f}_{p,i}^{j}&space;\\&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{r}_{i}^{j}=\mathbf{R}_{i}^{j}-\mathbf{f}_{i}^{j}=\begin{bmatrix}&space;\mathbf{R}_{f,i}^{j}&space;\\&space;\mathbf{R}_{p,i}^{j}&space;\\&space;\end{bmatrix}&space;-&space;\begin{bmatrix}&space;\mathbf{f}_{f,i}^{j}&space;\\&space;\mathbf{f}_{p,i}^{j}&space;\\&space;\end{bmatrix}" title="\mathbf{r}_{i}^{j}=\mathbf{R}_{i}^{j}-\mathbf{f}_{i}^{j}=\begin{bmatrix} \mathbf{R}_{f,i}^{j} \\ \mathbf{R}_{p,i}^{j} \\ \end{bmatrix} - \begin{bmatrix} \mathbf{f}_{f,i}^{j} \\ \mathbf{f}_{p,i}^{j} \\ \end{bmatrix}" /></a>

and he convergence criterion will be the following:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\|\mathbf{r}_{i}^{j}\|}{\|\mathbf{R}_{i}^{j}\|}&space;<&space;tolerance" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\|\mathbf{r}_{i}^{j}\|}{\|\mathbf{R}_{i}^{j}\|}&space;<&space;tolerance" title="\frac{\|\mathbf{r}_{i}^{j}\|}{\|\mathbf{R}_{i}^{j}\|} < tolerance" /></a>

where the tolerance indicates the convergence of our method and in practice it takes 
values between 10^6^ and 10^3^.

In the case that there are no external loadings applied and on the structure, 
but only displacements applied on the prescribed DOFs, 
we can assume that **R**<sub>f,i</sub> <sup>j</sup>
and so the external loading vector of the method becomes:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}_{i}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{0}&space;\\&space;\mathbf{R}_{p,i}^{j}&space;\\&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}_{i}^{j}&space;=&space;\begin{bmatrix}&space;\mathbf{0}&space;\\&space;\mathbf{R}_{p,i}^{j}&space;\\&space;\end{bmatrix}" title="\mathbf{R}_{i}^{j} = \begin{bmatrix} \mathbf{0} \\ \mathbf{R}_{p,i}^{j} \\ \end{bmatrix}" /></a>

**Incremental-iterative implementation**

At the first iteration (i=1)  of the first increment (j=1) the internal forces 
of the "free" DOFs are zero and the initial displacements of the prescribed DOFs 
have a specific value, so the incremental displacements vector of the the 
"free" DOFs becomes:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{u}_{f,1}=(\mathbf{K}_{ff,0})^{-1}[\mathbf{R}_{f,1}-\mathbf{K}_{fp,0}\delta\mathbf{u}_{p,1}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{u}_{f,1}=(\mathbf{K}_{ff,0})^{-1}[\mathbf{R}_{f,1}-\mathbf{K}_{fp,0}\delta\mathbf{u}_{p,1}]" title="\delta\mathbf{u}_{f,1}=(\mathbf{K}_{ff,0})^{-1}[\mathbf{R}_{f,1}-\mathbf{K}_{fp,0}\delta\mathbf{u}_{p,1}]" /></a>

It should be noted that in δ**u**<sub>p,1</sub>
the boundary conditions have been taken into account,
which means that in some DOFs the value is zero, but the vector's norm 
is not zero due to the initial displacements that have been applied to 
the other prescribed DOFs.

The external forces of the "prescribed" DOFs are calculated according to the
relationship:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}_{p,1}&space;=&space;\mathbf{K}_{pf,0}\delta\mathbf{u}_{f,1}&plus;\mathbf{K}_{pp,0}\delta\mathbf{u}_{p,1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}_{p,1}&space;=&space;\mathbf{K}_{pf,0}\delta\mathbf{u}_{f,1}&plus;\mathbf{K}_{pp,0}\delta\mathbf{u}_{p,1}" title="\mathbf{R}_{p,1} = \mathbf{K}_{pf,0}\delta\mathbf{u}_{f,1}+\mathbf{K}_{pp,0}\delta\mathbf{u}_{p,1}" /></a>

since the initial internal forces of the "prescribed" DOFs are 
also zero, **f**<sub>p,0</sub>=0.

The external forces of the "prescribed" DOFs are the corresponding forces 
that occur in the "pre" DOFs as the result of the displacements of the 
"free" DOFs and the applied displacements of the "pre" DOFs themselves.

For i>=2, the incremental displacements of the "prescribed" DOFs must be equal
to zero:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{u}_{p,2}&space;=&space;\mathbf{0}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{u}_{p,2}&space;=&space;\mathbf{0}" title="\delta\mathbf{u}_{p,2} = \mathbf{0}" /></a>

so, the displacements vector becomes:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{u}_{2}&space;=\begin{bmatrix}&space;\delta\mathbf{u}_{f,2}&space;\\&space;\mathbf{0}&space;\\&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{u}_{2}&space;=\begin{bmatrix}&space;\delta\mathbf{u}_{f,2}&space;\\&space;\mathbf{0}&space;\\&space;\end{bmatrix}" title="\delta\mathbf{u}_{2} =\begin{bmatrix} \delta\mathbf{u}_{f,2} \\ \mathbf{0} \\ \end{bmatrix}" /></a>

where the incremental displacements of the "free" DOFs are:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta\mathbf{u}_{f,2}&space;=&space;\mathbf{K}_{ff,1}^{-1}[-\mathbf{f}_{f,1}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta\mathbf{u}_{f,2}&space;=&space;\mathbf{K}_{ff,1}^{-1}[-\mathbf{f}_{f,1}]" title="\delta\mathbf{u}_{f,2} = \mathbf{K}_{ff,1}^{-1}[-\mathbf{f}_{f,1}]" /></a>

and the external loading vector of the "prescribed" DOFs is calculated 
from the following equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{R}_{p,2}&space;=&space;\mathbf{K}_{pf,1}\delta\mathbf{u}_{f,2}&plus;\mathbf{f}_{p,1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{R}_{p,2}&space;=&space;\mathbf{K}_{pf,1}\delta\mathbf{u}_{f,2}&plus;\mathbf{f}_{p,1}" title="\mathbf{R}_{p,2} = \mathbf{K}_{pf,1}\delta\mathbf{u}_{f,2}+\mathbf{f}_{p,1}" /></a>

## References
[1] Finite Element Procedures, K-J. Bathe, Prentice Hall, 2nd edition, 2014.

[2] Non-linear Finite Element Analysis of Solids and Structures, R. de Borst, M.A. Crisfield, J.J.C. Remmers, C.V. Verhoosel, Wiley, 2nd Edition, 2012.
