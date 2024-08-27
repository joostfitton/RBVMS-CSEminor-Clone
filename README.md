# RBVMS

 of the residual-based variational multiscale turbulence modeling (RBVMS) methodology at high Reynolds number.
 We show that the RBVMS formulation globally conserves angular momentum,
 a feature that is felt to be important for flows dominated by rotation, and that is not shared by standard stabilized 
formulations of fluid flow. Weak imposition of Dirichlet boundary conditions is employed to enhance the accuracy of the 
RBVMS framework in the presence of thin turbulent boundary layers near solid walls. 
Calculation of conservative boundary forces and torques is also presented for the case of weakly enforced boundary conditions.
 NURBS-based isogeometric analysis is employed for the spatial discretization, and 
mesh refinement is performed to assess the convergence characteristics of the proposed methodology. 
Numerical tests show that very accurate results are obtained on relatively coarse grids. To the best of the authors’ knowledge, 
this paper is the first to report large eddy simulation computations of this challenging test case.


Y. Bazilevs, I. Akkerman,
Large eddy simulation of turbulent Taylor–Couette flow using isogeometric analysis and the residual-based variational multiscale method,
Journal of Computational Physics,
Volume 229, Issue 9, 2010, Pages 3402-3414,ISSN 0021-9991,
https://doi.org/10.1016/j.jcp.2010.01.008.
[https://www.sciencedirect.com/science/article/pii/S0021999110000239][1]



## Install

`sudo apt-get install cmake build-essential openmpi metis`


<cd build
cmake ..>

## Run

### Isogeometric

### GMSH
