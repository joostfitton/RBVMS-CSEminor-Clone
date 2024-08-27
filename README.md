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



## Install

To install RBVMS on a Ubuntu-based machine make sure the appropriate packages are installed.
This can be done using the following command:

`sudo apt-get install cmake git build-essential openmpi metis`

The RBVMS code can be downloaded using the following command

` git clone git@github.com:IdoAkkerman/rbvms.git`

This will create a directory `rbvms` with the code in it. This code can be compiled using the following commands:

`
cd rbvms
mkdir build
cd build
cmake ..
make -j 16`

Here the number 16 is the number of files that will be compiled concurrently.
NB: This number you choose should depend on the number of cores you have available.

Depending on the machine compilation might take a few minutes.
After this hypre and MFEM will be installed in the `external` directory, and
the rbvms executable can be found in the `bin` directory.


## Run

`XXX/bin/rbvms `

### Isogeometric




### GMSH

GMSH is an open-source mesher. On Ubuntu-based machine the package can be installed
using the following command:

`sudo apt-get install gmsh`


