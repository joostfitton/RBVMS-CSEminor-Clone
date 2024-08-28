# RBVMS

Application that can simulate incompressible flow
using the residual-based variational multiscale (RBVMS) formulation.
This novel formulation is also suitable as a Large-Eddy Simulation model at turbulent Reynolds numbers.

Novel features include:
* Formulation combining numerical and physical arguments
* * Local Momentum conservation
* * Local Angular momentum conservation
* Wide range of discretisation possible, including
* * Tetrahedral
* * Isogeometric
* Weak imposition of Dirichlet boundary conditions
* Conservative boundary forces extraction
* Symmetry BC???

For details on the formulation see:

Y. Bazilevs, I. Akkerman,
Large eddy simulation of turbulent Taylor–Couette flow using isogeometric analysis and the residual-based variational multiscale method,
Journal of Computational Physics,
Volume 229, Issue 9, 2010, Pages 3402-3414,ISSN 0021-9991,
https://doi.org/10.1016/j.jcp.2010.01.008.

The parallel application is build on [mfem](https://github.com/mfem/mfem), a
lightweight, general, scalable C++ library for finite element methods.

## Install

To install RBVMS on a Ubuntu-based machine make sure the appropriate packages are installed.
This can be done using the following command:

```
sudo apt-get install cmake git build-essential openmpi metis
```

The RBVMS code can be downloaded using the following command

```
git clone git@github.com:IdoAkkerman/rbvms.git
```

This will create a directory `rbvms` with the code in it. This code can be compiled using the following commands:

```
cd rbvms
mkdir build
cd build
cmake ..
make -j 16
```

Here the number 16 is the number of files that will be compiled concurrently.
NB: This number you choose should depend on the number of cores you have available.

Depending on the machine compilation might take a few minutes.
After this hypre and MFEM will be installed in the `build/external` directory, and
the rbvms executable can be found in the `bin` directory.


## Run

```
XXX/bin/rbvms
```

### Isogeometric




### GMSH

GMSH is an open-source mesher. On Ubuntu-based machine the package can be installed
using the following command:

```
sudo apt-get install gmsh
```

# License

MFEM is distributed under the terms of the BSD-3 license. All new contributions must be made under this license. See LICENSE and NOTICE for details.

SPDX-License-Identifier: BSD-3-Clause
LLNL Release Number: LLNL-CODE-806117
DOI: 10.11578/dc.20171025.1248
