//--------------------------------------------------------------
// Solution function for the lid driven cavity case
//
// To compile run, for example:
// mpicc -shared -o libldc.so -fPIC ldc.c
// gcc   -shared -o libldc.so -fPIC ldc.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>

double sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   double x = coord[0];
   double y = coord[1];
   double eps = 1.0e-8;

   if ( (y > 1.0-eps) &&
        ( (x > 0.0+eps) && (x < 1.0-eps) ) )
   {
      u[0] = 1.0;
      u[1] = 0.0;
   }
   else
   {
      u[0] = 0.0;
      u[1] = 0.0;
   }
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = 0.0;
}

