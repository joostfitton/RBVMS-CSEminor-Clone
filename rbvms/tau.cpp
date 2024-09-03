// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "tau.hpp"

using namespace mfem;
using namespace RBVMS;

/* Evaluate the stabilisation parameter at @a ip.
    - tau[0] approximates the momentum operator and
      is used for velocity reconstruction.
    - tau[1] approximates the pressure poisson operator and
      is used for pressure reconstruction.*/
void Tau::Eval(Vector &tau,
               ElementTransformation &T,
               const IntegrationPoint &ip)
{
   real_t mu = c_mu.Eval(T, ip);
   c_adv.Eval(u, T, ip);

   MultAAt(T.InverseJacobian(),Gij);

   tau[0] = Ct/(dt*dt);
   for (int j = 0; j < dim; j++)
      for (int i = 0; i < dim; i++)
      {
         tau[0]  += Gij(i,j)*u[i]*u[j];
      }

   for (int j = 0; j < dim; j++)
      for (int i = 0; i < dim; i++)
      {
         tau[0]  += Cd*Gij(i,j)*Gij(i,j)*mu*mu;
      }

   tau[0]  = 1.0/sqrt(tau[0]);
   tau[1]  = 1.0/(tau[0]*Gij.Trace());
}
