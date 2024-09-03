// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_NAVSTO_HPP
#define RBVMS_NAVSTO_HPP

#include "mfem.hpp"
#include "tau.hpp"

using namespace mfem;

namespace RBVMS
{

/** Residual-based Variational multiscale integrator
    for incompressible Navier-Stokes flow
*/
class IncNavStoIntegrator
{
private:

   // Physical parameters
   Coefficient &c_mu;
   VectorCoefficient &c_force;

   /// Numerical parameters
   real_t dt = -1.0;
   Tau &tau_m;

   /// Dimension data
   int dim = -1;
   Array2D<int> hmap;

   /// Physical values
   mutable Vector u, dudt, f, grad_p, res_m, up;
   mutable DenseMatrix flux;

   /// Solution & Residual vector
   mutable DenseMatrix elf_u, elf_du, elv_u;

   /// Shape function data
   mutable Vector sh_u, ushg_u, sh_p, dupdu;
   mutable DenseMatrix shg_u, shh_u, shg_p, grad_u, hess_u;

public:
   /// Constructor
   IncNavStoIntegrator(Coefficient &mu_,
                       VectorCoefficient &force_,
                       Tau &tau);


   void SetTimeStep(real_t &dt_) {dt = dt_; };

   /// Assemble the local energy
   real_t GetElementEnergy(const Array<const FiniteElement *>&el,
                           ElementTransformation &Tr,
                           const Array<const Vector *> &elfun,
                           const Array<const Vector *> &elrate) const;

   /// Assemble the local residual vectors
   void AssembleElementVector(const Array<const FiniteElement *> &el,
                              ElementTransformation &Tr,
                              const Array<const Vector *> &elsol,
                              const Array<const Vector *> &elrate,
                              const Array<Vector *> &elvec) const;

   /// Assemble the local gradient matrices
   void AssembleElementGrad(const Array<const FiniteElement*> &el,
                            ElementTransformation &Tr,
                            const Array<const Vector *> &elsol,
                            const Array<const Vector *> &elrate,
                            const Array2D<DenseMatrix *> &elmats) const;
};

} // namespace RBVMS

#endif
