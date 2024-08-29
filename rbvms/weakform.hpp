// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

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
class IncNavStoIntegrator : public BlockNonlinearFormIntegrator
{
private:

   // Physical parameters
   Coefficient &c_mu;
   VectorCoefficient &c_force;

   /// Numerical parameters
   Tau &tau_m, &tau_c, &tau_b;

   /// Dimension data
   int dim = -1;
   Array2D<int> hmap;

   /// Physical values
   Vector u, f, grad_p, res, up;
   DenseMatrix flux;

   /// Solution & Residual vector
   DenseMatrix elf_u, elv_u;

   /// Shape function data
   Vector sh_u, ushg_u, sh_p;
   DenseMatrix shg_u, shh_u, shg_p, grad_u, hess_u;

public:
   /// Constructor
   IncNavStoIntegrator(Coefficient &mu_,
                       VectorCoefficient &force_,
                       Tau &tau_m, Tau &tau_c, Tau &tau_b);

   /// Assemble the local energy
   virtual real_t GetElementEnergy(const Array<const FiniteElement *>&el,
                                   ElementTransformation &Tr,
                                   const Array<const Vector *> &elfun);

   /// Assemble the local residual vectors
   virtual void AssembleElementVector(const Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      const Array<const Vector *> &elfun,
                                      const Array<Vector *> &elvec);

   /// Assemble the local gradient matrices
   virtual void AssembleElementGrad(const Array<const FiniteElement*> &el,
                                    ElementTransformation &Tr,
                                    const Array<const Vector *> &elfun,
                                    const Array2D<DenseMatrix *> &elmats);
};

} // namespace RBVMS

#endif
