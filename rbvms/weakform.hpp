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

    This is a specialized ParBlockNonlinearForm that includes timestepping
    interpolation, vis:
       add(x0,dt,dx,x);   // x = x0 + dt*dx

    Both x and dx need to be passed to the FormIntegrator.
    To avoid defining these interfaces this class is subsumed.
*/

class IncNavStoForm : public ParBlockNonlinearForm
{
private:

   // Physical parameters
   Coefficient &c_mu;
   VectorCoefficient &c_force;

   /// Numerical parameters
   real_t dt;
   Tau &tau_m, &tau_c, &tau_b;

   /// Dimension data
   int dim = -1;
   Array2D<int> hmap;

   /// Physical values
   mutable Vector u, dudt, f, grad_p, res, up;
   mutable DenseMatrix flux;

   /// Solution & Residual vector
   mutable Vector x0, x;
   mutable DenseMatrix elf_u, elf_du, elv_u;

   /// Shape function data
   mutable Vector sh_u, ushg_u, sh_p;
   mutable DenseMatrix shg_u, shh_u, shg_p, grad_u, hess_u;

   ///
   mutable BlockVector dxs;
   mutable BlockVector dxs_true;

public:
   /// Constructor
   IncNavStoForm(Array<ParFiniteElementSpace *> &pfes,
                 Coefficient &mu_,
                 VectorCoefficient &force_,
                 Tau &tau_m, Tau &tau_c, Tau &tau_b);

   // Set the solution
   void SetSolution(const real_t dt,
                    const Vector &u0);

   //------------------------------------------------
   // NonlinearForm memberfunctions
   //------------------------------------------------
   /// Specialized version of GetEnergy() for BlockVectors
   //real_t GetEnergyBlocked(const BlockVector &bx) const;

   /// Block T-Vector to Block T-Vector
   void Mult(const Vector &x, Vector &y) const;

   /// Specialized version of Mult() for BlockVector%s
   /// Block L-Vector to Block L-Vector
   void MultBlocked(const BlockVector &bx,
                    const BlockVector &dbx,
                          BlockVector &by) const;

   /// Return the local block gradient matrix for the given true-dof vector x
   const BlockOperator& GetLocalGradient(const Vector &x) const;

   /// Specialized version of GetGradient() for BlockVector
   void ComputeGradientBlocked(const BlockVector &bx,
                               const BlockVector &dbx) const;

   //------------------------------------------------
   // Integrator memberfunctions
   //------------------------------------------------
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
