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

/** This class defines the time-dependent integrator for the
    Residual-based Variational multiscale formulation
    for incompressible Navier-Stokes flow.
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
   mutable Vector u, dudt, f, grad_p, res_m, up, nor;
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

   /// Set the timestep size @a dt_
   void SetTimeAndStep(const real_t &t, const real_t &dt_)
   {
      dt = dt_;
      c_mu.SetTime(t);
      c_force.SetTime(t);
   };

   /// Find the local CFL-number
   real_t GetElementCFL(const Array<const FiniteElement *>&el,
                        ElementTransformation &Tr) const;

   /// Assemble the local energy
   real_t GetElementEnergy(const Array<const FiniteElement *>&el,
                           ElementTransformation &Tr,
                           const Array<const Vector *> &elfun,
                           const Array<const Vector *> &elrate) const;

   /// Assemble the element interior residual vectors
   void AssembleElementVector(const Array<const FiniteElement *> &el,
                              ElementTransformation &Tr,
                              const Array<const Vector *> &elsol,
                              const Array<const Vector *> &elrate,
                              const Array<Vector *> &elvec) const;

   /// Assemble the element interior gradient matrices
   void AssembleElementGrad(const Array<const FiniteElement*> &el,
                            ElementTransformation &Tr,
                            const Array<const Vector *> &elsol,
                            const Array<const Vector *> &elrate,
                            const Array2D<DenseMatrix *> &elmats) const;


   /// Assemble the outflow boundary residual vectors
   void AssembleOutflowVector(const Array<const FiniteElement *> &el1,
                              const Array<const FiniteElement *> &el2,
                              FaceElementTransformations &Tr,
                              const Array<const Vector *> &elfun,
                              const Array<Vector *> &elvect);

   /// Assemble the outflow boundary gradient matrices
   void AssembleOutflowGrad(const Array<const FiniteElement *>&el1,
                            const Array<const FiniteElement *>&el2,
                            FaceElementTransformations &Tr,
                            const Array<const Vector *> &elfun,
                            const Array2D<DenseMatrix *> &elmats);


   /// Assemble the weak Dirichlet BC boundary residual vectors
   void AssembleWeakDirBCVector(const Array<const FiniteElement *> &el1,
                                const Array<const FiniteElement *> &el2,
                                FaceElementTransformations &Tr,
                                const Array<const Vector *> &elfun,
                                const Array<Vector *> &elvect);

   /// Assemble the weak Dirichlet BC boundary gradient matrices
   void AssembleWeakDirBCGrad(const Array<const FiniteElement *>&el1,
                              const Array<const FiniteElement *>&el2,
                              FaceElementTransformations &Tr,
                              const Array<const Vector *> &elfun,
                              const Array2D<DenseMatrix *> &elmats);
};

} // namespace RBVMS

#endif
