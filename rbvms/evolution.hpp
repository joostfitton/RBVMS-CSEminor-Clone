// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_EVOLUTION_HPP
#define RBVMS_EVOLUTION_HPP

#include "mfem.hpp"
#include "weakform.hpp"

using namespace std;
using namespace mfem;

namespace RBVMS
{

class ParTimeDepBlockNonlinForm;

class Evolution : public TimeDependentOperator
{
private:
   ParTimeDepBlockNonlinForm &form;
   Solver &solver;
   Vector dudt;

public:
   Evolution(ParTimeDepBlockNonlinForm &form,
             Solver &solver);

   virtual void ImplicitSolve(const real_t dt,
                              const Vector &x,
                              Vector &k) override;

   ~Evolution() {}
};


/** Residual-based Variational multiscale integrator
    for incompressible Navier-Stokes flow

    This is a specialized ParBlockNonlinearForm that includes timestepping
    interpolation, vis:
       add(x0,dt,dx,x);   // x = x0 + dt*dx

    Both x and dx need to be passed to the FormIntegrator.
    To avoid defining these interfaces this class is subsumed.
*/

class ParTimeDepBlockNonlinForm : public ParBlockNonlinearForm
{
private:
   RBVMS::IncNavStoIntegrator &integrator;


   /// Numerical parameters
   real_t dt;

   /// Solution & Residual vector
   mutable Vector x0, x;
   mutable BlockVector dxs;
   mutable BlockVector dxs_true;

public:
   /// Constructor
   ParTimeDepBlockNonlinForm(Array<ParFiniteElementSpace *> &pfes,
                             RBVMS::IncNavStoIntegrator &integrator);

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


   virtual BlockOperator &GetGradient(const Vector &x) const;

   /// Return the local block gradient matrix for the given true-dof vector x
   const BlockOperator& GetLocalGradient(const Vector &x) const;

   /// Specialized version of GetGradient() for BlockVector
   void ComputeGradientBlocked(const BlockVector &bx,
                               const BlockVector &dbx) const;
};

} // namespace RBVMS

#endif

