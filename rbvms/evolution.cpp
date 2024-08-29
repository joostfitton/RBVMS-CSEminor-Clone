
// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#include "evolution.hpp"

using namespace mfem;
using namespace RBVMS;

//--------------------------------------------------------------------------
Evolution::Evolution(Array<ParFiniteElementSpace *> spaces,
                     Array<Array<int> *> ess_bdr,
                     Solver &solver,
                     IncNavStoIntegrator *integrator)
   : TimeDependentOperator(0, 0.0, IMPLICIT),
     form(spaces),solver(solver), integrator(integrator)
{
   // Find size
   width = 0;
   for (int i=0; i<spaces.Size(); ++i)
   {
      width += spaces[i]->GetTrueVSize();
   }
   height = width;

   // Initialize formulation
   form.AddDomainIntegrator(integrator);
   Array<Vector *> rhs(2);
   rhs = nullptr;
   form.SetEssentialBC(ess_bdr, rhs);

   // Initialize solver
   solver.SetOperator(form);
}

Evolution::~Evolution()
{

}

void Evolution::ImplicitSolve(const real_t dt,
                              const Vector &u0, Vector &dudt)
{
   integrator->SetSolution(dt, u0);
   Vector zero;
   solver.Mult(zero, dudt);
}
