
// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#include "evolution.hpp"

using namespace mfem;
using namespace RBVMS;

Evolution::Evolution(IncNavStoForm &form,
                     Solver &solver)
   : TimeDependentOperator(form.Width(), 0.0, IMPLICIT),
     form(form), solver(solver), dudt(form.Width())
{
   solver.SetOperator(form);
   dudt = 0.0; 
}

void Evolution::ImplicitSolve(const real_t dt,
                              const Vector &u0, Vector &dudt_)
{
   form.SetSolution(dt, u0);
   Vector zero;
   solver.Mult(zero, dudt);
   dudt_ = dudt;
}
