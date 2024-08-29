
// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#include "evolution.hpp"

using namespace mfem;
using namespace RBVMS;

//--------------------------------------------------------------------------
Evolution::Evolution(ParBlockNonlinearForm &form,
                     Solver &solver)
   : TimeDependentOperator(form->Width()),0,0, IMPLICIT),
     form(form),solver(solver)
{
   zero.SetSize(0);
}

Evolution::~Evolution()
{

}

void Evolution::ImplicitSolve(const double dt,
                              const Vector &u0, Vector &dudt)
{
  // form->SetParameters(dt, &u0);
   solver.Mult(zero, dudt);
}
