// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#ifndef RBVMS_EVOLUTION_HPP
#define RBVMS_EVOLUTION_HPP

#include "mfem.hpp"
#include "weakform.hpp"

using namespace std;
using namespace mfem;

namespace RBVMS
{

class Evolution : public TimeDependentOperator
{
private:
   ParBlockNonlinearForm form;
   Solver &solver;
   IncNavStoIntegrator *integrator;

public:
   Evolution(Array<ParFiniteElementSpace *> spaces,
             Array<Array<int> *> ess_bdr,
             Solver &solver,
             IncNavStoIntegrator *integrator);

   virtual void ImplicitSolve(const real_t dt,
                              const Vector &x,
                              Vector &k) override;

   virtual ~Evolution();
};

} // namespace RBVMS

#endif

