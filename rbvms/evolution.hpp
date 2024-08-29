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
   Vector zero;

public:
   ParBlockNonlinearForm &form;
   Solver &solver;

   Evolution(ParBlockNonlinearForm &form,
             Solver &solver);

   virtual void ImplicitSolve(const double dt,
                              const Vector &x,
                              Vector &k) override;

   virtual ~Evolution();
};

} // namespace RBVMS

#endif

