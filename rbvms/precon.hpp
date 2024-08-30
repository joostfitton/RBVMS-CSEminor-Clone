// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_PRECON_HPP
#define RBVMS_PRECON_HPP

#include "mfem.hpp"
#include "tau.hpp"

using namespace mfem;

namespace RBVMS
{

// Custom block preconditioner for the Jacobian
class JacobianPreconditioner : public
   BlockLowerTriangularPreconditioner //BlockDiagonalPreconditioner
{
protected:
   Array<Solver *> prec;
public:
   JacobianPreconditioner(Array<int> &offsets, Array<Solver *> p)
      : BlockLowerTriangularPreconditioner (offsets), prec(p)
   { MFEM_VERIFY(offsets.Size()-1 == p.Size(), ""); };

   virtual void SetOperator(const Operator &op);

   virtual ~JacobianPreconditioner();
};

} // namespace RBVMS

#endif
