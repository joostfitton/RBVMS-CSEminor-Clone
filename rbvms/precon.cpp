// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "precon.hpp"

using namespace mfem;
using namespace RBVMS;

// Set the diagonal and off-diagonal operators
void JacobianPreconditioner::SetOperator(const Operator &op)
{
    BlockOperator *jacobian = (BlockOperator *)&op;

    // Set up the preconditioner for the first block (AMG)
    if (prec[0] == nullptr)
    {
        prec[0] = new HypreBoomerAMG(); // Use HYPRE AMG for the first block
    }

    // Set up the preconditioner for the remaining blocks (ILU)
    for (int i = 1; i < prec.Size(); ++i)
    {
        if (prec[i] == nullptr)
        {
            prec[i] = new HypreILU(); // Use HYPRE ILU for other blocks
        }
    }

    // Configure the blocks
    for (int i = 0; i < prec.Size(); ++i)
    {
        prec[i]->SetOperator(jacobian->GetBlock(i, i));
        SetDiagonalBlock(i, prec[i]);

        for (int j = i + 1; j < prec.Size(); ++j)
        {
            SetBlock(j, i, const_cast<Operator *>(&jacobian->GetBlock(j, i)));
        }
    }
}

{
   for (int i = 0; i < prec.Size(); ++i)
   {
      delete prec[i];
   }
}
