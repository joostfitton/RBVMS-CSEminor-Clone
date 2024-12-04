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
void ExtendedPreconditioner::SetOperator(const Operator &op)
{
   BlockOperator *jacobian = (BlockOperator *) &op;

   for (int i = 0; i < prec.Size(); ++i)
   {
      if (!prec[i])
      {
         // Choose preconditioner type for each block
         if (i == 0)
         {
            prec[i] = new HypreAMG(); // Use AMG for the first block
         }
         else
         {
            prec[i] = new HypreILU(); // Use ILU for other blocks
         }
      }

      prec[i]->SetOperator(jacobian->GetBlock(i, i));
      SetDiagonalBlock(i, prec[i]);

      for (int j = i + 1; j < prec.Size(); ++j)
      {
         SetBlock(j, i, const_cast<Operator *>(&jacobian->GetBlock(j, i)));
      }
   }
}

// Destructor
ExtendedPreconditioner::~ExtendedPreconditioner()
{
   for (int i = 0; i < prec.Size(); ++i)
   {
      delete prec[i];
   }
}
