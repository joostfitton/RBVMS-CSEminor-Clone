// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#include "monitor.hpp"

using namespace mfem;


void GeneralResidualMonitor::MonitorResidual(int it, real_t norm,
                                             const Vector &r, bool final)
{
   if (it == 0)
   {
      norm0 = norm;
   }

   if ((print_level > 0 &&  it%print_level == 0) || final)
   {
      mfem::out << prefix << " iteration " << std::setw(2) << it
                << " : ||r|| = " << norm
                << ",  ||r||/||r_0|| = " << 100*norm/norm0<<" % \n";
   }
}

void SystemResidualMonitor::MonitorResidual(int it, real_t norm,
                                            const Vector &r, bool final)
{
   if (dc && (it > 0))
   {
      if (rank > 1)
      {
         for (int i = 0; i < nvar; ++i)
         {
            pgf[i]->Distribute(xp->GetBlock(i));
         }
      }
      dc->SetCycle(it);
      dc->Save();
   }

   Vector vnorm(nvar);

   for (int i = 0; i < nvar; ++i)
   {
      Vector r_i(r.GetData() + bOffsets[i], bOffsets[i+1] - bOffsets[i]);
      if ( rank == 1 )
      {
         vnorm[i] = r_i.Norml2();
      }
      else
      {
         vnorm[i] = sqrt(InnerProduct(MPI_COMM_WORLD, r_i, r_i));
      }
      if (it == 0) { norm0[i] = vnorm[i]; }
   }

   bool print = (print_level > 0 &&  it%print_level == 0) || final;
   if (print)
   {
      mfem::out << prefix << " iteration " << std::setw(3) << it <<"\n"
                << " ||r||  \t"<< "||r||/||r_0||  \n";
      for (int i = 0; i < nvar; ++i)
      {
         mfem::out <<vnorm[i]<<"\t"<< 100*vnorm[i]/norm0[i]<<" % \n";
      }
   }
}


void JacobianPreconditioner::SetOperator(const Operator &op)
{
   BlockOperator *jacobian = (BlockOperator *) &op;

   for (int i = 0; i < prec.Size(); ++i)
   {
      prec[i]->SetOperator(jacobian->GetBlock(i,i));
      SetDiagonalBlock(i, prec[i]);
   }

   SetBlock(1,0, const_cast<Operator*>(&jacobian->GetBlock(1,0)));
}

JacobianPreconditioner::~JacobianPreconditioner()
{
   for (int i = 0; i < prec.Size(); ++i)
   {
      delete prec[i];
   }
}

