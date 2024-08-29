// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#ifndef RBVMS_MONITOR_HPP
#define RBVMS_MONITOR_HPP

#include "mfem.hpp"
#include "tau.hpp"

namespace mfem
{

/**


*/


class GeneralResidualMonitor : public IterativeSolverMonitor
{
public:
   GeneralResidualMonitor(const std::string& prefix_, int print_lvl)
      : prefix(prefix_)
   {
      print_level = print_lvl;
      rank = 1;
   }

   GeneralResidualMonitor(MPI_Comm comm,
                          const std::string& prefix_, int print_lvl)
      : prefix(prefix_)
   {
#ifndef MFEM_USE_MPI
      print_level = print_lvl;
#else
      MPI_Comm_rank(comm, &rank);
      if (rank == 0)
      {
         print_level = print_lvl;
      }
      else
      {
         print_level = -1;
      }
#endif
   }

   virtual void MonitorResidual(int it, real_t norm, const Vector &r, bool final);

private:
   const std::string prefix;
   int rank, print_level;
   mutable real_t norm0;
};

class SystemResidualMonitor : public IterativeSolverMonitor
{
public:
   SystemResidualMonitor(const std::string& prefix_,
                         int print_lvl,
                         Array<int> &offsets,
                         DataCollection *dc_ = nullptr)
      : prefix(prefix_), bOffsets(offsets), dc(dc_)
   {
      print_level = print_lvl;
      nvar = bOffsets.Size()-1;
      norm0.SetSize(nvar);
      rank = 1;
   }

   SystemResidualMonitor(MPI_Comm comm,
                         const std::string& prefix_,
                         int print_lvl,
                         Array<int> &offsets)
      : prefix(prefix_), bOffsets(offsets), dc(nullptr), xp(nullptr)
   {
#ifndef MFEM_USE_MPI
      print_level = print_lvl;
      rank = 1;
#else
      MPI_Comm_rank(comm, &rank);
      if (rank == 0)
      {
         print_level = print_lvl;
      }
      else
      {
         print_level = -1;
      }
#endif
      nvar = bOffsets.Size()-1;
      norm0.SetSize(nvar);
   }
   SystemResidualMonitor(MPI_Comm comm,
                         const std::string& prefix_,
                         int print_lvl,
                         Array<int> &offsets,
                         DataCollection *dc_,
                         BlockVector *x,
                         Array<ParGridFunction *> pgf_)
      : prefix(prefix_), bOffsets(offsets), dc(dc_), xp(x), pgf(pgf_)
   {
#ifndef MFEM_USE_MPI
      print_level = print_lvl;
      rank = 1;
#else
      MPI_Comm_rank(comm, &rank);
      if (rank == 0)
      {
         print_level = print_lvl;
      }
      else
      {
         print_level = -1;
      }
#endif
      nvar = bOffsets.Size()-1;
      norm0.SetSize(nvar);
   }


   virtual void MonitorResidual(int it, real_t norm, const Vector &r, bool final);

private:
   const std::string prefix;
   int print_level, nvar, rank;
   mutable Vector norm0;
   // Offsets for extracting block vector segments
   Array<int> &bOffsets;
   DataCollection *dc;
   BlockVector *xp;
   Array<ParGridFunction *> pgf;
};

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


} // namespace mfem

#endif
