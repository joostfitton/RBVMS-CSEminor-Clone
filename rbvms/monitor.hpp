// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------


#ifndef RBVMS_MONITOR_HPP
#define RBVMS_MONITOR_HPP

#include "mfem.hpp"
#include "tau.hpp"

using namespace mfem;

namespace RBVMS
{

/**


*/

class GeneralResidualMonitor : public IterativeSolverMonitor
{
public:
   GeneralResidualMonitor(const std::string& prefix_, int print_lvl);

   GeneralResidualMonitor(MPI_Comm comm,
                          const std::string& prefix_, int print_lvl);

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
                         DataCollection *dc_ = nullptr);

   SystemResidualMonitor(MPI_Comm comm,
                         const std::string& prefix_,
                         int print_lvl,
                         Array<int> &offsets);

   SystemResidualMonitor(MPI_Comm comm,
                         const std::string& prefix_,
                         int print_lvl,
                         Array<int> &offsets,
                         DataCollection *dc_,
                         BlockVector *x,
                         Array<ParGridFunction *> pgf_);

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

} // namespace RBVMS

#endif
