// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_TAU_HPP
#define RBVMS_TAU_HPP

#include "mfem.hpp"

using namespace mfem;

namespace RBVMS
{

/** This Class defines the stabilisation parameter for reconstructing the
    small-scale velocity.
*/
class Tau: public VectorCoefficient
{
protected:
   /// The advection field
   VectorCoefficient &c_adv;

   /// The diffusion parameter field
   Coefficient &c_mu;
   real_t dt = -1.0;

   /// Coefficients
   real_t Cd, Ct;

   /// Dimension of the problem
   int dim;

   /// Velocity vector
   Vector u;

   // Metric tensor
   DenseMatrix Gij;

public:
   /** Construct a stabilized confection-diffusion integrator with:
       - @a adv the convection velocity
       - @a mu the diffusion coefficient*/
   Tau(VectorCoefficient &adv, Coefficient &mu,
       real_t Cd = 36.0, real_t Ct = 4.0)
      : VectorCoefficient(2), c_adv(adv), c_mu(mu),
        Cd(Cd), Ct(Ct)
   {
      dim = c_adv.GetVDim();
      u.SetSize(dim);
      Gij.SetSize(dim);
   };

   void SetTimeStep(const double &dt_) { dt = dt_; };

   /// Evaluate the coefficient at @a ip.
   virtual void Eval(Vector &V,
                     ElementTransformation &T,
                     const IntegrationPoint &ip) override;

};

} // namespace RBVMS

#endif
