// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "weakform.hpp"

using namespace mfem;
using namespace RBVMS;

// Constructor
IncNavStoIntegrator::IncNavStoIntegrator(Coefficient &mu,
                                         VectorCoefficient &force,
                                         VectorCoefficient &sol,
                                         Tau &tm)
   : c_mu(mu), c_force(force), c_sol(sol), tau_m(tm)
{
   dim = force.GetVDim();
   u.SetSize(dim);
   dudt.SetSize(dim);
   f.SetSize(dim);
   res_m.SetSize(dim);
   up.SetSize(dim);
   traction.SetSize(dim);
   grad_u.SetSize(dim);
   hess_u.SetSize(dim, (dim*(dim+1))/2);
   grad_p.SetSize(dim);
   nor.SetSize(dim);
   hmap.SetSize(dim,dim);

   if (dim == 2)
   {
      hmap(0,0) = 0;
      hmap(0,1) =  hmap(1,0) =  1;
      hmap(1,1) = 2;
   }
   else if (dim == 2)
   {
      hmap(0,0) = 0;
      hmap(0,1) = hmap(1,0) = 1;
      hmap(0,2) = hmap(2,0) = 2;
      hmap(1,1) = 3;
      hmap(1,2) = hmap(2,1) = 4;
      hmap(2,2) = 5;
   }
   else
   {
      mfem_error("Only implemented for 2D and 3D");
   }
}

// Get CFL
real_t IncNavStoIntegrator::GetElementCFL(const Array<const FiniteElement *>&el,
                                          ElementTransformation &Tr) const
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::GetElementCFL"
                 " has incorrect block finite element space size!");
   }
   Vector tau(3);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   real_t cfl = 0.0;
   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      tau_m.Eval(tau, Tr, ip);
      cfl = fmax(cfl, tau[2]);
   }

   return cfl;
}

// Get energy
real_t IncNavStoIntegrator::GetElementEnergy(
   const Array<const FiniteElement *>&el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate) const
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::GetElementEnergy"
                 " has incorrect block finite element space size!");
   }
   int dof_u = el[0]->GetDof();
   Vector tau(3);

   sh_u.SetSize(dof_u);
   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   real_t energy = 0.0;

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);

      real_t w = ip.weight * Tr.Weight();

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);

      energy += w*(u*u)/2;

      // Compute stability params
      tau_m.Eval(tau, Tr, ip);
   }

   return energy;
}

// Assemble the element interior residual vectors
void IncNavStoIntegrator::AssembleElementVector(
   const Array<const FiniteElement *> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array<Vector *> &elvec) const
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }
   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   Vector tau(3);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      real_t mu = c_mu.Eval(Tr, ip);
      c_force.Eval(f, Tr, ip);

      // Compute shape and interpolate
      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m += grad_p;        // Add pressure
      res_m -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int i = 0; i < dim; ++i)
         {
            for (int j = 0; j < dim; ++j)
            {
               res_m[j] -= mu*(hess_u(j,hmap(i,i)) +
                               hess_u(i,hmap(j,i)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }
      real_t res_c = grad_u.Trace();

      // Compute stability params
      tau_m.Eval(tau, Tr, ip);

      // Small scale reconstruction
      up.Set(-tau[0],res_m);
      u += up;
      p -= tau[1]*res_c;

      // Compute momentum weak residual
      flux.Diag(-p, dim);                         // Add pressure
      grad_u.Symmetrize();                        // Grad to strain
      flux.Add(2*mu,grad_u);                      // Add stress to flux
      AddMult_a_VVt(-1.0, u, flux);               // Add convection to flux
      AddMult_a_ABt(w, shg_u, flux, elv_u);       // Add flux term to rhs
      f -= dudt;                                  // Add Acceleration to force
      AddMult_a_VWt(-w, sh_u, f, elv_u);          // Add force + acc term to rhs

      // Compute continuity weak residual
      elvec[1]->Add(w*res_c, sh_p);               // Add Galerkin term
      shg_p.Mult(up, sh_p);                       // PSPG help term
      elvec[1]->Add(-w, sh_p);                    // Add PSPG term
   }
}

// Assemble the element interior gradient matrices
void IncNavStoIntegrator::AssembleElementGrad(
   const Array<const FiniteElement*> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array2D<DenseMatrix *> &elmats) const
{
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   bool hess = false;// = (el[0]->GetDerivType() == (int) FiniteElement::HESS);

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

   elmats(0,0)->SetSize(dof_u*dim, dof_u*dim);
   elmats(0,1)->SetSize(dof_u*dim, dof_p);
   elmats(1,0)->SetSize(dof_p, dof_u*dim);
   elmats(1,1)->SetSize(dof_p, dof_p);

   *elmats(0,0) = 0.0;
   *elmats(0,1) = 0.0;
   *elmats(1,0) = 0.0;
   *elmats(1,1) = 0.0;

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   Vector tau(3);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();

      real_t mu = c_mu.Eval(Tr, ip);

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);

      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m += grad_p;        // Add pressure
      res_m -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int i = 0; i < dim; ++i)
         {
            for (int j = 0; j < dim; ++j)
            {
               res_m[j] -= mu*(hess_u(j,hmap(i,i)) +
                               hess_u(i,hmap(j,i)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }

      // Compute stability params
      tau_m.Eval(tau, Tr, ip);

      // Small scale reconstruction
      up.Set(-tau[0],res_m);
      u += up;

      // Compute small scale jacobian
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         dupdu(j_u) = -tau[0]*(sh_u(j_u) + ushg_u(j_u)*dt);
      }

      // Recompute convective gradient
      MultAtB(elf_u, shg_u, grad_u);

      // Momentum - Velocity block (w,u)
      for (int i_u = 0; i_u < dof_u; ++i_u)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            // Diffusion
            real_t mat = 0.0;
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               mat += shg_u(i_u,dim_u)*shg_u(j_u,dim_u);
            }
            mat *= mu*dt;

            // Acceleration
            mat += sh_u(i_u)*sh_u(j_u);

            // Convection -- frozen convection
            mat -= ushg_u(i_u)*sh_u(j_u)*dt;           // Galerkin
            mat -= ushg_u(i_u)*dupdu(j_u);             // SUPG

            mat *= w;
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,0))(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += mat;
            }

            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               for (int j_dim = 0; j_dim < dim; ++j_dim)
               {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u)
                  += (mu + tau[1])*shg_u(i_u,i_dim)*shg_u(j_u,j_dim)*w*dt;
               }
            }
         }
      }

      // Momentum - Pressure block (w,p)
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,1))(j_u + dof_u * dim_u, i_p)
               += (shg_p(i_p,dim_u)*tau[0]*ushg_u(j_u)
                   - shg_u(j_u,dim_u)*sh_p(i_p))*w*dt;
            }
         }
      }

      // Continuity - Velocity block (q,u)
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(1,0))(i_p, j_u + dof_u * dim_u)
               += sh_p(i_p)*shg_u(j_u,dim_u)*w*dt;
               (*elmats(1,0))(i_p, j_u + dof_u * dim_u)
               -= shg_p(i_p, dim_u)*dupdu(j_u)*w;
            }
         }
      }

      // Continuity - Pressure block (w,p)
      AddMult_a_AAt(w*tau[0]*dt, shg_p, *elmats(1,1));
   }

}


// Assemble the outflow boundary residual vectors
void IncNavStoIntegrator
::AssembleOutflowVector(const Array<const FiniteElement *> &el1,
                        const Array<const FiniteElement *> &el2,
                        FaceElementTransformations &Tr,
                        const Array<const Vector *> &elsol,
                        const Array<Vector *> &elvec)
{
   int dof_u = el1[0]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   *elvec[0] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor);

      real_t w = ip.weight * 0.5;//* Tr.Weight();

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      real_t un = u*nor;
      AddMult_a_VWt(w*un, sh_u, u, elv_u);
   }
}

// Assemble the outflow boundary gradient matrices
void IncNavStoIntegrator
::AssembleOutflowGrad(const Array<const FiniteElement *>&el1,
                      const Array<const FiniteElement *>&el2,
                      FaceElementTransformations &Tr,
                      const Array<const Vector *> &elsol,
                      const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el1[0]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   elmats(0,0)->SetSize(dof_u*dim, dof_u*dim);
   *elmats(0,0) = 0.0;

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor);

      //real_t w = ip.weight * Tr.Weight(); //
      real_t w = ip.weight * 0.5;// instead???

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      real_t un = u*nor;

      // Momentum - Velocity block (w,u)
      for (int i_u = 0; i_u < dof_u; ++i_u)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            real_t mat = sh_u(i_u)*sh_u(j_u)*un*w*dt;

            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,0))(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += mat;
            }
         }
      }
   }
}


// Assemble the weak Dirichlet BC boundary residual vectors
void IncNavStoIntegrator
::AssembleWeakDirBCVector(const Array<const FiniteElement *> &el1,
                          const Array<const FiniteElement *> &el2,
                          FaceElementTransformations &Tr,
                          const Array<const Vector *> &elsol,
                          const Array<Vector *> &elvec)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      real_t mu = c_mu.Eval(*Tr.Elem1, eip);
      c_sol.Eval(up, *Tr.Elem1, eip);

      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      real_t w = ip.weight * Tr.Weight();

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      up -= u;
      up.Neg();
      real_t un = up*nor;

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Symmetrize();  // Grad to strain

      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);
      real_t p = sh_p*(*elsol[1]);

      Vector hn_vec(dim);
      //??T.InverseJacobian().Mult(nor,hn_vec);
      real_t Cb = 12.0;
      real_t lambda   = Cb*mu*hn_vec.Norml2();
      real_t lambda_n = 0.0;

      // Traction
      grad_u.Mult(nor, traction);
      traction *= -2*mu;              // Consistency
      traction.Add(p, nor);          // Pressure
      traction.Add(lambda,up);       // Penalty
      traction.Add(lambda_n*un,nor); // Penalty -- normal
      AddMult_a_VWt(w, sh_u, traction, elv_u);

      // Dual consistency
      MultVWt(nor,up, flux);
      flux.Symmetrize();
      AddMult_a_ABt(-w*2*mu, shg_u, flux, elv_u);

      // Continuity
      elvec[1]->Add(w*un, sh_p);
   }
}

// Assemble the weak Dirichlet BC boundary gradient matrices
void IncNavStoIntegrator
::AssembleWeakDirBCGrad(const
                        Array<const FiniteElement *>&el1,
                        const Array<const FiniteElement *>&el2,
                        FaceElementTransformations &Tr,
                        const Array<const Vector *> &elsol,
                        const Array2D<DenseMatrix *> &elmats)
{

}
