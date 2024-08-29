// This file is part of the RBVMS application. For more information and source code
// availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.

#include "weakform.hpp"

using namespace mfem;
using namespace RBVMS;

IncNavStoIntegrator::IncNavStoIntegrator(Coefficient &mu,
                                         VectorCoefficient &force,
                                         Tau &tm, Tau &tc, Tau &tb)
   : c_mu(mu), c_force(force), tau_m(tm), tau_c(tc), tau_b(tb)
{
   dim = force.GetVDim();
   u.SetSize(dim);
   f.SetSize(dim);
   res.SetSize(dim);
   up.SetSize(dim);
   grad_u.SetSize(dim);
   hess_u.SetSize(dim, (dim*(dim+1))/2);
   grad_p.SetSize(dim);
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

real_t IncNavStoIntegrator::GetElementEnergy(
   const Array<const FiniteElement *>&el,
   ElementTransformation &Tr,
   const Array<const Vector *>&elfun)
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::GetElementEnergy"
                 " has incorrect block finite element space size!");
   }
   int dof_u = el[0]->GetDof();

   sh_u.SetSize(dof_u);
   elf_u.UseExternalData(elfun[0]->GetData(), dof_u, dim);

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
   }

   return energy;
}

void IncNavStoIntegrator::AssembleElementVector(
   const Array<const FiniteElement *> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elfun,
   const Array<Vector *> &elvec)
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

   elf_u.UseExternalData(elfun[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

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

      el[0]->CalcPhysDShape(Tr, shg_u);
      shg_u.Mult(u, ushg_u);
      MultAtB(elf_u, shg_u, grad_u);

      if (hess)
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
      }
      else
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elfun[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elfun[1], grad_p);

      // Compute strong residual
      grad_u.Mult(u,res);   // Add convection
      res += grad_p;        // Add pressure
      res -= f;             // Subtract force
      for (int i = 0; i < dim; ++i)
      {
         for (int j = 0; j < dim; ++j)
         {
            res[j] -= mu*(hess_u(j,hmap(i,i)) +
                          hess_u(i,hmap(j,i))); // Add diffusion
         }
      }

      // Compute stability params
      real_t tm = tau_m.Eval(Tr, ip);
      real_t tc = tau_c.Eval(Tr, ip);

      // Compute momentum weak residual
      flux.Diag(-p + tc*grad_u.Trace(),dim);  // Add pressure & LSIC to flux
      grad_u.Symmetrize();                   // Grad to strain
      flux.Add(2*mu,grad_u);                 // Add stress to flux
      AddMult_a_VVt(-1.0, u, flux);          // Add convection to flux
      AddMult_a_VWt(tm, res, u,
                    flux);        // Add SUPG to flux --> check order u and res
      AddMult_a_ABt(w, shg_u, flux, elv_u);  // Add flux term to rhs
      AddMult_a_VWt(-w, sh_u, f,    elv_u);  // Add force term to rhs

      // Compute momentum weak residual
      elvec[1]->Add(w*grad_u.Trace(), sh_p); // Add Galerkin term
      shg_p.Mult(res, sh_p);                 // PSPG help term
      elvec[1]->Add(w*tm, sh_p);              // Add PSPG term  - sign looks worng?
   }
}

void IncNavStoIntegrator::AssembleElementGrad(
   const Array<const FiniteElement*> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elfun,
   const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   bool hess = false;// = (el[0]->GetDerivType() == (int) FiniteElement::HESS);

   elf_u.UseExternalData(elfun[0]->GetData(), dof_u, dim);

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
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      real_t mu = c_mu.Eval(Tr, ip);
      real_t tm = tau_m.Eval(Tr, ip);
      real_t tc = tau_c.Eval(Tr, ip);

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);

      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elfun[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elfun[1], grad_p);

      // u,u block
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
            mat *= mu;

            // Convection
            mat -= ushg_u(i_u)*sh_u(j_u);      // Galerkin
            mat += tm*ushg_u(i_u)*ushg_u(j_u);  // SUPG

            mat *= w;
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,0))(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += mat;
            }

            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               for (int j_dim = 0; j_dim < dim; ++j_dim)
               {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) +=
                     (mu + tc)*shg_u(i_u,j_dim)*shg_u(j_u,i_dim)*w;
               }
            }
         }
      }

      // u,p and p,u blocks
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,1))(j_u + dof_u * dim_u, i_p) += (shg_p(i_p, dim_u)*tm*ushg_u(j_u)
                                                            -shg_u(j_u,dim_u)*sh_p(i_p))*w;
               (*elmats(1,0))(i_p, j_u + dof_u * dim_u) +=  shg_u(j_u,dim_u)*sh_p(i_p)*w;
            }
         }
      }

      // p,p block
      AddMult_a_AAt(w*tm, shg_p, *elmats(1,1));
   }
}
