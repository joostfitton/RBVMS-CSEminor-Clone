// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
//   _____  ______      ____  __  _____
//   |  __ \|  _ \ \    / /  \/  |/ ____|
//   | |__) | |_) \ \  / /| \  / | (___
//   |  _  /|  _ < \ \/ / | |\/| |\___ \ 
//   | | \ \| |_) | \  /  | |  | |____) |
//   |_|  \_\____/   \/   |_|  |_|_____/
//
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "mfem.hpp"
#include "coefficients.hpp"
#include "weakform.hpp"
#include "evolution.hpp"
#include "precon.hpp"
#include "monitor.hpp"

using namespace std;
using namespace mfem;

extern void printInfo();

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Print relevant info
   printInfo();

   // 3. Parse command-line options.
   const char *mesh_file = "../../mfem/data/inline-quad.mesh";
   const char *ref_file  = "";
   int order = 1;
   int ref_levels = 0;

   real_t mu_param = 1.0;
   real_t penalty = -1.0;

   const char *lib_file = "libfun.so";

   int ode_solver_type = 35;
   real_t dt = 0.01;
   real_t t_final = 10.0;
   int vis_steps = 10;

   Array<int> strong_bdr;
   Array<int> weak_bdr;
   Array<int> outflow_bdr;

   OptionsParser args(argc, argv);

   // Mesh and discretization parameters
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order isoparametric space.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh.");

   // Problem parameters
   args.AddOption(&weak_bdr, "-wbc", "--weak-bdr",
                  "List of boundaries where Dirichelet BCs are enforced weakly."
                  "\n - Default: strong Dirichelet BCs");
   args.AddOption(&outflow_bdr, "-out", "--outflow-bdr",
                  "List of outflow boundaries.");
   args.AddOption(&lib_file, "-l", "--lib",
                  "Library file to use for IC, BC, force and solution function defintions.");
   args.AddOption(&mu_param, "-m", "--mu",
                  "Sets the diffusion parameters, should be positive.");

   // Solver parameters
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "Time integrator");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");

   // Parse parameters
   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (Mpi::Root()) { args.PrintOptions(cout); }

   // 4. Read the mesh from the given mesh file.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // Refine mesh
   {
      if (mesh.NURBSext && (strlen(ref_file) != 0))
      {
         mesh.RefineNURBSFromFile(ref_file);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
      if (Mpi::Root()) { mesh.PrintInfo(); }
   }

   // Partition mesh
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();

   // 5. Define a finite element space on the mesh.
   Array<FiniteElementCollection *> fecs(2);
   fecs[0] = new H1_FECollection(order, dim);
   fecs[1] = new H1_FECollection(order, dim);

   Array<ParFiniteElementSpace *> spaces(2);
   spaces[0] = new ParFiniteElementSpace(&pmesh, fecs[0],
                                         dim); //, Ordering::byVDIM);
   spaces[1] = new ParFiniteElementSpace(&pmesh, fecs[1]);

   // Report the degree of freedoms used
   {
      Array<int> tdof(num_procs),udof(num_procs),pdof(num_procs);
      tdof = 0;
      tdof[myid] = spaces[0]->TrueVSize();
      MPI_Reduce(tdof.GetData(), udof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

      tdof = 0;
      tdof[myid] = spaces[1]->TrueVSize();
      MPI_Reduce(tdof.GetData(), pdof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

      if (Mpi::Root())
      {
         mfem::out << "Number of finite element unknowns:\n";
         mfem::out << "\tVelocity = "<<spaces[0]->GlobalTrueVSize() << endl;
         mfem::out << "\tPressure = "<<spaces[1]->GlobalTrueVSize() << endl;
         mfem::out << "Number of finite element unknowns per partition:\n";
         mfem::out <<  "\tVelocity = "; udof.Print(mfem::out, num_procs);
         mfem::out <<  "\tPressure = "; pdof.Print(mfem::out, num_procs);
      }
   }

   // 6. Boundary conditions
   int amax = spaces[0]->GetMesh()->bdr_attributes.Max();
   Array<bool> bnd_flag(amax+1);
   bnd_flag = false;

   // Check weak boundaries
   if (Mpi::Root() && weak_bdr.Size() > 0 )
   {
      cout<<"Weak    = ";weak_bdr.Print();
   }
   for (int b = 0; b < weak_bdr.Size(); b++)
   {
      if ( weak_bdr[b] < 0 || weak_bdr[b] > amax )
      {
         mfem_error("Boundary out of range.");
      }
      if (bnd_flag[weak_bdr[b]])
      {
         mfem_error("Boundary specified more then once.");
      }
      bnd_flag[weak_bdr[b]] = true;
   }

   // Check outflow boundaries
   if (Mpi::Root() && outflow_bdr.Size() > 0)
   {
      cout<<"Outflow = "; outflow_bdr.Print();
   }
   for (int b = 0; b < outflow_bdr.Size(); b++)
   {
      if ( outflow_bdr[b] < 0 || outflow_bdr[b] > amax )
      {
         mfem_error("Boundary out of range.");
      }
      if (bnd_flag[outflow_bdr[b]])
      {
         mfem_error("Boundary specified more then once.");
      }
      bnd_flag[outflow_bdr[b]] = true;
   }

   // Assign strong boundaries
   strong_bdr.SetSize(amax - weak_bdr.Size() - outflow_bdr.Size());
   strong_bdr = -9;
   for (int b = 1, s = 0; b < amax+1; b++)
   {
      if (!bnd_flag[b])
      {
         strong_bdr[s] = b;
         s++;
      }
   }
   if (Mpi::Root() && strong_bdr.Size() > 0)
   {
      cout<<"Strong  = ";strong_bdr.Print();
   }

   // 7. Define the solution vector xp as a finite element grid function
   Array<int> bOffsets(3);
   bOffsets[0] = 0;
   bOffsets[1] = spaces[0]->TrueVSize();
   bOffsets[2] = spaces[1]->TrueVSize();
   bOffsets.PartialSum();

   BlockVector xp(bOffsets);

   // Define the gridfunctions
   ParGridFunction x_u(spaces[0]);
   ParGridFunction x_p(spaces[1]);

   LibVectorCoefficient sol(dim, lib_file, "sol_u");
   x_u.ProjectCoefficient(sol);
   x_p = 0.0;

   x_u.GetTrueDofs(xp.GetBlock(0));
   x_p.GetTrueDofs(xp.GetBlock(1));

   // Define the visualisation output
   VisItDataCollection visit_dc("navsto", &pmesh);
   visit_dc.RegisterField("u", &x_u);
   visit_dc.RegisterField("p", &x_p);
   visit_dc.SetCycle(0);
   visit_dc.Save();

   // 8. Define the time stepping algorithm
   // Set up the preconditioner
   RBVMS::JacobianPreconditioner jac_prec(bOffsets);

   // Set up the Jacobian solver
   RBVMS::GeneralResidualMonitor j_monitor(MPI_COMM_WORLD,"\t\tFGMRES", 25);
   FGMRESSolver j_gmres(MPI_COMM_WORLD);
   j_gmres.iterative_mode = false;
   j_gmres.SetRelTol(1e-3);
   j_gmres.SetAbsTol(1e-12);
   j_gmres.SetMaxIter(300);
   j_gmres.SetPrintLevel(-1);
   j_gmres.SetMonitor(j_monitor);
   j_gmres.SetPreconditioner(jac_prec);

   // Set up the newton solver
   Array<ParGridFunction *> pgf_array({&x_u, &x_p});
   RBVMS::SystemResidualMonitor newton_monitor(MPI_COMM_WORLD,
                                               "Newton", 1,
                                               bOffsets);
   NewtonSolver newton_solver(MPI_COMM_WORLD);
   newton_solver.iterative_mode = true;
   newton_solver.SetPrintLevel(-1);
   newton_solver.SetMonitor(newton_monitor);
   newton_solver.SetRelTol(1e-3);
   newton_solver.SetAbsTol(1e-8);
   newton_solver.SetMaxIter(10);
   newton_solver.SetSolver(j_gmres);

   // Define the physical parameters
   LibCoefficient mu(lib_file, "mu", false, mu_param);
   LibVectorCoefficient force(dim, lib_file, "force");

   // Define the stabilisation parameters
   VectorGridFunctionCoefficient adv(&x_u);
   // ElasticInverseEstimateCoefficient invEst(spaces[0]);
   RBVMS::Tau tau(adv, mu);

   // Define evolution
   RBVMS::IncNavStoIntegrator integrator(mu, force, sol, tau);
   RBVMS::ParTimeDepBlockNonlinForm form(spaces, integrator);

   form.SetStrongBC (strong_bdr);
   form.SetWeakBC   (weak_bdr);
   form.SetOutflowBC(outflow_bdr);

   RBVMS::Evolution evo(form, newton_solver);

   // Select the time integrator
   ODESolver *ode_solver = NULL;
   switch (ode_solver_type)
   {
      // Implicit (L-stable) methods
      case 11: ode_solver = new BackwardEulerSolver; break;
      case 12: ode_solver = new SDIRK23Solver(2); break;
      case 13: ode_solver = new SDIRK33Solver; break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 24: ode_solver = new SDIRK34Solver; break;

      case 35: ode_solver = new GeneralizedAlphaSolver(0.5); break;

      default:
         cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         mfem_error("Can not continue");
   }
   ode_solver->Init(evo);

   // 9. Actual time integration
   real_t t = 0.0;
   bool done = false;
   char dimName[] = "xyz";

   std::ostringstream filename;
   filename << "force" << ".dat";

   std::ofstream os;
   if (Mpi::Root()) { os.open(filename.str().c_str()); }
   auto line = [](int len)
   {
      for (int b=0; b<len; ++b) { cout<<"-"; }
      cout<<"\n";
   };

   for (int ti = 0; !done; )
   {
      real_t dt_real = min(dt, t_final - t);
      tau.SetTimeStep(dt_real);
      if (Mpi::Root())
      {
         line(80);
         cout<<"time step: " << ti << ", time: " << t << endl;
         line(80);
      }
      ode_solver->Step(xp, t, dt_real);
      ti++;
      done = (t >= t_final - 1e-8*dt);

      real_t cfl = evo.GetCFL();
      DenseMatrix bdrForce = evo.GetForce();
      if (Mpi::Root())
      {
         // Print to file
         int nbdr = bdrForce.Height();
         os << ti<<"\t"<<t<<"\t"<<cfl<<"\t";
         for (int b=0; b<nbdr; ++b)
            for (int v=0; v<dim; ++v)
            {
               os<<bdrForce(b,v)<<"\t";
            }
         os<<"\n"<< std::flush;

         // Print to screen
         line(80);
         cout<<"cfl = "<<cfl<<endl;
         line(80);

         auto pline = [](int len)
         {
            cout<<" +";
            for (int b=0; b<len; ++b) { cout<<"-"; }
            cout<<"+\n";
         };

         cout<<"\n";
         pline(10+13*nbdr);
         cout<<" | Boundary | ";
         for (int b=0; b<nbdr; ++b)
         {
            cout<<std::setw(10)<<b<<" | ";
         }
         cout<<"\n";
         pline(10+13*nbdr);
         for (int v=0; v<dim; ++v)
         {
            cout<<" | Force "<<dimName[v]<<"  | ";
            for (int b=0; b<nbdr; ++b)
            {
               cout<<std::defaultfloat<<std::setprecision(4)<<std::setw(10);
               cout<<bdrForce(b,v)<<" | ";
            }
            cout<<"\n";
         }
         pline(10+13*nbdr);
         cout<<"\n"<<std::flush;
      }

      if (done || ti % vis_steps == 0)
      {
         if (Mpi::Root())
         {
            line(80);
            cout << "Visit output: Cycle " << ti << "\t Time: " << t << endl;
            line(80);
         }
         x_u.Distribute(xp.GetBlock(0));
         x_p.Distribute(xp.GetBlock(1));

         visit_dc.SetCycle(ti);
         visit_dc.SetTime(t);
         visit_dc.Save();
      }
      x_u.Distribute(xp.GetBlock(0)); // Due to adv???!!!--> other mech!!
   }
   os.close();

   // 10. Free the used memory.
   for (int i = 0; i < fecs.Size(); ++i)
   {
      delete fecs[i];
   }
   for (int i = 0; i < spaces.Size(); ++i)
   {
      delete spaces[i];
   }

   return 0;
}

