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
   if (Mpi::Root())
   {
      printInfo();
   }

   // 3. Parse command-line options.
   const char *mesh_file = "../../mfem/data/inline-quad.mesh";
   const char *ref_file  = "";
   int order = 1;
   int ref_levels = 0;

   real_t mu_param = 1.0;
   real_t penalty = -1.0;

   const char *lib_file = "libfun.so";

   int ode_solver_type = 11;
   real_t dt = 0.01;
   real_t t_final = 10.0;
   int vis_steps = 1;

   OptionsParser args(argc, argv);

   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order isoparametric space.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh.");

   args.AddOption(&lib_file, "-l", "--lib",
                  "Library file to use for IC, BC, force and solution function defintions.");

   args.AddOption(&mu_param, "-m", "--mu",
                  "Sets the diffusion parameters, should be positive.");

   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "...");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");

   args.Parse();
   if (!args.Good())
   {
      if (myid == 0) { args.PrintUsage(cout); }
      return 1;
   }
   if (myid == 0) { args.PrintOptions(cout); }

   // 4. Read the mesh from the given mesh file.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // Refine
   {
      if (mesh.NURBSext && (strlen(ref_file) != 0))
      {
         mesh.RefineNURBSFromFile(ref_file);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
      if (myid == 0) { mesh.PrintInfo(); }
   }

   // Partition
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

      if (myid == 0)
      {
         mfem::out << "Number of finite element unknowns:\n";
         mfem::out << "\tVelocity = "<<spaces[0]->GlobalTrueVSize() << endl;
         mfem::out << "\tPressure = "<<spaces[1]->GlobalTrueVSize() << endl;
         mfem::out << "Number of finite element unknowns per partition:\n";
         mfem::out <<  "\tVelocity = "; udof.Print(mfem::out, num_procs);
         mfem::out <<  "\tPressure = "; pdof.Print(mfem::out, num_procs);
      }
   }

   // Mark all velocity boundary dofs as essential
   Array<Array<int> *> ess_bdr(2);
   Array<int> ess_bdr_u(spaces[0]->GetMesh()->bdr_attributes.Max());
   Array<int> ess_bdr_p(spaces[1]->GetMesh()->bdr_attributes.Max());

   ess_bdr_p = 0;
   ess_bdr_u = 1;

   ess_bdr[0] = &ess_bdr_u;
   ess_bdr[1] = &ess_bdr_p;

   // 6. Define the solution vector xp as a finite element grid function
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

   // 7. Define the time stepping algorithm
   // Set up the preconditioner
   Array<Solver *> sol_array({new HypreSmoother(),
            new HypreSmoother()});
   RBVMS::JacobianPreconditioner jac_prec(bOffsets, sol_array);

   // Set up the Jacobian solver
   RBVMS::GeneralResidualMonitor j_monitor(MPI_COMM_WORLD,"\t\t\t\tFGMRES", 25);
   FGMRESSolver j_gmres(MPI_COMM_WORLD);
   j_gmres.iterative_mode = false;
   j_gmres.SetRelTol(1e-2);
   j_gmres.SetAbsTol(1e-12);
   j_gmres.SetMaxIter(300);
   j_gmres.SetPrintLevel(-1);
   j_gmres.SetMonitor(j_monitor);
   j_gmres.SetPreconditioner(jac_prec);

   // Set up the newton solver
   Array<ParGridFunction *> pgf_array({&x_u, &x_p});
   RBVMS::SystemResidualMonitor newton_monitor(MPI_COMM_WORLD,"Newton", 1,
                                               bOffsets, nullptr, &xp,
                                               pgf_array);
   NewtonSolver newton_solver(MPI_COMM_WORLD);
   newton_solver.iterative_mode = true;
   newton_solver.SetPrintLevel(-1);
   newton_solver.SetMonitor(newton_monitor);
   newton_solver.SetRelTol(1e-4);
   newton_solver.SetAbsTol(1e-8);
   newton_solver.SetMaxIter(25);
   newton_solver.SetSolver(j_gmres);

   // 7. Define the weka form
   // Define the physical parameters
   LibCoefficient mu(lib_file, "mu", false, mu_param);
   LibVectorCoefficient force(dim, lib_file, "force");

   // Define the stabilisation parameters
   VectorGridFunctionCoefficient adv(&x_u);
   // ElasticInverseEstimateCoefficient invEst(spaces[0]);
   RBVMS::Tau tau(adv, mu);

   // Define evolution
   RBVMS::IncNavStoForm form(spaces,
                             mu, force,
                             tau, tau);
   Array<Vector *> rhs(2);
   rhs = nullptr; // Set all entries in the array
   form.SetEssentialBC(ess_bdr, rhs);

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

      default:
         cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         mfem_error("Can not continue");
   }
   cout << "ODE solver type: " << ode_solver_type << '\n';
   ode_solver->Init(evo);

   // 9. Actual time integration
   real_t t = 0.0;
   bool done = false;
   for (int ti = 0; !done; )
   {
      cout<<"----------------------------------------\n";
      cout << "time step: " << ti << ", time: " << t << endl;
      cout<<"----------------------------------------\n";

      real_t dt_real = min(dt, t_final - t);
      tau.SetTimeStep(dt_real);
      ode_solver->Step(xp, t, dt_real);
      ti++;
      done = (t >= t_final - 1e-8*dt);

      if (done || ti % vis_steps == 0)
      {
         cout<<"\n\n";
         cout << "Visit output: Cycle " << ti << "\t Time: " << t << endl;
         cout<<"\n\n";

         x_u.Distribute(xp.GetBlock(0));
         x_p.Distribute(xp.GetBlock(1));

         visit_dc.SetCycle(ti);
         visit_dc.SetTime(t);
         visit_dc.Save();
      }
   }

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

