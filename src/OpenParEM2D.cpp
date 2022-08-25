////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM2D - A fullwave 2D electromagnetic simulator.                  //
//    Copyright (C) 2022 Brian Young                                          //
//                                                                            //
//    This program is free software: you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Run:          mpirun -np N OpenParEM2D project.proj
//               mpirun -np N --quiet OpenParEM2D project.proj
//
// Description:  This code solves Maxwell's equations in 2D for waveguide/transmission line 
//               eigenvalue/eigenvector pairs for the propagating modes.
//
//               The code is an implementation of the formulation in the paper:
//                   Lee, Jin-Fa, "Finite Element Analysis of Lossy Dielectric Waveguides,"
//                   IEEE Trans. Microwave Theory Techniques, vol. 42, no. 6, June 1994, pp. 1025-1031.
//

#include "mfem.hpp"
#include "modes.hpp"
#include "results.hpp"
#include "convergence.hpp"
#include "frequencyPlan.hpp"
#include "fem2D.hpp"
#include "OpenParEMmaterials.hpp"
#include "project.h"
#include "mesh.hpp"
#include "fieldPoints.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <chrono>
#include <unistd.h>
#include <lapacke.h>

using namespace std;
using namespace mfem;

extern "C" int eigensolve (struct projectData *, int, double, int, PetscInt *, int, PetscInt *, double **, double **, int *, PetscMPIInt);
extern "C" void init_project (struct projectData *);
extern "C" char* get_project_name (const char *);
extern "C" int load_project_file(const char*, projectData*, const char*);
extern "C" void print_project (struct projectData *, struct projectData *, const char *indent);
extern "C" void free_project(projectData*);
extern "C" void matrixTest();
extern "C" void colMajorTest();

void print_copyright_notice ();
void print_license();

void help () {
   PetscPrintf(PETSC_COMM_WORLD,"usage: OpenParEM2D [-h] filename\n");
   PetscPrintf(PETSC_COMM_WORLD,"       -h          : Print this help text\n");
   PetscPrintf(PETSC_COMM_WORLD,"       filename    : Filename of an OpenParEM setup file.\n");
   PetscPrintf(PETSC_COMM_WORLD,"\nOpenParEM2D is a full-wave 2D electromagnetic solver.\n");
}

// check that the field points are within the mesh's bounding box
bool check_field_points (const char *filename, struct projectData *projData, Mesh *mesh, int order)
{
   double tol=1e-12;
   bool fail=false;

   Vector lowerLeft,upperRight;
   mesh->GetBoundingBox(lowerLeft,upperRight,max(order,1));

   int i=0;
   while (i < projData->field_points_count) {

      bool pointFail=false;

      if (projData->field_points_x[i] < lowerLeft.Elem(0)-tol) pointFail=true;
      if (projData->field_points_x[i] > upperRight.Elem(0)+tol) pointFail=true;

      if (projData->field_points_y[i] < lowerLeft.Elem(1)-tol) pointFail=true;
      if (projData->field_points_y[i] > upperRight.Elem(1)+tol) pointFail=true;

      if (! fail && pointFail) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR400: Project file \"%s\":\n",filename);
      }

      if (pointFail) {
         fail=true;
         PetscPrintf(PETSC_COMM_WORLD,"          field.point %g,%g falls outside of the mesh bounding box.\n",
                                      projData->field_points_x[i],projData->field_points_y[i]);
      }

      i++;
   }
   return fail;
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

void printSummary(double frequency, ResultDatabase *resultDatabase, int iteration)
{
   double NpTodB=20*log10(exp(1));
   double eps0=8.8541878176e-12;
   double ko=2*M_PI*frequency*sqrt(4.0e-7*M_PI*eps0);

   PetscPrintf(PETSC_COMM_WORLD,"      Results:\n");
   PetscPrintf(PETSC_COMM_WORLD,"         -------------------------------------------------------------------------------------------------------------------------\n");
   PetscPrintf(PETSC_COMM_WORLD,"          mode     frequency           er,eff      alpha(dB/m)  beta/1000,rad/m          beta/ko        Zii(real)        Zii(imag)\n");
   PetscPrintf(PETSC_COMM_WORLD,"         -------------------------------------------------------------------------------------------------------------------------\n");

   Result *result=resultDatabase->get_Result(frequency,iteration);
   long unsigned int mode=0;
   while (mode < result->get_modeCount()) {
      PetscPrintf(PETSC_COMM_WORLD,"      %5ld",mode+1);
      PetscPrintf(PETSC_COMM_WORLD,"%17.8g",frequency);

      PetscPrintf(PETSC_COMM_WORLD,"%17.8g",pow(imag(result->get_gamma(mode))/ko,2));
      PetscPrintf(PETSC_COMM_WORLD,"%17.8g",real(result->get_gamma(mode))*NpTodB+result->get_alpha_perturbation(mode)/real(result->get_Pz(mode))/2*NpTodB);
      PetscPrintf(PETSC_COMM_WORLD,"%17.8g",imag(result->get_gamma(mode))/1000);
      PetscPrintf(PETSC_COMM_WORLD,"%17.8g",imag(result->get_gamma(mode))/ko);

      long unsigned int offset=0;
      if (result->get_modalImpedanceCalculation()) offset=mode;
      else offset=mode+mode*result->get_modeCount();

      if (result->get_Z(offset) != complex<double>(DBL_MAX,DBL_MAX)) {
         PetscPrintf(PETSC_COMM_WORLD,"%17.8g",real(result->get_Z(offset)));
         PetscPrintf(PETSC_COMM_WORLD,"%17.8g",imag(result->get_Z(offset)));
      }

      PetscPrintf(PETSC_COMM_WORLD,"\n");

      mode++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"         -------------------------------------------------------------------------------------------------------------------------\n");

   return;
}

void show_memory (int show, string space)
{
   double vm_usage, resident_set;

   if (show) {
      process_mem_usage(vm_usage, resident_set);
      PetscPrintf(PETSC_COMM_WORLD,"%smemory (kB): %g  page size (kB): %g\n",space.c_str(),vm_usage,resident_set);
   }
}

void exit_job_on_error (PetscMPIInt rank, chrono::system_clock::time_point job_start_time, const char *lockfile)
{
   PetscPrintf(PETSC_COMM_WORLD,"Job Complete\n");

   chrono::system_clock::time_point job_end_time=chrono::system_clock::now();
   chrono::duration<double> elapsed = job_end_time - job_start_time;
   PetscPrintf(PETSC_COMM_WORLD,"Elapsed time: %g s\n",elapsed.count());

   // remove the lock - not 100% safe
   if (rank == 0) {
      if (std::filesystem::exists(lockfile)) {
         std::filesystem::remove(lockfile);
      }
   }

   MPI_Abort(PETSC_COMM_WORLD,1);
   //PetscFinalize();
   SlepcFinalize();
   exit (1);
}

int main(int argc, char *argv[])
{
   int printHelp=0;
   struct projectData projData,defaultData;
   BoundaryDatabase boundaryDatabase;
   BorderDatabase borderDatabase;
   ResultDatabase resultDatabase;
   ConvergenceDatabase *convergenceDatabase;
   FieldPointDatabase fieldPointDatabase;
   FrequencyPlan frequencyPlan;
   double *alphaList,*betaList;
   int iteration;
   double frequency,lastFrequency,betaScale;
   bool refineMesh,restartMesh;
   int matrixSize;
   complex<double> Pz,Zo,ZoPV,ZoPI,ZoVI;
   //double loss_adder;
   int use_initial_guess;   // integer since this is going to the C-language eigensolve
   const char *projFile;
   //const char *device_config = "cpu";
   //bool herm_conv = true;
   meshMaterialList meshMaterials;
   MaterialDatabase localMaterialDatabase;
   MaterialDatabase materialDatabase;
   PetscMPIInt size,rank;
   int is_locked=0;

   // Initialize Petsc and MPI
   SlepcInitializeNoArguments();
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   chrono::system_clock::time_point job_start_time=chrono::system_clock::now();

   // parse inputs
   if (argc <= 1) printHelp=1;
   else {
      if (strcmp(argv[1],"-h") == 0) printHelp=1;
      else projFile=argv[1];
   }
   if (printHelp) {help(); /*PetscFinalize()*/ SlepcFinalize();; exit(1);}

   print_copyright_notice ();

   // create a lock file - this method is not guaranteed
   char *rootProjName=get_project_name (projFile);
   stringstream ssLock;
   ssLock << "." << rootProjName << ".lock";

   if (rank == 0) {
      if (std::filesystem::exists(ssLock.str().c_str())) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR401: Project \"%s\" is locked.\n",projFile);
         MPI_Abort(PETSC_COMM_WORLD,1);
         //PetscFinalize();
         SlepcFinalize();
         exit(1);
      }

      // assume that the file does not exist

      ofstream lock;
      lock.open(ssLock.str().c_str(),ofstream::out);
      if (lock.is_open()) {
         lock << "locked" << endl;
         lock.close();
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR402: Cannot open \"%s\" for writing.\n",ssLock.str().c_str());
         //PetscFinalize();
         SlepcFinalize();
         exit(1);
      }

      is_locked=1;
   }
   MPI_Bcast(&is_locked,1,MPI_INT,0,PETSC_COMM_WORLD);

   PetscPrintf(PETSC_COMM_WORLD,"Setting up ...\n");

   // load the project file
   init_project (&defaultData);
   init_project (&projData);
   PetscPrintf(PETSC_COMM_WORLD,"   loading project file \"%s\"\n",projFile);
   if (load_project_file (projFile, &projData, "   ")) {
      if (projData.debug_show_project) {print_project (&projData,&defaultData,"      ");}
      exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   }
   if (projData.debug_show_project) {print_project (&projData,&defaultData,"      ");}

   if (projData.output_show_license) {
      print_license();
      exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   }

   // load materials
   bool fail=false;
   if (strlen(projData.materials_global_name) == 0) {
      if (strlen(projData.materials_local_name) == 0) {
         PetscPrintf(PETSC_COMM_WORLD,"   ERROR403: No material databases are defined.\n"); fail=true;
      } else {
         if (materialDatabase.load(projData.materials_local_path,projData.materials_local_name, projData.materials_check_limits)) fail=true;
      }
   } else {
      if (materialDatabase.load(projData.materials_global_path,projData.materials_global_name, projData.materials_check_limits)) {
         fail=true;

         // go ahead and try the local to get error messages
         if (strlen(projData.materials_local_name) != 0) {
            if (localMaterialDatabase.load(projData.materials_local_path,projData.materials_local_name, projData.materials_check_limits)) fail=true;
         }
      } else {
         if (strlen(projData.materials_local_name) != 0) {
            if (localMaterialDatabase.load(projData.materials_local_path,projData.materials_local_name, projData.materials_check_limits)) {
               fail=true;
            } else {
               if (materialDatabase.merge(&localMaterialDatabase,materialDatabase.get_indent())) fail=true;
            }
         }
      }
   }
   if (fail) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   if (projData.debug_show_materials) {materialDatabase.print("   ");}

   show_memory (projData.debug_show_memory, "   ");

   // get the boundary database indicating impedance integration paths and boundary conditions
   if (boundaryDatabase.load(projData.mode_definition_file,projData.solution_check_closed_loop)) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   if (projData.debug_show_mode_definitions) boundaryDatabase.print();

   double eps0=8.8541878176e-12;

   // various files

   stringstream ssTemp;
   ssTemp << "temp_" << projData.project_name;

   stringstream ssParaView;
   ssParaView << "ParaView_" << projData.project_name;

   stringstream ssResults;
   ssResults << projData.project_name << "_results.csv";

   stringstream ssFields;
   ssFields << projData.project_name << "_fields.csv";

   stringstream ssTests;
   ssTests << projData.project_name << "_prototype_test_cases.csv";

   // file removals to prevent stale data

   if (rank == 0) {
      if (std::filesystem::exists(ssTemp.str().c_str())) {
        std::filesystem::remove_all(ssTemp.str().c_str());
      }

      if (std::filesystem::exists(ssResults.str().c_str())) {
        std::filesystem::remove(ssResults.str().c_str());
      }

      if (std::filesystem::exists(ssFields.str().c_str())) {
        std::filesystem::remove(ssFields.str().c_str());
      }

      if (std::filesystem::exists(ssTests.str().c_str())) {
        std::filesystem::remove(ssTests.str().c_str());
      }

      if (std::filesystem::exists(ssParaView.str().c_str())) {
        std::filesystem::remove_all(ssParaView.str().c_str());
      }

      // create the temp directory
      if (! std::filesystem::create_directory(ssTemp.str().c_str())) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR404: Failed to create results directory \"%s\".\n",ssTemp.str().c_str());
         exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
      }
   }
   MPI_Barrier(PETSC_COMM_WORLD);

   // Enable hardware devices such as GPUs, and programming models such as
   // CUDA, OCCA, RAJA and OpenMP based on command line options.
   //Device device(device_config);
   //if (myid == 0) { device.Print(); }

   PetscPrintf(PETSC_COMM_WORLD,"Loading mesh and assigning materials ...\n");
   if (!projData.materials_check_limits) {
      PetscPrintf(PETSC_COMM_WORLD,"    Skipping limit checks on material values\n");
   }

   // get the mapping from the region number to the material name
   if (meshMaterials.loadGMSH(projData.mesh_file)) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   //meshMaterials.print();

   // read in the mesh to a serial mesh structure
   Mesh mesh(projData.mesh_file, 1, 1);
   //   mesh.ScaleElements(0.001);    // hard coded for now to convert from mm to m.  - ToDo - Generalize, but ScaleElements may be broken.
   int dim=mesh.Dimension();

   int i=0;
   while (i < projData.mesh_uniform_refinement_count) {
      mesh.UniformRefinement();
      i++;
   }

   // check field points paths for scale
   if (check_field_points (projFile, &projData, &mesh, projData.mesh_order)) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   if (boundaryDatabase.check_scale(&mesh,projData.mesh_order)) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());

   // mark the boundaries
   boundaryDatabase.mark_boundaries (&mesh, &borderDatabase);
   mesh.SetAttributes(); // recalulates the support data structures

   // 2D solver, so require dim=2
   if (! (dim == 2)) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR405: Mesh must be 2-dimensional.\n");
      exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   }

   // construct coefficient vectors to hold material properties per FEM element
   Vector ko2RePermittivity(mesh.attributes.Max());
   Vector ko2ImPermittivity(mesh.attributes.Max());
   Vector InvPermeability(mesh.attributes.Max());  // inverse of the permeability for use in (12)
   Vector InvOmegaMu(mesh.attributes.Max());
   Vector OmegaMu(mesh.attributes.Max());
   Vector InvOmegaMuEps(mesh.attributes.Max());

   // assign materials

   if (mesh.attributes.Max() != meshMaterials.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR406: Mesh file \"%s\" does not include the correct number of regions for material definitions.\n",projData.mesh_file);
      PetscPrintf(PETSC_COMM_WORLD,"       The $PhysicalNames block should have %d entries, but only %d were found\n.",
                                           mesh.attributes.Max()+1,meshMaterials.size()+1);
      exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   }

   alphaList=(double *)malloc (projData.solution_modes*sizeof(double));
   betaList=(double *)malloc (projData.solution_modes*sizeof(double));

   // loop for all frequencies

   show_memory (projData.debug_show_memory, "");

   // set up the frequency plan
   if (frequencyPlan.assemble(&projData)) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
   if (projData.debug_show_frequency_plan) frequencyPlan.print();
   lastFrequency=0;  // for scaling the initial guess

   // set up pmesh when not using adaptive mesh refinement
   ParMesh *pmesh=NULL;
   if (! frequencyPlan.is_refining()) {
      pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
      pmesh->ReorientTetMesh();    // re-orient the mesh in case of a tet mesh
   }

   // cannot use initial guess on the first iteration since there is no data
   use_initial_guess=0;

   // loop
   while (frequencyPlan.get_frequency(projData.refinement_frequency,&frequency,&refineMesh,&restartMesh)) {
      PetscPrintf(PETSC_COMM_WORLD,"Frequency: %g\n",frequency);

      double ko=2*M_PI*frequency*sqrt(4.0e-7*M_PI*eps0);
      double ko2=pow(ko,2);

      // for initial guess scaling with frequency
      if (lastFrequency == 0) lastFrequency=frequency;
      betaScale=frequency/lastFrequency;

      // set up pmesh when using adaptive mesh refinement
      if (refineMesh && restartMesh) {
         if (pmesh != NULL) delete pmesh;
         pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
         pmesh->ReorientTetMesh();

         // new mesh means that the existing solution is invalid
         use_initial_guess=0;
      }

      if (refineMesh) {
         if (restartMesh) {PetscPrintf(PETSC_COMM_WORLD,"   Refining mesh starting with the initial mesh.\n");}
         else {PetscPrintf(PETSC_COMM_WORLD,"   Refining mesh starting with the last mesh.\n");}
      }

      int j=0;
      while (j < pmesh->attributes.Max()) {
         Material *useMaterial=materialDatabase.get(meshMaterials.get_name(meshMaterials.get_index(j)));

         if (useMaterial != NULL) {

            // permittivity
            complex<double> e=useMaterial->get_eps(projData.solution_temperature,frequency,materialDatabase.get_tol(),materialDatabase.get_indent());
            if (e == complex<double>(-DBL_MAX,0)) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());

            ko2RePermittivity[j]=ko2*real(e)/eps0;
            ko2ImPermittivity[j]=ko2*imag(e)/eps0;

            // permeability
            double mu=useMaterial->get_mu(projData.solution_temperature,frequency,materialDatabase.get_tol(),materialDatabase.get_indent());
            if (mu == -DBL_MAX) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());

            InvPermeability[j]=1.0/(mu/(4e-7*M_PI));    // 1/relative permeability
            InvOmegaMu[j]=1.0/(2*M_PI*frequency*mu);    // 1/(w*permeability)
            OmegaMu[j]=2*M_PI*frequency*mu;    // w*permeability

            InvOmegaMuEps[j]=1/(2*M_PI*frequency*mu*real(e));

         } else {
            PetscPrintf(PETSC_COMM_WORLD,"ERROR407: Material \"%s\" for region %d is not present in the material database.\n",
                                         meshMaterials.get_name(j).c_str(),j+1);
            exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
         }

         j++;
      }

      PWConstCoefficient ko2Re_e(ko2RePermittivity);
      PWConstCoefficient ko2Im_e(ko2ImPermittivity);
      PWConstCoefficient Inv_mu(InvPermeability);
      PWConstCoefficient Inv_w_mu(InvOmegaMu);
      PWConstCoefficient w_mu(OmegaMu);
      PWConstCoefficient Inv_w_mu_e(InvOmegaMuEps);

      show_memory (projData.debug_show_memory, "   ");

      // set up for convergence testing
      convergenceDatabase=new ConvergenceDatabase();
      convergenceDatabase->initialize(projData.solution_modes,projData.refinement_required_passes,projData.refinement_variable,projData.refinement_tolerance);

      iteration=0;
      bool iterate=true;
      while (iterate) {
         chrono::system_clock::time_point solve_start_time=chrono::system_clock::now();

         if (refineMesh) {PetscPrintf(PETSC_COMM_WORLD,"   Iteration %d ...\n",iteration+1);}
         else {
            PetscPrintf(PETSC_COMM_WORLD,"   Using existing mesh ...\n");
            iterate=false;
         }
         PetscPrintf(PETSC_COMM_WORLD,"      mesh size: %d\n",GetGlobalNE(pmesh));
         show_memory (projData.debug_show_memory, "      ");

         double h_min,h_max,kappa_min,kappa_max;
         pmesh->GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
         PetscPrintf(PETSC_COMM_WORLD,"      mesh worst element aspect ratio: %g",kappa_max);
         if (kappa_max < 5) PetscPrintf(PETSC_COMM_WORLD," < target: 5\n");
         else PetscPrintf(PETSC_COMM_WORLD," > target: 5\n");

         // create the finite element spaces and matrices on the parallel mesh
         fem2D *fem=new fem2D(&projData,pmesh,projData.mesh_order,frequency,iteration,&ko2Re_e,&ko2Im_e,&Inv_mu,&w_mu,ssTemp.str());

         // scale beta to approximate the change due to a shift in frequency
         if (use_initial_guess && projData.solution_active_mode_count > 0) {
            betaList[0]*=betaScale*2;  // ToDo: Find a way to eliminate the scale factor of 2x
            betaScale=1;   // will update when the frequency changes
         }

         // solve the eigenvalue problem, Ax=kBx
         PetscInt *ess_tdof_ND=fem->get_ess_tdof_ND(&boundaryDatabase,&borderDatabase);
         PetscInt *ess_tdof_H1=fem->get_ess_tdof_H1(&boundaryDatabase,&borderDatabase);
         if (projData.debug_skip_solve || eigensolve (&projData,use_initial_guess,frequency,
                                                      fem->get_ess_tdof_size_ND(),ess_tdof_ND,
                                                      fem->get_ess_tdof_size_H1(),ess_tdof_H1,
                                                      &alphaList,&betaList,&matrixSize,rank) != 0) 
         {
            // failed or canceled - stop iterating
            iterate=false;
         } else {

            // process the eigenvalue solutions
            fem->set_t_size();
            fem->set_z_size();
            if (fem->buildFields(alphaList,betaList)) exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
            fem->sort();
            fem->calculateImpedanceMatrix(projData.solution_impedance_definition,&boundaryDatabase,&borderDatabase);
            fem->calculatePerturbationalLoss(&boundaryDatabase,&borderDatabase,&materialDatabase);
            fem->calculateFieldPoints(&fieldPointDatabase);
            fem->saveParaView();
            chrono::system_clock::time_point solve_end_time=chrono::system_clock::now();
            Result *result=fem->updateResults(&resultDatabase,convergenceDatabase,solve_start_time,solve_end_time);

            chrono::duration<double> elapsed = solve_end_time - solve_start_time;
            PetscPrintf(PETSC_COMM_WORLD,"         solve elapsed time: %g s\n",elapsed.count());

            // stop if no solutions were found
            if (projData.solution_active_mode_count == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"ERROR408: Cannot continue.  Modify setup to increase chances of finding a solution.\n");
               PetscPrintf(PETSC_COMM_WORLD,"          - check frequency\n");
               PetscPrintf(PETSC_COMM_WORLD,"          - check dimensions\n");
               PetscPrintf(PETSC_COMM_WORLD,"          - reduce solution.tolerance\n");
               PetscPrintf(PETSC_COMM_WORLD,"          - decrease mesh.order\b");
               PetscPrintf(PETSC_COMM_WORLD,"          - apply uniform mesh refinement by setting mesh.uniform_refinement.count\n");
               exit_job_on_error (rank,job_start_time,ssLock.str().c_str());
            }

            // print a summary to the console
            printSummary(frequency, &resultDatabase, iteration);

            // show convergence
            if (iteration > 0) convergenceDatabase->show_progress();

            if (iteration < projData.refinement_iteration_min-1) convergenceDatabase->set_not_converged();

            if (refineMesh && convergenceDatabase->is_converged()) iterate=false;

            // refine the mesh using all modes
            if (refineMesh && ! convergenceDatabase->is_converged() && projData.solution_active_mode_count > 0) {
               fem->ZZrefineMeshes(convergenceDatabase);
            }

            // write the fields as the initial guess for the next calculation
            if (projData.solution_use_initial_guess) fem->writeInitialGuess();

            PetscPrintf(PETSC_COMM_WORLD,"      Finished\n");

            chrono::system_clock::time_point iteration_end_time=chrono::system_clock::now();
            elapsed = iteration_end_time - job_start_time;
            PetscPrintf(PETSC_COMM_WORLD,"      cummulative elapsed time: %g s\n",elapsed.count());

            result->set_end_refine_time(iteration_end_time);
         }

         if (ess_tdof_ND) {PetscFree(ess_tdof_ND); ess_tdof_ND=nullptr;}
         if (ess_tdof_H1) {PetscFree(ess_tdof_H1); ess_tdof_H1=nullptr;}

         // update stats with convergence and information
         resultDatabase.update_convergence(frequency,convergenceDatabase->is_converged(),convergenceDatabase->get_last_error());

         // save the results to a results csv file
         if (rank == 0) resultDatabase.save(ssResults.str().c_str());

         // initial guess
         if (projData.solution_use_initial_guess) use_initial_guess=true;

         delete fem;
         ++iteration;

         if (iteration == projData.refinement_iteration_max) iterate=false;
      }

      if (refineMesh && convergenceDatabase->is_converged()) PetscPrintf(PETSC_COMM_WORLD,"   Converged\n");
      if (refineMesh && ! convergenceDatabase->is_converged()) {PetscPrintf(PETSC_COMM_WORLD,"   NOT CONVERGED\n");}

      delete convergenceDatabase;
      lastFrequency=frequency;
   }

   // save the results and field point data as test cases
   if (projData.test_create_cases && rank == 0) {
      fieldPointDatabase.normalize();
      fieldPointDatabase.save(ssFields.str().c_str());

      resultDatabase.save_as_test(ssTests.str().c_str(),projData.project_name);
      fieldPointDatabase.save_as_test(ssTests.str().c_str(),projData.project_name);
   }

   PetscPrintf(PETSC_COMM_WORLD,"Job Complete\n");

   MPI_Barrier(PETSC_COMM_WORLD);
   if (rank == 0) {
      if (! projData.debug_tempfiles_keep && std::filesystem::exists(ssTemp.str().c_str())) {
        std::filesystem::remove_all(ssTemp.str().c_str());
      }
   }

   delete pmesh;
   free_project (&projData);
   if (alphaList) free(alphaList);
   if (betaList) free(betaList);

   chrono::system_clock::time_point job_end_time=chrono::system_clock::now();
   chrono::duration<double> elapsed = job_end_time - job_start_time;
   PetscPrintf(PETSC_COMM_WORLD,"Elapsed time: %g s\n",elapsed.count());

   show_memory (projData.debug_show_memory, "");

   // remove the lock - not 100% safe
   if (rank == 0) {
      if (std::filesystem::exists(ssLock.str().c_str())) {
        std::filesystem::remove(ssLock.str().c_str());
      }
   }

   //PetscFinalize();
   SlepcFinalize();

   return 0;
}

// last ERROR useds is 408
