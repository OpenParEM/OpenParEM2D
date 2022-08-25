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

#include "fem2D.hpp"

void saveMat (ParBilinearForm &blf, string name, string directory, char *uniquifier, bool transpose) {
   stringstream ss;
   ss << directory << "/" << name << "." << uniquifier;

   HypreParMatrix *mat=blf.ParallelAssemble();
   if (transpose) mat=mat->Transpose();
   mat->Print(ss.str().c_str());

   delete mat;
}

void saveMat (ParMixedBilinearForm &blf, string name, string directory, char *uniquifier, bool transpose) {
   stringstream ss;
   ss << directory << "/" << name << "." << uniquifier;

   HypreParMatrix *mat=blf.ParallelAssemble();
   if (transpose) mat=mat->Transpose();
   mat->Print(ss.str().c_str());

   delete mat;
}

int GetGlobalNE (ParMesh *pmesh)
{
   int localNE=pmesh->GetNE();
   int globalNE=0;

   MPI_Allreduce (&localNE,&globalNE,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

   return globalNE;
}

void printImpedanceVector (complex<double> *z, long unsigned int n)
{
   long unsigned int i=0;
   while (i < n) {
      PetscPrintf(PETSC_COMM_WORLD,"               %3ld",i+1);
      if (z[i] != complex<double>(DBL_MAX,DBL_MAX)) PetscPrintf(PETSC_COMM_WORLD," (%g,%g)\n",real(z[i]),imag(z[i]));
      else PetscPrintf(PETSC_COMM_WORLD," not defined\n");
      i++;
   }
}

//--------------------------------------------------------------------------------------------------------------------
// Fields
//--------------------------------------------------------------------------------------------------------------------

Fields::~Fields()
{
   if (grid_Et_re) {delete grid_Et_re; grid_Et_re=nullptr;}
   if (grid_Et_im) {delete grid_Et_im; grid_Et_im=nullptr;}
   if (grid_Ez_re) {delete grid_Ez_re; grid_Ez_re=nullptr;}
   if (grid_Ez_im) {delete grid_Ez_im; grid_Ez_im=nullptr;}

   if (grid_Ht_re) {delete grid_Ht_re; grid_Ht_re=nullptr;}
   if (grid_Ht_im) {delete grid_Ht_im; grid_Ht_im=nullptr;}
   if (grid_Hz_re) {delete grid_Hz_re; grid_Hz_re=nullptr;}
   if (grid_Hz_im) {delete grid_Hz_im; grid_Hz_im=nullptr;}

   if (grid_Pz_re) {delete grid_Pz_re; grid_Pz_re=nullptr;}
   if (grid_Pz_im) {delete grid_Pz_im; grid_Pz_im=nullptr;}

   // owned by above?  Deleting causes double delete error in Valgrind
   //if (EfieldRe) {delete EfieldRe; EfieldRe=nullptr;}
   //if (EfieldIm) {delete EfieldIm; EfieldIm=nullptr;}
   //if (HfieldRe) {delete HfieldRe; HfieldRe=nullptr;}
   //if (HfieldIm) {delete HfieldIm; HfieldIm=nullptr;}
}

// field = false for Efield
// field = true  for Hfield
// allocates eVecRe and eVecIm, and these need to be freed elsewhere
int Fields::loadEigenvector (fem2D *fem, long unsigned int mode, bool field, double **eVecRe, double **eVecIm)
{
   char filename[128];
   string line;
   int lineNumber=0;
   vector<string> tokens;
   bool foundVec=false,foundEndVec=false;
   size_t allocated, current, blockSize=256;
   size_t length;

   if (field) { // Hfield
      sprintf (filename,"temp_%s/Hfield_mode_%ld.dat",fem->projData->project_name,mode+1);
   } else {     // Efield
      sprintf (filename,"temp_%s/Efield_mode_%ld.dat",fem->projData->project_name,mode+1);
   }

   ifstream eVecFile;
   eVecFile.open(filename,ifstream::in);

   if (eVecFile.is_open()) {

      allocated=blockSize;
      *eVecRe=(double *)malloc(allocated*sizeof(double));
      *eVecIm=(double *)malloc(allocated*sizeof(double));
      current=0;

      if (*eVecRe == nullptr || *eVecIm == nullptr) {
         if (eVecRe != nullptr) free(eVecRe);
         if (eVecIm != nullptr) free(eVecIm);
         PetscPrintf(PETSC_COMM_WORLD,"ERROR100: Failed to allocate memory.\n");
         return 1;
      }

      while (getline(eVecFile,line)) {
         lineNumber++;
         split_on_space (&tokens,line);

         if (foundVec) {
            if (tokens.size() > 0) if (tokens[0].compare("];") == 0) foundEndVec=true;

            if (! foundEndVec) {
               bool savedValue=false;

               if (tokens.size() == 1) { // real part only
                  try {
                     (*eVecRe)[current]=stod(tokens[0]);
                  } catch(const std::out_of_range&) {
                     // assume that the number is less than DBL_MIN and set to 0
                     (*eVecRe)[current]=0;
                  }
                  (*eVecIm)[current]=0;

                  savedValue=true;
                  current++;
               } else if (tokens.size() == 2) {  // ?
                 PetscPrintf(PETSC_COMM_WORLD,"ERROR107: Unsupported formatting in file \"%s\" at line %d\n",filename,lineNumber);
                 return 1; 
               } else if (tokens.size() == 3) {  // complex

                  try {
                     (*eVecRe)[current]=stod(tokens[0]);
                  } catch(const std::out_of_range&) {
                     // assume that the number is less than DBL_MIN and set to 0
                     (*eVecRe)[current]=0;
                  }
                  length=tokens[2].length();

                  try {
                     (*eVecIm)[current]=stod(tokens[2].substr(0,length-1));
                  } catch(const std::out_of_range&) {
                     // assume that the number is less than DBL_MIN and set to 0
                     (*eVecIm)[current]=0;
                  }
                  if (tokens[1].compare("-") == 0) (*eVecIm)[current]=-(*eVecIm)[current];

                  savedValue=true;
                  current++;
               } else {
                 PetscPrintf(PETSC_COMM_WORLD,"ERROR108: Unsupported formatting in file \"%s\" at line %d\n",filename,lineNumber);
               }

               if (savedValue) {
                  if (current == allocated) {
                     allocated+=blockSize;
                     *eVecRe=(double *)realloc(*eVecRe,allocated*sizeof(double));
                     *eVecIm=(double *)realloc(*eVecIm,allocated*sizeof(double));

                     if (*eVecRe == nullptr || *eVecIm == nullptr) {
                        if (eVecRe != nullptr) free(eVecRe);
                        if (eVecIm != nullptr) free(eVecIm);
                        PetscPrintf(PETSC_COMM_WORLD,"ERROR101: Failed to allocate memory.\n");
                        return 1;
                     }
                  }
               }

            }

         } else {
            if (tokens.size() > 0) {if (tokens[0].substr(0,3).compare("Vec") == 0) foundVec=true;}
         }
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR109: Unable to open file \"%s\" for reading.\n",filename);
      return 1;
   }

   // check that the expected number of lines loaded
   if (current != fem->t_size+fem->z_size) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR110: Unexpected line load count of %ld vs. expected %ld in file \"%s\".\n",current,fem->t_size+fem->z_size,filename);
      return 1;
   }

   // scale the eigenvector by the largest Et vector component
   // This is equivalent to scaling the amplitude and phase-shifting the solution.
   // Makes for nicer plots so that the real part dominates the imaginary part.

   if (current > 0) {
      double x1,x2,y1,y2;
      double mag2=0;
      double largest=0;
      int index=0;

      // get the scale
      x2=0;
      y2=0;
      if (field) {
         // must run for E-fields before running for H-fields
         x2=xScale;
         y2=yScale;
      } else {
         // largest Et component
         size_t i=0;
         while (i < fem->t_size) {
            mag2=(*eVecRe)[i]*(*eVecRe)[i]+(*eVecIm)[i]*(*eVecIm)[i];
            if (mag2 > largest) {largest=mag2; index=i;}
            i++;
         }
         x2=(*eVecRe)[index];
         y2=(*eVecIm)[index];

         // save the scale factor for use with H fields
         xScale=x2;
         yScale=y2;
      }

      // scale all
      if (fabs(x2) > 0 || fabs(y2) > 0) {
         size_t i=0;
         while (i < current) {
            x1=(*eVecRe)[i];
            y1=(*eVecIm)[i];

            // (x1+jy1)/(x2+jy2)
            (*eVecRe)[i]=(x1*x2+y1*y2)/(x2*x2+y2*y2);
            (*eVecIm)[i]=(x2*y1-x1*y2)/(x2*x2+y2*y2);

            i++;
         }
      }
   }

   return 0;
}

bool Fields::build (fem2D *fem, long unsigned int mode)
{
   bool HFIELD=true;
   bool EFIELD=false;

   HYPRE_BigInt *offset_ND=fem->fespace_ND->GetTrueDofOffsets();
   HYPRE_BigInt *offset_H1=fem->fespace_H1->GetTrueDofOffsets();

   // must load E before H to get the scale factor

   // build the Efield
   if (! loadEigenvector (fem,mode,EFIELD,&EfieldRe,&EfieldIm)) {

      Vector EtRe=Vector(EfieldRe+offset_ND[0],offset_ND[1]-offset_ND[0]);
      Vector EtIm=Vector(EfieldIm+offset_ND[0],offset_ND[1]-offset_ND[0]);
      Vector EzRe=Vector(EfieldRe+fem->t_size+offset_H1[0],offset_H1[1]-offset_H1[0]);
      Vector EzIm=Vector(EfieldIm+fem->t_size+offset_H1[0],offset_H1[1]-offset_H1[0]);

      grid_Et_re=new ParGridFunction(fem->fespace_ND);
      grid_Et_im=new ParGridFunction(fem->fespace_ND);
      grid_Ez_re=new ParGridFunction(fem->fespace_H1);
      grid_Ez_im=new ParGridFunction(fem->fespace_H1);

      grid_Et_re->Distribute(EtRe);
      grid_Et_im->Distribute(EtIm);
      grid_Ez_re->Distribute(EzRe);
      grid_Ez_im->Distribute(EzIm);

   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR102: Failed to load Efield.\n");
      return true;
   }

   // build the Hfield
   if (! loadEigenvector (fem,mode,HFIELD,&HfieldRe,&HfieldIm)) {

      Vector HtRe=Vector(HfieldRe+offset_ND[0],offset_ND[1]-offset_ND[0]);
      Vector HtIm=Vector(HfieldIm+offset_ND[0],offset_ND[1]-offset_ND[0]);
      Vector HzRe=Vector(HfieldRe+fem->t_size+offset_H1[0],offset_H1[1]-offset_H1[0]);
      Vector HzIm=Vector(HfieldIm+fem->t_size+offset_H1[0],offset_H1[1]-offset_H1[0]);

      grid_Ht_re=new ParGridFunction(fem->fespace_ND);
      grid_Ht_im=new ParGridFunction(fem->fespace_ND);
      grid_Hz_re=new ParGridFunction(fem->fespace_H1);
      grid_Hz_im=new ParGridFunction(fem->fespace_H1);

      grid_Ht_re->Distribute(HtRe);
      grid_Ht_im->Distribute(HtIm);
      grid_Hz_re->Distribute(HzRe);
      grid_Hz_im->Distribute(HzIm);

   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR103: Failed to load Hfield.\n");
      return true;
   }
   return false;
}

void Fields::updateGrids()
{
   if (grid_Et_re) grid_Et_re->Update();
   if (grid_Et_im) grid_Et_im->Update();
   if (grid_Ez_re) grid_Ez_re->Update();
   if (grid_Ez_im) grid_Ez_im->Update();

   if (grid_Ht_re) grid_Ht_re->Update();
   if (grid_Ht_im) grid_Ht_im->Update();
   if (grid_Hz_re) grid_Hz_re->Update();
   if (grid_Hz_im) grid_Hz_im->Update();

   if (grid_Pz_re) grid_Pz_re->Update();
   if (grid_Pz_im) grid_Pz_im->Update();
}

// save the results for ParaView
void Fields::saveParaView(fem2D *fem, long unsigned int mode)
{
   if (!fem->projData->project_save_fields) return;

   stringstream ss;
   ss << fem->projData->project_name << "_frequency_" << fem->frequency << "_mode_" << mode+1;

   stringstream ssParaView;
   ssParaView << "ParaView_" << fem->projData->project_name;

   ParaViewDataCollection *pd=nullptr;
   pd=new ParaViewDataCollection(ss.str(),fem->pmesh);
   pd->SetOwnData(false);
   pd->SetPrefixPath(ssParaView.str());
   pd->RegisterField("Et_re",grid_Et_re);
   pd->RegisterField("Ez_re",grid_Ez_re);
   pd->RegisterField("Et_im",grid_Et_im);
   pd->RegisterField("Ez_im",grid_Ez_im);
   pd->RegisterField("Ht_re",grid_Ht_re);
   pd->RegisterField("Hz_re",grid_Hz_re);
   pd->RegisterField("Ht_im",grid_Ht_im);
   pd->RegisterField("Hz_im",grid_Hz_im);
   pd->RegisterField("Pz_re",grid_Pz_re);
   pd->RegisterField("Pz_im",grid_Pz_im);
   pd->SetLevelsOfDetail(fem->order);
   pd->SetDataFormat(VTKFormat::ASCII);
   pd->SetHighOrderOutput(true);
   pd->Save();
   delete pd;
}

// write as initial guess for the next iteration
void Fields::writeInitialGuess (fem2D *fem, long unsigned int mode)
{
   stringstream ss;
   ss << "_" << mode+1;

   string slash="/";

   string name_Et_re="initial_guess_Et_re";
   grid_Et_re->Save((fem->temporaryDirectory+slash+name_Et_re+ss.str()).c_str());

   string name_Et_im="initial_guess_Et_im";
   grid_Et_im->Save((fem->temporaryDirectory+slash+name_Et_im+ss.str()).c_str());

   string name_Ez_re="initial_guess_Ez_re";
   grid_Ez_re->Save((fem->temporaryDirectory+slash+name_Ez_re+ss.str()).c_str());

   string name_Ez_im="initial_guess_Ez_im";
   grid_Ez_im->Save((fem->temporaryDirectory+slash+name_Ez_im+ss.str()).c_str());
}

// calculate Pz from the E and H fields using the Poynting vector, so P=1/2 surface_integral Et x Ht*
complex<double> Fields::calculatePz(fem2D *fem)
{
   //**************************************************************************************************
   // real part
   //**************************************************************************************************

   // Et_re x Ht_re
   ParDiscreteLinearOperator Et_rexHt_re(fem->fespace_ND,fem->fespace_L2);
   VectorGridFunctionCoefficient Et_re_coef(grid_Et_re);
   Et_rexHt_re.AddDomainInterpolator(new ScalarCrossProductInterpolator(Et_re_coef));
   Et_rexHt_re.Assemble();
   Et_rexHt_re.Finalize();

   ParGridFunction grid_Pz1_re=ParGridFunction(fem->fespace_L2);
   Et_rexHt_re.Mult(*grid_Ht_re,grid_Pz1_re);

   GridFunctionCoefficient Pz1_re_coef(&grid_Pz1_re);

   // Et_im x Ht_im
   ParDiscreteLinearOperator Et_imxHt_im(fem->fespace_ND,fem->fespace_L2);
   VectorGridFunctionCoefficient Et_im_coef(grid_Et_im);
   Et_imxHt_im.AddDomainInterpolator(new ScalarCrossProductInterpolator(Et_im_coef));
   Et_imxHt_im.Assemble();
   Et_imxHt_im.Finalize();

   ParGridFunction grid_Pz2_re=ParGridFunction(fem->fespace_L2);
   Et_imxHt_im.Mult(*grid_Ht_im,grid_Pz2_re);

   GridFunctionCoefficient Pz2_re_coef(&grid_Pz2_re);

   // sum the two components

   SumCoefficient Pz_re_coef(Pz1_re_coef,Pz2_re_coef,0.5,0.5);

   // put it on the grid function

   if (grid_Pz_re) delete grid_Pz_re;
   grid_Pz_re=new ParGridFunction(fem->fespace_L2);
   grid_Pz_re->ProjectCoefficient(Pz_re_coef);

   //**************************************************************************************************
   // imaginary part
   //**************************************************************************************************

   // Et_im x Ht_re
   ParDiscreteLinearOperator Et_imxHt_re(fem->fespace_ND,fem->fespace_L2);
   Et_imxHt_re.AddDomainInterpolator(new ScalarCrossProductInterpolator(Et_im_coef));
   Et_imxHt_re.Assemble();
   Et_imxHt_re.Finalize();

   ParGridFunction grid_Pz1_im=ParGridFunction(fem->fespace_L2);
   Et_imxHt_re.Mult(*grid_Ht_re,grid_Pz1_im);

   GridFunctionCoefficient Pz1_im_coef(&grid_Pz1_im);

   // Et_re x Ht_im
   ParDiscreteLinearOperator Et_rexHt_im(fem->fespace_ND,fem->fespace_L2);
   Et_rexHt_im.AddDomainInterpolator(new ScalarCrossProductInterpolator(Et_re_coef));
   Et_rexHt_im.Assemble();
   Et_rexHt_im.Finalize();

   ParGridFunction grid_Pz2_im=ParGridFunction(fem->fespace_L2);
   Et_rexHt_im.Mult(*grid_Ht_im,grid_Pz2_im);

   GridFunctionCoefficient Pz2_im_coef(&grid_Pz2_im);

   // sum the two components

   SumCoefficient Pz_im_coef(Pz1_im_coef,Pz2_im_coef,0.5,-0.5);

   // put it on the grid function
   if (grid_Pz_im) delete grid_Pz_im;
   grid_Pz_im=new ParGridFunction(fem->fespace_L2);
   grid_Pz_im->ProjectCoefficient(Pz_im_coef);

   //**************************************************************************************************
   // integegte over the cross section
   //**************************************************************************************************

   ConstantCoefficient one(1.0);
   ParLinearForm Pz_lf(fem->fespace_L2);
   Pz_lf.AddDomainIntegrator(new DomainLFIntegrator(one));
   Pz_lf.Assemble();

   return complex<double>(Pz_lf(*grid_Pz_re),Pz_lf(*grid_Pz_im));
}

// get field values given a list of (x,y) coordinates
// set fieldValues == nullptr if calculating for Ez or Hz
bool getFieldValues (int numPoints, ParMesh *pmesh, DenseMatrix *points, ParGridFunction *grid_t_re, ParGridFunction *grid_t_im,
                     vector<complex<double>> *fieldValuesX, vector<complex<double>> *fieldValuesY)
{
   if ((int)fieldValuesX->size() != numPoints || (fieldValuesY && (int)fieldValuesY->size() != numPoints)) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: fieldValues not properly allocated.\n");
      return true;
   }

   PetscMPIInt rank,size;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // set up a datatype for data transfer
   // use an int for an array location then two doubles for a complex value
   // copied from Hsolve.c - ToDo - Make this reusable to avoid copying code.

   MPI_Datatype mpi_complex_int_type;
   int lengths[3]={1,1,1};

   MPI_Aint displacements[3];
   struct mpi_complex_int dummy;
   MPI_Aint base_address;
   MPI_Get_address(&dummy,&base_address);
   MPI_Get_address(&dummy.real,&displacements[0]);
   MPI_Get_address(&dummy.imag,&displacements[1]);
   MPI_Get_address(&dummy.location,&displacements[2]);
   displacements[0]=MPI_Aint_diff(displacements[0],base_address);
   displacements[1]=MPI_Aint_diff(displacements[1],base_address);
   displacements[2]=MPI_Aint_diff(displacements[2],base_address);

   MPI_Datatype types[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_INT};
   MPI_Type_create_struct(3,lengths,displacements,types,&mpi_complex_int_type);
   MPI_Type_commit(&mpi_complex_int_type);

   // get element indices for each point
   Array<int> elementIndices(numPoints);
   Array<IntegrationPoint> integrationPoints(numPoints);
   findPoints(pmesh, *points, elementIndices, integrationPoints);

   // space to hold transfer lengths and displacements
   int *counts_recv=(int *)malloc(size*sizeof(int));
   int *displacements_recv=(int *)malloc(size*sizeof(int));

   // space for data to gather
   struct mpi_complex_int *sendDataX=nullptr;
   struct mpi_complex_int *sendDataY=nullptr;
   struct mpi_complex_int *recvDataX=nullptr;
   struct mpi_complex_int *recvDataY=nullptr;

   sendDataX=(struct mpi_complex_int *)malloc(numPoints*sizeof(struct mpi_complex_int));
   if (fieldValuesY) sendDataY=(struct mpi_complex_int *)malloc(numPoints*sizeof(struct mpi_complex_int));
   recvDataX=(struct mpi_complex_int *)malloc(size*numPoints*sizeof(struct mpi_complex_int));                    // allocate more because a point may land in two triangles (edge)
   if (fieldValuesY) recvDataY=(struct mpi_complex_int *)malloc(size*numPoints*sizeof(struct mpi_complex_int));  // allocate more because a point may land in two triangles (edge)

   // field data
   // Will only be valid for this rank's section of the mesh
   // Store it to support gathering the data across ranks
   int count_from=0;
   int i=0;
   while (i < numPoints) {
      if (elementIndices[i] >= 0) {
         Vector tRe,tIm;
         grid_t_re->GetVectorValue(elementIndices[i],integrationPoints[i],tRe);
         grid_t_im->GetVectorValue(elementIndices[i],integrationPoints[i],tIm);

         sendDataX[count_from].location=i;
         sendDataX[count_from].real=tRe.Elem(0);
         sendDataX[count_from].imag=tIm.Elem(0);

         if (fieldValuesY) {
            sendDataY[count_from].location=i;
            sendDataY[count_from].real=tRe.Elem(1);
            sendDataY[count_from].imag=tIm.Elem(1);
         }

         count_from++;
      }
      i++;
   }

   // gather the transfer lengths 
   if (MPI_Gather(&count_from,1,MPI_INT,counts_recv,1,MPI_INT,0,PETSC_COMM_WORLD)) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR111: Failed to gather transfer lengths.\n");
      return true;
   }

   // calculate displacements
   if (rank == 0) {
      displacements_recv[0]=0;
      i=1;
      while (i < size) {
         displacements_recv[i]=displacements_recv[i-1]+counts_recv[i-1];
         i++;
      }
   }

   // get the total number of points
   // Could be greater than numPoints due to duplicated points on a triangle boundary.

   int totalPoints=0;
   if (rank == 0) {
      i=0;
      while (i < size) {
         totalPoints+=counts_recv[i];
         i++;
      }
   }

   if (MPI_Bcast(&totalPoints,1,MPI_INT,0,PETSC_COMM_WORLD)) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR116: Failed to broadcast totalPoints.\n");
      return true;
   }

   // send all to rank 0
   if (MPI_Gatherv(sendDataX,count_from,mpi_complex_int_type,recvDataX,counts_recv,displacements_recv,mpi_complex_int_type,0,PETSC_COMM_WORLD)) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR112: Failed to gather sendDataX.\n");
      return true;
   }

   if (fieldValuesY && MPI_Gatherv(sendDataY,count_from,mpi_complex_int_type,recvDataY,counts_recv,displacements_recv,mpi_complex_int_type,0,PETSC_COMM_WORLD)) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR113: Failed to gather sendDataY.\n");
      return true;
   }

   // cleanup
   if (counts_recv) {free(counts_recv); counts_recv=nullptr;}
   if (displacements_recv) {free(displacements_recv); displacements_recv=nullptr;}
   if (sendDataX) {free(sendDataX); sendDataX=nullptr;}
   if (fieldValuesY && sendDataY) {free(sendDataY); sendDataY=nullptr;}

   // broadcast

   if (MPI_Bcast(recvDataX,totalPoints,mpi_complex_int_type,0,PETSC_COMM_WORLD)) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR114: Failed to broadcast recvDataX.\n");
      return true;
   }

   if (fieldValuesY && MPI_Bcast(recvDataY,totalPoints,mpi_complex_int_type,0,PETSC_COMM_WORLD)) {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR115: Failed to broadcast recvDataY.\n");
      return true;
   }

   // initialize data to enable detection of errors
   i=0;
   while (i < numPoints) {
      (*fieldValuesX)[i]=complex<double>(DBL_MAX,DBL_MAX);
      if (fieldValuesY) (*fieldValuesY)[i]=complex<double>(DBL_MAX,DBL_MAX);
      i++;
   }

   // transfer the results 
   i=0;
   while (i < totalPoints) {
      if (recvDataX[i].location >= 0) {
         (*fieldValuesX)[recvDataX[i].location]=complex<double>(recvDataX[i].real,recvDataX[i].imag);
         if (fieldValuesY) (*fieldValuesY)[recvDataY[i].location]=complex<double>(recvDataY[i].real,recvDataY[i].imag);
      }
      i++;
   }

   // cleanup
   if (recvDataX) {free(recvDataX); recvDataX=nullptr;}
   if (fieldValuesY && recvDataY) {free(recvDataY); recvDataY=nullptr;}

   // check for missing points
   i=0;
   while (i < numPoints) {
      if ((*fieldValuesX)[i] == complex<double>(DBL_MAX,DBL_MAX) || (fieldValuesY && (*fieldValuesY)[i] == complex<double>(DBL_MAX,DBL_MAX))) {
         PetscPrintf(PETSC_COMM_WORLD,"WARNING: Integration point (%g,%g) is not detected inside the mesh.\n",
                                      points->Elem(0,i),points->Elem(1,i));
      }
      i++;
   }

   MPI_Type_free(&mpi_complex_int_type);

   return false;
}

// field = false => electric field
// field = true  => magnetic field
complex<double> Fields::calculateLineIntegral (fem2D *fem, Boundary *boundary, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase, bool field)
{
   long unsigned int max=-1;
   double x1,y1,x2,y2;
   double length, nhatx, nhaty;
   ParGridFunction *grid_t_re, *grid_t_im;
   vector<complex<double>> fieldValuesX,fieldValuesY;
   double xi1,yi1,xi2,yi2,segmentLength;
   complex<double> valueX,valueY;
   bool has_printed_integrator=false;
   bool has_printed_numerical=false;

   fieldValuesX.resize(fem->pointsCount);
   fieldValuesY.resize(fem->pointsCount);

   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // set the field type
   if (field) {
      grid_t_re=grid_Ht_re;
      grid_t_im=grid_Ht_im;
   } else {
      grid_t_re=grid_Et_re;
      grid_t_im=grid_Et_im;
   }

   // sum the results over the entire path
   complex<double> integral=complex<double>(0,0);

   // see if this boundary is on a mesh border
   struct EdgeAttribute test_attribute;
   test_attribute.boundary=boundary->get_attribute();
   test_attribute.path=-1;
   test_attribute.segment=-1;
   long unsigned int boundary_index=borderDatabase->exists_test_boundary_only(test_attribute);
   bool on_border=true;
   if (boundary_index == max) on_border=false;

   // loop through the paths
   long unsigned int iPath=0;
   while (iPath < boundary->get_path_size()) {

      bool reverse=boundary->get_reverse(iPath);

      vector<Path *> *pathList=boundaryDatabase->get_pathList();
      Path *path=(*pathList)[boundary->get_path(iPath)];

      // loop through the path segments

      long unsigned int limit=path->get_points_size();
      if (path->is_closed()) limit++;

      x2=path->get_point_x(0);
      y2=path->get_point_y(0);

      long unsigned int j=1;
      while (j < limit) {

         x1=x2;
         y1=y2;

         if (path->is_closed()) {
            if (j < limit-1) {
               x2=path->get_point_x(j);
               y2=path->get_point_y(j);
            } else {
               x2=path->get_point_x(0);
               y2=path->get_point_y(0);
            }
         } else {
            x2=path->get_point_x(j);
            y2=path->get_point_y(j);
         }

         // get a unit vector in the direction of the supplied line
         length=sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
         nhatx=(x2-x1)/length;
         nhaty=(y2-y1)/length;

         if (on_border) {

            // integrate with an MFEM integrator

            has_printed_integrator=true;
            if (! has_printed_integrator) {
               PetscPrintf(PETSC_COMM_WORLD,"            calculating integral with MFEM integrator\n");
               has_printed_integrator=true;
            }

            // restrict to the needed border elements

            // knock out all
            fem->border_attributes=0;

            // add back in this border, path, and segment
            long unsigned int k=0;
            while (k < borderDatabase->get_size()) {
               Border *border=borderDatabase->get_border(k);
               if ((int)border->get_mode(boundary->get_mode()).boundary == boundary->get_attribute() && 
                        border->get_mode(boundary->get_mode()).path == iPath &&
                        border->get_mode(boundary->get_mode()).segment == j-1) {
                  fem->border_attributes[k-1]=1;
               }
               k++;
            }

            // unit vector for this segment
            Vector *unit=new Vector(2);
            unit->Elem(0)=nhatx;
            unit->Elem(1)=nhaty;
            VectorCoefficient *unitVec=new VectorConstantCoefficient(*unit);

            ParLinearForm lf(fem->fespace_ND);
            lf.AddBoundaryIntegrator(new VectorFEDomainLFIntegrator(*unitVec),fem->border_attributes);
            lf.Assemble();

            double I_re=lf(*grid_t_re);
            double I_im=lf(*grid_t_im);

            if (reverse) integral-=complex<double>(I_re,I_im);
            else         integral+=complex<double>(I_re,I_im);

            delete unit;
            delete unitVec;

         } else {

            // numerically integrate

            has_printed_numerical=true;
            if (! has_printed_numerical) {
               PetscPrintf(PETSC_COMM_WORLD,"            calculating integral with numerical integration\n");
               has_printed_numerical=true;
            }

            // get points along the integration line
            DenseMatrix points(2,fem->pointsCount);
            double theta=atan2(y2-y1,x2-x1);
            int i=0;
            while (i < fem->pointsCount) {
               points.Elem(0,i)=x1+cos(theta)*length*(double)i/((double)(fem->pointsCount-1));
               points.Elem(1,i)=y1+sin(theta)*length*(double)i/((double)(fem->pointsCount-1));
               i++;
            }

            // get the field values
            if (getFieldValues (fem->pointsCount,fem->pmesh,&points,grid_t_re,grid_t_im,&fieldValuesX,&fieldValuesY)) {
               return complex<double>(DBL_MAX,DBL_MAX);  // fatal error
            }

            // calculate the integral
            // ToDo - break this up into sections by ranks to parallel process; not really worth it for small pointsCount

            i=0;
            while (i < fem->pointsCount-1) {

               xi1=points.Elem(0,i);
               yi1=points.Elem(1,i);

               xi2=points.Elem(0,i+1);
               yi2=points.Elem(1,i+1);

               segmentLength=sqrt((xi2-xi1)*(xi2-xi1)+(yi2-yi1)*(yi2-yi1));

               // enable the calculation to continue in case a point lies just outside the mesh
               if (fieldValuesX[i] == complex<double>(DBL_MAX,DBL_MAX) || fieldValuesX[i+1] == complex<double>(DBL_MAX,DBL_MAX) ||
                   fieldValuesY[i] == complex<double>(DBL_MAX,DBL_MAX) || fieldValuesY[i+1] == complex<double>(DBL_MAX,DBL_MAX)) {i++; continue;}

               valueX=(fieldValuesX[i]+fieldValuesX[i+1])/2;
               valueY=(fieldValuesY[i]+fieldValuesY[i+1])/2;

               if (reverse) integral-=(valueX*nhatx+valueY*nhaty)*segmentLength;
               else         integral+=(valueX*nhatx+valueY*nhaty)*segmentLength;

               i++;
            }
         }

         j++;
      }
      iPath++;
   }

   return integral;
}

void Fields::calculateFieldPoints (fem2D *fem, int mode, FieldPointDatabase *fieldPointDatabase)
{
   FieldPoint fieldPoint;
   vector<complex<double>> Etx,Ety,Etz,Htx,Hty,Htz;

   // skip if there are no frequency points to calculate
   if (fem->projData->field_points_count <= 0) return;

   // intermediate storage to call getFieldValues
   Etx.resize(fem->projData->field_points_count);
   Ety.resize(fem->projData->field_points_count);
   Etz.resize(fem->projData->field_points_count);

   Htx.resize(fem->projData->field_points_count);
   Hty.resize(fem->projData->field_points_count);
   Htz.resize(fem->projData->field_points_count);


   // transfer the (x,y) locations to a new data structure
   DenseMatrix points(2,fem->projData->field_points_count);
   int i=0;
   while (i < fem->projData->field_points_count) {
      points.Elem(0,i)=fem->projData->field_points_x[i];
      points.Elem(1,i)=fem->projData->field_points_y[i];
      i++;
   }

   // This setup is inefficient because each call to getFieldValues re-calculates all the information about the mesh.
   // It is assumed that the field point calculation is only used for validation work, so the loss of efficiency is acceptable.

   if (getFieldValues(fem->projData->field_points_count,fem->pmesh,&points,grid_Et_re,grid_Et_im,&Etx,&Ety)) return;
   if (getFieldValues(fem->projData->field_points_count,fem->pmesh,&points,grid_Ez_re,grid_Ez_im,&Etz,nullptr)) return;
   if (getFieldValues(fem->projData->field_points_count,fem->pmesh,&points,grid_Ht_re,grid_Ht_im,&Htx,&Hty)) return;
   if (getFieldValues(fem->projData->field_points_count,fem->pmesh,&points,grid_Hz_re,grid_Hz_im,&Htz,nullptr)) return;

   i=0;
   while (i < fem->projData->field_points_count) {

      fieldPoint.set_frequency(fem->frequency);
      fieldPoint.set_mode(mode);
      fieldPoint.set_x(fem->projData->field_points_x[i]);
      fieldPoint.set_y(fem->projData->field_points_y[i]);

      fieldPoint.set_Ex_re(real(Etx[i]));
      fieldPoint.set_Ex_im(imag(Etx[i]));
      fieldPoint.set_Ey_re(real(Ety[i]));
      fieldPoint.set_Ey_im(imag(Ety[i]));
      fieldPoint.set_Ez_re(real(Etz[i]));
      fieldPoint.set_Ez_im(imag(Etz[i]));

      fieldPoint.set_Hx_re(real(Htx[i]));
      fieldPoint.set_Hx_im(imag(Htx[i]));
      fieldPoint.set_Hy_re(real(Hty[i]));
      fieldPoint.set_Hy_im(imag(Hty[i]));
      fieldPoint.set_Hz_re(real(Htz[i]));
      fieldPoint.set_Hz_im(imag(Htz[i]));

      fieldPointDatabase->push(&fieldPoint);

      i++;
   }
}

// perturbational loss calculated over each boundary segment
double Fields::calculatePerturbationalLoss(fem2D *fem, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase, MaterialDatabase *materialDatabase)
{
   double tolerance=1e-12;
   string indent="   ";
   double perturbationalLoss=0;

   // loop over the border attributes
   int i=0;
   while (i <= (int)fem->border_attributes.Size()) {

      // find the boundary associated with this border
      Border *border=borderDatabase->get_border(i);
      if (border->has_boundary()) {
         long unsigned int boundary_index=border->get_boundary().boundary;
         Boundary *boundary=boundaryDatabase->get_boundary(boundary_index);
         if (boundary->is_surface_impedance()) {

            // get the surface impedance
            double Rs=0;
            Material *material=materialDatabase->get(boundary->get_material());
            if (material) {
               Rs=material->get_Rs(fem->projData->solution_temperature, fem->frequency, tolerance, indent);
            } else {
               PetscPrintf(PETSC_COMM_WORLD,"ERROR105: Material \"%s\" does not exist in the materials database.  Defaulting to PEC.\n",
                                            boundary->get_material().c_str());
            }

            // knock out all boundaries and modes
            fem->border_attributes=0;

            // put this one back
            fem->border_attributes[i-1]=1;

            // integrate the Ht contribution

            // real
            VectorGridFunctionCoefficient Ht_re_coef(grid_Ht_re);
            ParLinearForm Ht_re_lf(fem->fespace_ND);
            Ht_re_lf.AddBoundaryIntegrator(new VectorFEDomainLFIntegrator(Ht_re_coef),fem->border_attributes);
            Ht_re_lf.Assemble();
            double Ht_re_mag2=Ht_re_lf(*grid_Ht_re);

            // imag
            VectorGridFunctionCoefficient Ht_im_coef(grid_Ht_im);
            ParLinearForm Ht_im_lf(fem->fespace_ND);
            Ht_im_lf.AddBoundaryIntegrator(new VectorFEDomainLFIntegrator(Ht_im_coef),fem->border_attributes);
            Ht_im_lf.Assemble();
            double Ht_im_mag2=Ht_im_lf(*grid_Ht_im);

            // integrate the Hz contribution

            // real
            GridFunctionCoefficient Hz_re_coef(grid_Hz_re);
            ParLinearForm Hz_re_lf(fem->fespace_H1);
            Hz_re_lf.AddBoundaryIntegrator(new DomainLFIntegrator(Hz_re_coef),fem->border_attributes);
            Hz_re_lf.Assemble();
            double Hz_re_mag2=Hz_re_lf(*grid_Hz_re);

            // imag
            GridFunctionCoefficient Hz_im_coef(grid_Hz_im);
            ParLinearForm Hz_im_lf(fem->fespace_H1);
            Hz_im_lf.AddBoundaryIntegrator(new DomainLFIntegrator(Hz_im_coef),fem->border_attributes);
            Hz_im_lf.Assemble();
            double Hz_im_mag2=Hz_im_lf(*grid_Hz_im);

            perturbationalLoss+=0.5*Rs*(Ht_re_mag2+Ht_im_mag2+Hz_re_mag2+Hz_im_mag2);

         }
      }
      i++;
   }

   return perturbationalLoss;
}

void Fields::print()
{
   PetscPrintf(PETSC_COMM_WORLD,"      Fields:\n");
   PetscPrintf(PETSC_COMM_WORLD,"         this=%p\n",this);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Et_re=%p\n",grid_Et_re);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Et_im=%p\n",grid_Et_im);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Ez_re=%p\n",grid_Ez_re);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Ez_im=%p\n",grid_Ez_im);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Ht_re=%p\n",grid_Ht_re);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Ht_im=%p\n",grid_Ht_im);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Hz_re=%p\n",grid_Hz_re);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Hz_im=%p\n",grid_Hz_im);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Pz_re=%p\n",grid_Pz_re);
   PetscPrintf(PETSC_COMM_WORLD,"         grid_Pz_im=%p\n",grid_Pz_im);
   PetscPrintf(PETSC_COMM_WORLD,"         EfieldRe=%p\n",EfieldRe);
   PetscPrintf(PETSC_COMM_WORLD,"         EfieldIm=%p\n",EfieldIm);
   PetscPrintf(PETSC_COMM_WORLD,"         HfieldRe=%p\n",HfieldRe);
   PetscPrintf(PETSC_COMM_WORLD,"         HfieldIm=%p\n",HfieldIm);
   PetscPrintf(PETSC_COMM_WORLD,"         xScale=%g\n",xScale);
   PetscPrintf(PETSC_COMM_WORLD,"         yScale=%g\n",yScale);
}

void Fields::save_grid_E(long unsigned int i)
{
   char filename[128];

   sprintf (filename,"%s_mode_%ld","grid_Et_re",i);
   grid_Et_re->Save(filename);

   sprintf (filename,"%s_mode_%ld","grid_Et_im",i);
   grid_Et_im->Save(filename);

   sprintf (filename,"%s_mode_%ld","grid_Ez_re",i);
   grid_Ez_re->Save(filename);

   sprintf (filename,"%s_mode_%ld","grid_Ez_im",i);
   grid_Ez_im->Save(filename);

}

void Fields::saveAsOne_grid_E(long unsigned int i)
{
   char filename[128];

   sprintf (filename,"%s_mode_%ld_asOne","grid_Et_re",i);
   grid_Et_re->SaveAsOne(filename);

   sprintf (filename,"%s_mode_%ld_asOne","grid_Et_im",i);
   grid_Et_im->SaveAsOne(filename);

   sprintf (filename,"%s_mode_%ld_asOne","grid_Ez_re",i);
   grid_Ez_re->SaveAsOne(filename);

   sprintf (filename,"%s_mode_%ld_asOne","grid_Ez_im",i);
   grid_Ez_im->SaveAsOne(filename);

}


//--------------------------------------------------------------------------------------------------------------------
// Mode
//--------------------------------------------------------------------------------------------------------------------

Mode::Mode()
{
   modeNumber=0;
   Pzavg=complex<double>(0,0);
   alpha=0;
   beta=0;
   perturbationalLoss=0;
   fields=new Fields();
}

bool Mode::buildFields(fem2D *fem)
{
   bool fail=false;
   fail=fields->build(fem,modeNumber);
   if (! fail) calculatePz(fem);
   return fail;
}

void Mode::updateGrids(fem2D *fem)
{
   fields->updateGrids();
}

void Mode::saveParaView(fem2D *fem)
{
   fields->saveParaView(fem,modeNumber);
}

void PrintError (struct mpi_double_int_int *a, int i)
{
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   cout << "[" << rank << "]  i=" << i << "  position=" << a->location << "  value=" << a->value << "  rank=" << a->rank << endl;
}

void PrintErrors (struct mpi_double_int_int *A, int lenA)
{
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   cout << "[" << rank << "] PrintErrors:" << endl;
   cout << "[" << rank << "]   A=" << A << endl;
   cout << "[" << rank << "]  lenA=" << lenA << endl;
   int i=0;
   while (i < lenA) {
      PrintError(&A[i],i);
      i++;
   }
}

// merge the A array into B
// B is allocated and must be freed elsewhere
int MergeErrors (struct mpi_double_int_int *A, int lenA, struct mpi_double_int_int **B, int *lenB)
{
   int size,rank;

   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // set up a datatype for data transfer
   // use an int for an array location, an int for the rank, and a double for a value

   MPI_Datatype mpi_double_int_int_type;
   int lengths[3]={1,1,1};

   MPI_Aint displacements[3];
   struct mpi_double_int_int dummy;
   MPI_Aint base_address;
   MPI_Get_address(&dummy,&base_address);
   MPI_Get_address(&dummy.value,&displacements[0]);
   MPI_Get_address(&dummy.location,&displacements[1]);
   MPI_Get_address(&dummy.rank,&displacements[2]);
   displacements[0]=MPI_Aint_diff(displacements[0],base_address);
   displacements[1]=MPI_Aint_diff(displacements[1],base_address);
   displacements[2]=MPI_Aint_diff(displacements[2],base_address);

   MPI_Datatype types[3]={MPI_DOUBLE,MPI_INT,MPI_INT};
   MPI_Type_create_struct(3,lengths,displacements,types,&mpi_double_int_int_type);
   MPI_Type_commit(&mpi_double_int_int_type);

   // space to hold transfer lengths and displacements
   int *counts_recv,*displacements_recv;
   counts_recv=(int *)malloc(size*sizeof(int));
   displacements_recv=(int *)malloc(size*sizeof(int));

   // gather the transfer lengths
   if (MPI_Gather(&lenA,1,MPI_INT,counts_recv,1,MPI_INT,0,PETSC_COMM_WORLD)) {
      PetscPrintf (PETSC_COMM_WORLD,"ERROR119: Failed to gather data.\n");
      return 1;
   }

   // calculate the total length of the data
   *lenB=0;
   if (rank == 0) {
      int i=0;
      while (i < size) {
         (*lenB)+=counts_recv[i];
         i++;
      }
   }

   if (MPI_Bcast(lenB,1,MPI_INT,0,PETSC_COMM_WORLD)) {
      PetscPrintf (PETSC_COMM_WORLD,"ERROR120: Failed to broadcast total.\n");
      return 1;
   }

   // calculate displacements
   if (rank == 0) {
      displacements_recv[0]=0;
      int i=1;
      while (i < size) {
         displacements_recv[i]=displacements_recv[i-1]+counts_recv[i-1];
         i++;
      }
   }

   // space for the received data
   *B=(struct mpi_double_int_int *)malloc((*lenB)*sizeof(struct mpi_double_int_int));

   // send all to rank 0
   if (MPI_Gatherv(A,lenA,mpi_double_int_int_type,*B,counts_recv,displacements_recv,mpi_double_int_int_type,0,PETSC_COMM_WORLD)) {
      PetscPrintf (PETSC_COMM_WORLD,"ERROR121: Failed to gather data.\n");
      return 1;
   }

   if (counts_recv) {free(counts_recv); counts_recv=nullptr;}
   if (displacements_recv) {free(displacements_recv); displacements_recv=nullptr;}

   if (MPI_Bcast(*B,*lenB,mpi_double_int_int_type,0,PETSC_COMM_WORLD)) {
      PetscPrintf (PETSC_COMM_WORLD,"ERROR122: Failed to broadcast data.\n");
      return 1;
   }

   MPI_Type_free(&mpi_double_int_int_type);

   return 0;
}

// refine on a partial component of the tangengial electric field
// ToDo: write a custom integrator to enable adaptive meshing on the total electric field
bool Mode::ZZrefineMesh (fem2D *fem, ConvergenceDatabase *convergenceDatabase)
{
   long unsigned int i;
   int j;
   bool use_initial_guess=true;
   int rank;

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // skip refinement if the mode is converged to avoid unnecessary mesh creation
   if (! fem->projData->refinement_refine_converged_modes && convergenceDatabase->is_converged(modeNumber)) {
      if (fem->projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"         mode %ld skipped\n",modeNumber+1);
      return use_initial_guess;
   }

   // continuing

   if (fem->projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"         mode %ld:\n",modeNumber+1);

   // flux
   L2_FECollection *fec_L2_flux=new L2_FECollection(fem->order-1,fem->pmesh->Dimension());
   ParFiniteElementSpace *fespace_L2_flux=new ParFiniteElementSpace(fem->pmesh,fec_L2_flux);

   // smoothed flux
   H1_FECollection *fec_H1_flux=new H1_FECollection(fem->order,fem->pmesh->Dimension());
   ParFiniteElementSpace *fespace_H1_flux=new ParFiniteElementSpace(fem->pmesh,fec_H1_flux);

   // estimator
   L2ZienkiewiczZhuEstimator *estimator;
   CurlCurlIntegrator CCinteg=CurlCurlIntegrator ();
   estimator=new L2ZienkiewiczZhuEstimator(CCinteg,*fields->get_grid_Et_re(),fespace_L2_flux,fespace_H1_flux);

   // local errors
   Vector localErrors;
   localErrors=estimator->GetLocalErrors();
   //cout << "estimator->GetTotalError()=" << estimator->GetTotalError() << endl;

   // element centers and local error limits
   DenseMatrix centers(2,localErrors.Size());
   Vector center(2);
   j=0;
   while (j < localErrors.Size()) {
      fem->pmesh->GetElementCenter (j,center);
      centers.Elem(0,j)=center.Elem(0);
      centers.Elem(1,j)=center.Elem(1);
      j++;
   }

   // get element indices for each point
   Array<int> element_indices(localErrors.Size());
   Array<IntegrationPoint> integrationPoints(localErrors.Size());
   findPoints(fem->pmesh, centers, element_indices, integrationPoints);

   // store the local errors in a struct for merging
   struct mpi_double_int_int *localErr=(struct mpi_double_int_int *)malloc(localErrors.Size()*sizeof(struct mpi_double_int_int));
   j=0;
   while (j < localErrors.Size()) {
      localErr[j].rank=rank;
      localErr[j].location=element_indices[j];
      localErr[j].value=localErrors[j];
      j++;
   }

   // merge the local errors for global mesh refinement
   struct mpi_double_int_int *globalErrors=nullptr;
   int globalLen;
   if (MergeErrors (localErr,localErrors.Size(),&globalErrors,&globalLen)) return false;

   if (localErr) {free(localErr); localErr=nullptr;}

   // bubble sort - ToDo - replace with something more efficient and preferrably parallel
   if (fem->projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"            sorting %d elements ...\n",globalLen);
   bool found=true;
   while (found) {
      found=false;
      j=0;
      while (j < globalLen-1) {
         if (globalErrors[j].value < globalErrors[j+1].value) {
            struct mpi_double_int_int temp;
            temp=globalErrors[j];
            globalErrors[j]=globalErrors[j+1];
            globalErrors[j+1]=temp;

            found=true;
         }
         j++;
      }
   }

   // get the global refinement count

   // set the count by keeping the elements that have an error that is as a fraction of the max error
   long unsigned int globalRefineCount=0;
   while (globalRefineCount < (long unsigned int)globalLen) {
      if (globalErrors[globalRefineCount].value < fem->projData->mesh_refinement_cutoff*globalErrors[0].value) break;
      globalRefineCount++;
   }

   // cap the count
   if (globalRefineCount > globalLen*fem->projData->mesh_refinement_fraction/fem->meshScale)
         globalRefineCount=globalLen*fem->projData->mesh_refinement_fraction/fem->meshScale;

   if (globalRefineCount == 0) globalRefineCount=1;

   // count the number of local refinements
   long unsigned int localRefineCount=0;
   i=0;
   while (i < globalRefineCount) {
      if (globalErrors[i].rank == rank) localRefineCount++;
      i++;
   }

   // refine locally
   Array<int> localRefineList(localRefineCount);
   long unsigned int index=0;
   i=0;
   while (i < globalRefineCount) {
      if (globalErrors[i].rank == rank) {
         localRefineList[index]=globalErrors[i].location;
         index++;
      }
      i++;
   }

   if (localRefineList.Size() == 1) {
      if (fem->projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"            refining %d element ...\n",globalRefineCount);
   } else {
      if (fem->projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"            refining %d elements ...\n",globalRefineCount);
   }

   fem->pmesh->GeneralRefinement(localRefineList);

   // cleanup
   if (globalErrors) {free(globalErrors); globalErrors=nullptr;}

   if (fem->projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"            new mesh size: %d\n",GetGlobalNE(fem->pmesh));

   delete fec_H1_flux;
   delete fec_L2_flux;
   delete estimator;

   // see if the mesh has grown too much for more refinement
   // This can happen when the starting mesh is very small.
   if (GetGlobalNE(fem->pmesh) > 1.5*fem->startingMeshSize) { // 1.5
      use_initial_guess=false;  // older way
      use_initial_guess=true;
      if (fem->projData->output_show_refining_mesh) PetscPrintf(PETSC_COMM_WORLD,"            stopping refinement due to mesh excess mesh growth\n");
   }

   return use_initial_guess;
}

void Mode::calculatePz(fem2D *fem)
{
   Pzavg=fields->calculatePz(fem);
}

void Mode::calculateVoltages (fem2D *fem, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase)
{
   // loop through the modes 
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {

      // loop through the boundary definitions 
      bool found=false;
      long unsigned int j=0;
      while (j < boundaryDatabase->get_boundary_size()) {

         // pick off voltage definitions with for the given mode
         if (boundaryDatabase->get_boundary(j)->is_mode_voltage() && boundaryDatabase->get_boundary(j)->get_mode()-1 == (int)i) {

            // calculate line integrals
            voltage.push_back(-fields->calculateLineIntegral(fem,boundaryDatabase->get_boundary(j),boundaryDatabase,borderDatabase,false));
            found=true;
            break;
         }
         j++;
      }
      if (!found) voltage.push_back(complex<double>(DBL_MAX,DBL_MAX));
      i++;
   }
}

void Mode::calculateCurrents (fem2D *fem, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase)
{
   // loop through the modes
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {

      // loop through the mode definitions
      bool found=false;
      long unsigned int j=0;
      while (j < boundaryDatabase->get_boundary_size()) {

         // pick off current definitions for the given mode
         if (boundaryDatabase->get_boundary(j)->is_mode_current() && boundaryDatabase->get_boundary(j)->get_mode()-1 == (int)i) {

            // calculate line integrals
            current.push_back(fields->calculateLineIntegral(fem,boundaryDatabase->get_boundary(j),boundaryDatabase,borderDatabase,true));
            found=true;
            break;
         }
         j++;
      }
      if (!found) current.push_back(complex<double>(DBL_MAX,DBL_MAX));
      i++;
   } 
}

complex<double> Mode::get_current (long unsigned int i)
{
   if (i >= current.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Mode::get_current out of range.\n");
      return complex<double>(DBL_MAX,DBL_MAX);
   }

   return current[i];
}

complex<double> Mode::get_voltage (long unsigned int i)
{
   if (i >= voltage.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Mode::get_voltage out of range.\n");
      return complex<double>(DBL_MAX,DBL_MAX);
   }

   return voltage[i];
}

void Mode::calculatePerturbationalLoss(fem2D *fem, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase, MaterialDatabase *materialDatabase)
{
   perturbationalLoss=fields->calculatePerturbationalLoss(fem,boundaryDatabase,borderDatabase,materialDatabase);
}

void Mode::print()
{
   PetscPrintf(PETSC_COMM_WORLD,"Mode:\n");
   PetscPrintf(PETSC_COMM_WORLD,"   this=%p\b",this);
   PetscPrintf(PETSC_COMM_WORLD,"   modeNumber=%ld\n",modeNumber);

   unsigned long int i=0;
   while (i < current.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"   current[%ld]=(%g,%g)\n",i,real(current[i]),imag(current[i]));
      i++;
   }

   i=0;
   while (i < voltage.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"   voltage[%ld]=(%g,%g)\n",i,real(voltage[i]),imag(voltage[i]));
      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"   Pzavg=(%g,%g)\n",real(Pzavg),imag(Pzavg));
   PetscPrintf(PETSC_COMM_WORLD,"   alpha=%g\n",alpha);
   PetscPrintf(PETSC_COMM_WORLD,"   beta=%g\n",beta);
   PetscPrintf(PETSC_COMM_WORLD,"   perturbationalLoss=%g\n",perturbationalLoss);
   PetscPrintf(PETSC_COMM_WORLD,"   fields=%p\n",fields);
   PetscPrintf(PETSC_COMM_WORLD,"   fields:\n");
   fields->print();
}

void Mode::save_grid_E (long unsigned int i)
{
   if (fields == nullptr) return;
   fields->save_grid_E(i);
}

void Mode::saveAsOne_grid_E (long unsigned int i)
{
   if (fields == nullptr) return;
   fields->saveAsOne_grid_E(i);
}

void Mode::calculateFieldPoints (fem2D *fem, FieldPointDatabase *fieldPointDatabase)
{
   fields->calculateFieldPoints(fem,modeNumber,fieldPointDatabase);
}

Mode::~Mode()
{
   if (fields) {delete fields; fields=nullptr;}
}

//--------------------------------------------------------------------------------------------------------------------
// ModeDatabase
//--------------------------------------------------------------------------------------------------------------------

bool ModeDatabase::buildFields(fem2D *fem, double *alpha, double *beta)
{
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {
      Mode *mode=new Mode();
      modeList.push_back(mode);

      mode->set_modeNumber(i);    // zero based
      mode->set_alpha(alpha[i]);
      mode->set_beta(beta[i]);
      if (mode->buildFields(fem)) return true;

      i++;
   }
   return false;
}

void ModeDatabase::saveParaView(fem2D *fem)
{
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {
      modeList[i]->saveParaView(fem);
      i++;
   }
}

bool ModeDatabase::ZZrefineMeshes (fem2D *fem, ConvergenceDatabase *convergenceDatabase)
{
   bool use_initial_guess=false;

   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {
      use_initial_guess=modeList[i]->ZZrefineMesh(fem,convergenceDatabase);

      // update fespaces
      fem->fespace_ND->Update();
      fem->fespace_H1->Update();
      fem->fespace_L2->Update();

      // update the grids
      long unsigned int j=0;
      while (j < (long unsigned int)fem->projData->solution_active_mode_count) {
         modeList[j]->updateGrids(fem);
         j++;
      }

      if (! use_initial_guess) break;

      i++;
   }

   return use_initial_guess;
}

void ModeDatabase::writeInitialGuess (fem2D *fem)
{
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {
      modeList[i]->writeInitialGuess (fem);
      i++;
   }
}

void ModeDatabase::calculatePz(fem2D *fem)
{
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {
      modeList[i]->calculatePz(fem);
      i++;
   }
}

void ModeDatabase::calculateImpedanceMatrix (fem2D *fem, const char *definition, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase)
{
   complex<double>two=complex<double>(2,0);
   long unsigned int n=fem->projData->solution_active_mode_count;

   if (strcmp(fem->projData->solution_impedance_definition,"none") == 0) return;

   // get the voltages and currents
   long unsigned int i=0;
   while (i < n) {
      modeList[i]->calculateCurrents(fem,boundaryDatabase,borderDatabase);
      modeList[i]->calculateVoltages(fem,boundaryDatabase,borderDatabase);
      i++;
   }

   // modal calculation
   // Assumes that the voltage and/or current paths are set up to make impedance calculations on the mode as a whole.
   // This means that the resulting impedance matrix is diagonal since the modes are orthogonal.
   if (is_modal_impedance(fem->projData->solution_impedance_calculation)) {

      // allocate memory that must be freed elsewhere
      // space for an nx1 array, where n is the number of modes
      if (Zvi) {free(Zvi); Zvi=nullptr;}
      if (Zpv) {free(Zpv); Zpv=nullptr;}
      if (Zpi) {free(Zpi); Zpi=nullptr;}
      Zvi=(complex<double> *) malloc(n*sizeof(complex<double>));
      Zpv=(complex<double> *) malloc(n*sizeof(complex<double>));
      Zpi=(complex<double> *) malloc(n*sizeof(complex<double>));

      if (fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Mode\n");
      }

      int mode=-1;
      long unsigned int i=0;
      while (i < n) {

         long unsigned int k=0;
         while (k < modeList.size()) {
            if (modeList[k]->get_modeNumber() == i) {mode=(int)k; break;}
            k++;
         }
         if (mode == -1) {
            PetscPrintf(PETSC_COMM_WORLD,"ASSERT: ModeDatabase::calculateImpedanceMatrix did not find a mode number.\n");
            return;
         }

         complex<double> voltage=modeList[mode]->get_voltage(mode);
         bool validVoltage=true;
         if (voltage == complex<double>(DBL_MAX,DBL_MAX)) validVoltage=false;

         complex<double> current=modeList[mode]->get_current(mode);
         bool validCurrent=true;
         if (current == complex<double>(DBL_MAX,DBL_MAX)) validCurrent=false;

         complex<double> Pzavg=modeList[mode]->get_Pzavg();

         if (fem->projData->debug_show_impedance_details) {
            PetscPrintf(PETSC_COMM_WORLD,"              %3ld",i+1);
            if (validVoltage) {PetscPrintf(PETSC_COMM_WORLD,"  voltage (V): (%g,%g)\n", real(voltage),imag(voltage));}
            else              {PetscPrintf(PETSC_COMM_WORLD,"  voltage (V): not defined\n");}
            if (validCurrent) {PetscPrintf(PETSC_COMM_WORLD,"                   current (I): (%g,%g)\n",real(current),imag(current));}
            else              {PetscPrintf(PETSC_COMM_WORLD,"                   current (I): not defined\n");}
                               PetscPrintf(PETSC_COMM_WORLD,"                   Pz (Pz,avg): (%g,%g)\n", real(Pzavg),imag(Pzavg));
         }

         if (validVoltage && validCurrent) Zvi[i]=voltage/current;
         else Zvi[i]=complex<double>(DBL_MAX,DBL_MAX);

         if (validVoltage) Zpv[i]=0.5*voltage*conj(voltage)/conj(Pzavg);
         else Zpv[i]=complex<double>(DBL_MAX,DBL_MAX);

         if (validCurrent) Zpi[i]=2*Pzavg/(current*conj(current));
         else Zpi[i]=complex<double>(DBL_MAX,DBL_MAX);

         i++;
      }

      if (fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Impedance (VI):\n");
         printImpedanceVector (Zvi,n);
      }

      if (fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Impedance (PV):\n");
         printImpedanceVector (Zpv,n);
      }

      if (fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Impedance (PI):\n");
         printImpedanceVector (Zpi,n);
      }

   // line calculation
   // Assumes that the voltage and/or current paths are set up for individual lines.
   // The resulting impedance matrix is dense, and modal impedances must be computed from the matrix.
   // Since matrix calculations are used, it is required that all modes have valid voltages and/or currents.
   } else if (is_line_impedance(fem->projData->solution_impedance_calculation)) {

      // allocate memory that must be freed elsewhere
      // space for an nxn matrix, where n is the number of modes
      if (Zvi) {free(Zvi); Zvi=nullptr;}
      if (Zpv) {free(Zpv); Zpv=nullptr;}
      if (Zpi) {free(Zpi); Zpi=nullptr;}
      Zvi=(complex<double> *) malloc(n*n*sizeof(complex<double>));
      Zpv=(complex<double> *) malloc(n*n*sizeof(complex<double>));
      Zpi=(complex<double> *) malloc(n*n*sizeof(complex<double>));

      // temporary memory and some pointers
      complex<double> *V=Zpv;
      complex<double> *VT=(complex<double> *) malloc(n*n*sizeof(complex<double>));
      complex<double> *I=Zpi;
      complex<double> *Ivi=Zvi;
      complex<double> *IT=(complex<double> *) malloc(n*n*sizeof(complex<double>));
      complex<double> *Ppv=(complex<double> *) malloc(n*n*sizeof(complex<double>));
      complex<double> *Ppi=(complex<double> *) malloc(n*n*sizeof(complex<double>));

      i=0;
      while (i < n*n) {
         Zvi[i]=complex<double>(DBL_MAX,DBL_MAX);
         Zpv[i]=complex<double>(DBL_MAX,DBL_MAX);
         Zpi[i]=complex<double>(DBL_MAX,DBL_MAX);
         i++;
      }

      int mode=-1;
      i=0;
      while (i < (long unsigned int)n) {  // row

         long unsigned int k=0;
         while (k < modeList.size()) {
            if (modeList[k]->get_modeNumber() == i) {mode=(int)k; break;}
            k++;
         }
         if (mode == -1) {
            PetscPrintf(PETSC_COMM_WORLD,"ASSERT: ModeDatabase::calculateImpedanceMatrix did not find a mode number.\n");
            return;
         }

         // use column major format
         long unsigned int j=0;
         while (j < (long unsigned int)n) {  // column
            V[mode+j*n]=modeList[mode]->get_voltage(j);
            VT[mode+j*n]=modeList[mode]->get_voltage(j);

            I[mode+j*n]=modeList[mode]->get_current(j);
            Ivi[mode+j*n]=modeList[mode]->get_current(j);
            IT[mode+j*n]=modeList[mode]->get_current(j);

            if (j == (long unsigned int)mode) Ppv[mode+j*n]=modeList[mode]->get_Pzavg();
            else Ppv[mode+j*n]=complex<double>(0,0);

            if (j == (long unsigned int)mode) Ppi[mode+j*n]=modeList[mode]->get_Pzavg();
            else Ppi[mode+j*n]=complex<double>(0,0);

            j++;
         }

         i++;
      }

      // check for valid voltages and currents
      bool validCurrent=true;
      bool validVoltage=true;
      i=0;
      while (i < n*n) {
         if (I[i] == complex<double>(DBL_MAX,DBL_MAX)) validCurrent=false;
         if (V[i] == complex<double>(DBL_MAX,DBL_MAX)) validVoltage=false;
         i++;
      }

      // normalize the voltages and currents by the power
      lapack_int p,q;
      p=0;
      while (p < (int)n) {
         q=0;
         while (q < (int)n) {
            V[p+q*n]/=sqrt(Ppv[p+p*n]);
            VT[p+q*n]/=sqrt(Ppv[p+p*n]);
            I[p+q*n]/=sqrt(Ppi[p+p*n]);
            Ivi[p+q*n]/=sqrt(Ppi[p+p*n]);
            IT[p+q*n]/=sqrt(Ppi[p+p*n]);
            q++;
         }
         Ppv[p+p*n]=complex<double>(1,0);
         Ppi[p+p*n]=complex<double>(1,0);
         p++;
      }

      // show voltage, current, and power
      if (fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Voltage matrix (V):");
         if (validVoltage) {
            PetscPrintf(PETSC_COMM_WORLD,"\n");
            matrixPrint((lapack_complex_double *)V,n);
         } else PetscPrintf(PETSC_COMM_WORLD," not defined\n");

         PetscPrintf(PETSC_COMM_WORLD,"            Current matrix (I):");
         if (validCurrent) {
            PetscPrintf(PETSC_COMM_WORLD,"\n");
            matrixPrint((lapack_complex_double *)I,n);
         } else PetscPrintf(PETSC_COMM_WORLD," not defined");

         PetscPrintf(PETSC_COMM_WORLD,"            Average propagated power (Pz,avg):\n");
         matrixDiagonalPrint((lapack_complex_double *)Ppv,n);
      }

      // VI definition

      if (validVoltage && validCurrent) {
         matrixInverse((lapack_complex_double *)Ivi,n);
         matrixMultiply((lapack_complex_double *)V,(lapack_complex_double *)Ivi,n);
      } 

      if (strcmp(fem->projData->solution_impedance_definition,"VI") == 0 || fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Impedance matrix (VI):");
         if (validVoltage && validCurrent) {
            PetscPrintf(PETSC_COMM_WORLD,"\n");
            matrixPrint((lapack_complex_double *)Ivi,n);
         } else PetscPrintf(PETSC_COMM_WORLD," not defined\n");
      }

      // PV definition

      if (validVoltage) {
         matrixInverse((lapack_complex_double *)V,n);
         matrixTranspose((lapack_complex_double *)VT,n);
         matrixConjugate((lapack_complex_double *)VT,n);
         matrixInverse((lapack_complex_double *)VT,n);

         matrixTranspose((lapack_complex_double *)Ppv,n);
         matrixConjugate((lapack_complex_double *)Ppv,n);

         matrixMultiply((lapack_complex_double *)VT,(lapack_complex_double *)Ppv,n);
         matrixMultiply((lapack_complex_double *)Ppv,(lapack_complex_double *)V,n);
         matrixScale((lapack_complex_double *)V,(lapack_complex_double *)&two,n);
         matrixInverse((lapack_complex_double *)V,n);
      }

      if (strcmp(fem->projData->solution_impedance_definition,"PV") == 0 || fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Impedance matrix (PV):");
         if (validVoltage) {
            PetscPrintf(PETSC_COMM_WORLD,"\n");
            matrixPrint((lapack_complex_double *)V,n);
         } else PetscPrintf(PETSC_COMM_WORLD," not defined\n");
      }

      // PI definition

      if (validCurrent) {
         matrixInverse((lapack_complex_double *)I,n);
         matrixTranspose((lapack_complex_double *)IT,n);
         matrixConjugate((lapack_complex_double *)IT,n);
         matrixInverse((lapack_complex_double *)IT,n);
         matrixMultiply((lapack_complex_double *)IT,(lapack_complex_double *)Ppi,n);
         matrixMultiply((lapack_complex_double *)Ppi,(lapack_complex_double *)I,n);
         matrixScale((lapack_complex_double *)I,(lapack_complex_double *)&two,n);
      }

      if (strcmp(fem->projData->solution_impedance_definition,"PI") == 0 || fem->projData->debug_show_impedance_details) {
         PetscPrintf(PETSC_COMM_WORLD,"            Impedance matrix (PI):");
         if (validCurrent) {
            PetscPrintf(PETSC_COMM_WORLD,"\n");
            matrixPrint((lapack_complex_double *)I,n);
         } else PetscPrintf(PETSC_COMM_WORLD," not defined\n");
      }

      // clean up
      free(VT);
      free(IT);
      free(Ppv);
      free(Ppi);

   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Invalid impedance calculation selection.\n");
   }
}


void ModeDatabase::calculatePerturbationalLoss(fem2D *fem, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase, MaterialDatabase *materialDatabase)
{
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {
      modeList[i]->calculatePerturbationalLoss(fem,boundaryDatabase,borderDatabase,materialDatabase);
      i++;
   }
}

void ModeDatabase::calculateFieldPoints (fem2D *fem, FieldPointDatabase *fieldPointDatabase)
{
   long unsigned int i=0;
   while (i < (long unsigned int)fem->projData->solution_active_mode_count) {
      modeList[i]->calculateFieldPoints (fem,fieldPointDatabase);
      i++;
   }
}

// sort degenerate modes by power
void ModeDatabase::sort (fem2D *fem)
{
   PetscMPIInt rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   long unsigned int i=0;
   while (i < modeList.size()-1) {
      bool found=true;
      while (found) {
         found=false;
         long unsigned int j=i+1;
         while (j < modeList.size()) {
            complex<double> gamma_i=modeList[i]->get_gamma();
            complex<double> gamma_j=modeList[j]->get_gamma();

            if (! modeList[i]->get_fields()) continue;
            if (! modeList[j]->get_fields()) continue;

            complex<double> scale_i=modeList[i]->get_fields()->get_scale();
            complex<double> scale_j=modeList[j]->get_fields()->get_scale();

            // scale Pzavg by the scale factor used in Fields::loadEigenvector
            // This is an absolute comparison, and the modes use unique scale factors.
            // To get a valid comparison, the scaling has to be reversed.
            // Note that for use in Z calculations, the calculation cancels the scaling,
            // so unique scale factors per mode is ok there.
            double Pzavg_i=abs(modeList[i]->get_Pzavg()*scale_i*scale_i);
            double Pzavg_j=abs(modeList[j]->get_Pzavg()*scale_j*scale_j);

            // check if degenerate
            if (abs((gamma_i-gamma_j)/gamma_i) < 1e-4) {

               // swap if j has more power than i
               if (Pzavg_i < Pzavg_j) {

                  Mode *temp=modeList[i];
                  modeList[i]=modeList[j];
                  modeList[j]=temp;

                  long unsigned int mode=modeList[i]->get_modeNumber();
                  modeList[i]->set_modeNumber(modeList[j]->get_modeNumber());
                  modeList[j]->set_modeNumber(mode);

                  found=true;
               }
            }
            j++;
         }
      }
      i++;
   }
}

void ModeDatabase::print()
{
   PetscPrintf(PETSC_COMM_WORLD,"ModeDatabase:\n");
   PetscPrintf(PETSC_COMM_WORLD,"   this=%p\n",this);
   PetscPrintf(PETSC_COMM_WORLD,"   modeList.size=%ld\n",modeList.size());

   long unsigned int i=0;
   while (i < (long unsigned int)modeList.size()) {
      modeList[i]->print();
      i++;
   }
}

void ModeDatabase::save_grid_E()
{
   long unsigned int i=0;
   while (i < (long unsigned int)modeList.size()) {
      modeList[i]->save_grid_E(i);
      i++;
   }
}

void ModeDatabase::saveAsOne_grid_E()
{
   long unsigned int i=0;
   while (i < (long unsigned int)modeList.size()) {
      modeList[i]->saveAsOne_grid_E(i);
      i++;
   }
}

ModeDatabase::~ModeDatabase()
{
   if (Zvi) {free(Zvi); Zvi=nullptr;}
   if (Zpv) {free(Zpv); Zpv=nullptr;}
   if (Zpi) {free(Zpi); Zpi=nullptr;}

   long unsigned int i=0;
   while (i < modeList.size()) {
      delete modeList[i];
      i++;
   }
}

//--------------------------------------------------------------------------------------------------------------------
// fem2D
//--------------------------------------------------------------------------------------------------------------------

fem2D::fem2D(struct projectData *projData_, ParMesh *pmesh_, int order_, double frequency_, int iteration_,
             PWConstCoefficient *ko2Re_e, PWConstCoefficient *ko2Im_e, PWConstCoefficient *Inv_mu, // for E-field
             PWConstCoefficient *w_mu,                                                             // for H-field
             string temporaryDirectory_) 
{
   projData=projData_;
   pmesh=pmesh_;
   order=order_;
   frequency=frequency_;
   iteration=iteration_;
   temporaryDirectory=temporaryDirectory_;

   // Array for holding border element selections
   border_attributes.SetSize(pmesh->bdr_attributes.Max());

   //*****************************************************************************************************
   // finite element spaces
   //*****************************************************************************************************

   fec_ND=new ND_FECollection(order,pmesh->Dimension());
   fespace_ND=new ParFiniteElementSpace(pmesh,fec_ND);

   fec_H1=new H1_FECollection(order,pmesh->Dimension());
   fespace_H1=new ParFiniteElementSpace(pmesh,fec_H1);

   fec_L2=new L2_FECollection(order-1,pmesh->Dimension());
   fespace_L2=new ParFiniteElementSpace(pmesh,fec_L2);

   //*****************************************************************************************************
   // matrices for the solution of the eigenvalue problem for the E-field and gamma
   //*****************************************************************************************************

   // bilinear forms - Implements the integrals over the finite element spaces.
   // These are defined as (13) in the paper.  The constants from (12) are pulled in so that MFEM can integrate
   // material changes over the cross section.  The dielectric constant can be complex, while the relative permeability
   // is only real.
   //
   // There are 5 equations in (12).
   // Tt (1st instance) - complex eps constant - Done in real and imaginary parts using real math since MFEM doesn't support complex.
   // Tt (2nd instance) - real mur constant
   // G  - real mur constant
   // Sz - real mur constant
   // Tz - complex eps constant - Done in real and imaginary parts using real math since MFEM doesn't support complex.
   // St - real mur constant
   //
   // So there are a total of 8 integrals to calculate.

   // Tt

   // instance with 1/mur

   ParBilinearForm Tt_mur=ParBilinearForm(fespace_ND);
   Tt_mur.AddDomainIntegrator(new VectorFEMassIntegrator (*Inv_mu));
   Tt_mur.Assemble();
   Tt_mur.Finalize();
   saveMat (Tt_mur,"Tt_mur_mat",temporaryDirectory,projData->project_name,false);

   // instance with ko^2*Re(er)

   ParBilinearForm Tt_eps_re=ParBilinearForm(fespace_ND);
   Tt_eps_re.AddDomainIntegrator(new VectorFEMassIntegrator (*ko2Re_e));
   Tt_eps_re.Assemble();
   Tt_eps_re.Finalize();
   saveMat (Tt_eps_re,"Tt_eps_re_mat",temporaryDirectory,projData->project_name,false);

   // instance with ko^2*Im(er)

   ParBilinearForm Tt_eps_im=ParBilinearForm(fespace_ND);
   Tt_eps_im.AddDomainIntegrator(new VectorFEMassIntegrator (*ko2Im_e));
   Tt_eps_im.Assemble();
   Tt_eps_im.Finalize();
   saveMat (Tt_eps_im,"Tt_eps_im_mat",temporaryDirectory,projData->project_name,false);

   // G

   ParMixedBilinearForm G=ParMixedBilinearForm(fespace_H1,fespace_ND);
   G.AddDomainIntegrator(new MixedVectorGradientIntegrator(*Inv_mu));
   G.Assemble();
   G.Finalize();
   saveMat (G,"G_mat",temporaryDirectory,projData->project_name,false);

   // GT - transpose of G

   //ToDo - put back to just being a transpose of G?
   ParMixedBilinearForm GT=ParMixedBilinearForm(fespace_H1,fespace_ND);
   GT.AddDomainIntegrator(new MixedVectorGradientIntegrator(*Inv_mu));
   GT.Assemble();
   GT.Finalize();
   saveMat (GT,"GT_mat",temporaryDirectory,projData->project_name,true);

   // Sz

   ParBilinearForm Sz=ParBilinearForm(fespace_H1);
   Sz.AddDomainIntegrator(new DiffusionIntegrator (*Inv_mu));
   Sz.Assemble();
   Sz.Finalize();
   saveMat (Sz,"Sz_mat",temporaryDirectory,projData->project_name,false);

   // Tz

   // instance with ko^2*Re(er)

   ParBilinearForm Tz_eps_re=ParBilinearForm(fespace_H1);
   Tz_eps_re.AddDomainIntegrator(new MassIntegrator (*ko2Re_e));
   Tz_eps_re.Assemble();
   Tz_eps_re.Finalize();
   saveMat (Tz_eps_re,"Tz_eps_re_mat",temporaryDirectory,projData->project_name,false);

   // instance with ko^2*Im(er)

   ParBilinearForm Tz_eps_im=ParBilinearForm(fespace_H1);
   Tz_eps_im.AddDomainIntegrator(new MassIntegrator (*ko2Im_e));
   Tz_eps_im.Assemble();
   Tz_eps_im.Finalize();
   saveMat (Tz_eps_im,"Tz_eps_im_mat",temporaryDirectory,projData->project_name,false);

   // St

   ParBilinearForm St=ParBilinearForm(fespace_ND);
   St.AddDomainIntegrator(new CurlCurlIntegrator (*Inv_mu));
   St.Assemble();
   St.Finalize();
   saveMat (St,"St_mat",temporaryDirectory,projData->project_name,false);

   //*****************************************************************************************************
   // matrices for the solution of the H-field
   //*****************************************************************************************************

   // Mt

   ParBilinearForm Mt=ParBilinearForm(fespace_ND);
   Mt.AddDomainIntegrator(new VectorFEMassIntegrator(*w_mu));
   Mt.Assemble();
   Mt.Finalize();
   saveMat (Mt,"Mt_mat",temporaryDirectory,projData->project_name,false);

   // for cross product with zHat
   DenseMatrix *m=new DenseMatrix(2,2);
   m->Elem(0,0)=0;
   m->Elem(0,1)=-1;
   m->Elem(1,0)=1;
   m->Elem(1,1)=0;
   MatrixCoefficient *mmat=new MatrixConstantCoefficient (*m);

   // Cz
   ParMixedBilinearForm Cz=ParMixedBilinearForm(fespace_H1,fespace_ND);
   Cz.AddDomainIntegrator(new MixedVectorGradientIntegrator(*mmat));
   Cz.Assemble();
   Cz.Finalize();
   saveMat (Cz,"Cz_mat",temporaryDirectory,projData->project_name,false);

   // Zt

   ParBilinearForm Zt=ParBilinearForm(fespace_ND);
   Zt.AddDomainIntegrator(new VectorFEMassIntegrator(*mmat));
   Zt.Assemble();
   Zt.Finalize();
   saveMat (Zt,"Zt_mat",temporaryDirectory,projData->project_name,false);

   // clean up for zHat x vector
   delete mmat;
   delete m;

   // Mz

   ParBilinearForm Mz=ParBilinearForm(fespace_H1);
   Mz.AddDomainIntegrator(new MassIntegrator(*w_mu));
   Mz.Assemble();
   Mz.Finalize();
   saveMat (Mz,"Mz_mat",temporaryDirectory,projData->project_name,false);

   // Ct

   ParMixedBilinearForm Ct=ParMixedBilinearForm(fespace_ND,fespace_H1);
   Ct.AddDomainIntegrator(new MixedScalarCurlIntegrator());
   Ct.Assemble();
   Ct.Finalize();
   saveMat (Ct,"Ct_mat",temporaryDirectory,projData->project_name,false);

}

// for Et
PetscInt* fem2D::get_ess_tdof_ND(BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase)
{
   int i;
   PetscInt *ess_tdof;

   Array<int> ess_tdof_list_ND;

   // set all to PEC
   // surface impedances are computed with PEC fields, then a perturbational calculation is used to get the loss
   border_attributes=1;

   // loop through the borders to find and eliminate the PMC boundaries
   long unsigned int j=0;
   while (j < borderDatabase->get_size()) {
      Border *border=borderDatabase->get_border(j);
      if (border->has_boundary()) {
         if (boundaryDatabase->get_boundary(border->get_boundary().boundary)->is_perfect_magnetic_conductor()) {
            border_attributes[j-1]=0;  // offset down by 1 to align with MFEM convention
         }
      }
      j++;
   }

   // Get dofs for which border_attributes == 1.
   // The rows and columns of these will be zero'ed out in eigensolve.
   // The natural boundary condition is PMC, so not zeroing out the row and column provides for a PMC boundary.
   // Zeroing out the row and column provides for a PEC boundary.
   fespace_ND->GetEssentialTrueDofs(border_attributes, ess_tdof_list_ND);
   ess_tdof_size_ND=ess_tdof_list_ND.Size();

   // convert to a PetscInt array for passing to the C-language eigensolve

   PetscMalloc(ess_tdof_size_ND*sizeof(PetscInt),&ess_tdof);
   HYPRE_BigInt *offset=fespace_ND->GetTrueDofOffsets();

   i=0;
   while (i < ess_tdof_size_ND) {
      ess_tdof[i]=ess_tdof_list_ND[i]+offset[0];
      i++;
   }

   return ess_tdof;
}

// for Ez
PetscInt* fem2D::get_ess_tdof_H1(BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase)
{
   int i;
   PetscInt *ess_tdof;

   Array<int> ess_tdof_list_H1;

   // set all to PEC
   // surface impedances are computed with PEC fields, then a perturbational calculation is used to get the loss
   border_attributes=1;

   // loop through the borders to find and eliminate the PMC boundaries
   long unsigned int j=0;
   while (j < borderDatabase->get_size()) {
      Border *border=borderDatabase->get_border(j);
      if (border->has_boundary()) {
         if (boundaryDatabase->get_boundary(border->get_boundary().boundary)->is_perfect_magnetic_conductor()) {
            border_attributes[j-1]=0;  // offset down by 1 to align with MFEM convention
         }
      }
      j++;
   }

   // Get dofs for which border_attributes == 1.
   // The rows and columns of these will be zero'ed out in eigensolve.
   // The natural boundary condition is PMC, so not zeroing out the row and column provides for a PMC boundary.
   // Zeroing out the row and column provides for a PEC boundary.
   fespace_H1->GetEssentialTrueDofs(border_attributes, ess_tdof_list_H1);
   ess_tdof_size_H1=ess_tdof_list_H1.Size();

   // convert to a PetscInt array for passing to the C-language eigensolve

   PetscMalloc(ess_tdof_size_H1*sizeof(PetscInt),&ess_tdof);
   HYPRE_BigInt *offset=fespace_H1->GetTrueDofOffsets();

   i=0;
   while (i < ess_tdof_size_H1) {
      ess_tdof[i]=ess_tdof_list_H1[i]+offset[0];
      i++;
   }

   return ess_tdof;
}

bool fem2D::buildFields (double *alpha, double *beta)
{
   return modeDatabase.buildFields(this,alpha,beta);
}

void fem2D::calculateImpedanceMatrix (const char *definition, BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase)
{
   modeDatabase.calculateImpedanceMatrix(this,definition,boundaryDatabase,borderDatabase);
}

void fem2D::calculatePerturbationalLoss(BoundaryDatabase *boundaryDatabase, BorderDatabase *borderDatabase, MaterialDatabase *materialDatabase)
{
   modeDatabase.calculatePerturbationalLoss(this,boundaryDatabase,borderDatabase,materialDatabase);
}

void fem2D::calculateFieldPoints (FieldPointDatabase *fieldPointDatabase)
{
   modeDatabase.calculateFieldPoints(this,fieldPointDatabase);
}

void fem2D::saveParaView()
{
   modeDatabase.saveParaView(this);
}

bool fem2D::ZZrefineMeshes (ConvergenceDatabase *convergenceDatabase)
{
   startingMeshSize=GetGlobalNE(pmesh);

   // see how many modes will be refined and calculate a scale factor
   if (projData->refinement_refine_converged_modes) meshScale=projData->solution_active_mode_count;
   else {
      meshScale=projData->solution_active_mode_count-convergenceDatabase->is_converged_count();
      if (meshScale == 0) meshScale=1;  // should not happen
   }

   PetscPrintf(PETSC_COMM_WORLD,"      Refining mesh ...\n");

   return modeDatabase.ZZrefineMeshes(this,convergenceDatabase);
}

void fem2D::writeInitialGuess ()
{
   modeDatabase.writeInitialGuess(this);
}

Result* fem2D::updateResults(ResultDatabase *resultDatabase, ConvergenceDatabase *convergenceDatabase,
                          chrono::system_clock::time_point solve_start_time, chrono::system_clock::time_point solve_end_time)
{
   // results

   Result *result=new Result();

   result->set_active();
   result->set_iteration(iteration);
   result->set_frequency(frequency);
   result->set_modeCount(projData->solution_active_mode_count);
   result->set_modalImpedanceCalculation(is_modal_impedance(projData->solution_impedance_calculation));

   long unsigned int i=0;
   while (i < modeDatabase.get_size()) {
      if (modeDatabase.get_mode(i)->get_modeNumber() == i) {
         result->push_gamma(complex<double>(modeDatabase.get_mode(i)->get_alpha(),modeDatabase.get_mode(i)->get_beta()));
         result->push_alpha_perturbation(modeDatabase.get_mode(i)->get_perturbationalLoss());
         result->push_Pz(modeDatabase.get_mode(i)->get_Pzavg());
      }
      i++;
   }

   int Zdim=projData->solution_active_mode_count;
   complex<double> *Z=nullptr;
   if (strcmp(projData->solution_impedance_definition,"none") != 0) {
      if (strcmp(projData->solution_impedance_definition,"VI") == 0) Z=modeDatabase.get_Zvi();
      if (strcmp(projData->solution_impedance_definition,"PV") == 0) Z=modeDatabase.get_Zpv();
      if (strcmp(projData->solution_impedance_definition,"PI") == 0) Z=modeDatabase.get_Zpi();

      if (is_modal_impedance(projData->solution_impedance_calculation)) {
         int i=0;
         while (i < Zdim) {
            result->push_Z(Z[i]);
            i++;
         }
      } else {
         int i=0;
         while (i < Zdim*Zdim) {
            result->push_Z(Z[i]);
            i++;
         }
      }
   }

   run_statistics *run_stats=new run_statistics;
   run_stats->set_error(0);
   run_stats->set_meshSize(GetGlobalNE(pmesh));
   run_stats->set_matrixSize(t_size+z_size);
   run_stats->set_start_time(solve_start_time);
   run_stats->set_end_solve_time(solve_end_time);
   run_stats->set_end_refine_time(solve_end_time);  // update later if there is refinement
   run_stats->set_converged(false);

   result->set_run_stats(run_stats);

   resultDatabase->push(result);

   // convergence

   i=0;
   while (i < modeDatabase.get_size()) {
      complex<double> z=complex<double>(DBL_MAX,DBL_MAX);
      if (Z) {
         if (is_modal_impedance(projData->solution_impedance_calculation)) z=Z[i];
         if (is_line_impedance(projData->solution_impedance_calculation)) z=Z[i*Zdim+i];
      }
      convergenceDatabase->push(i,z,modeDatabase.get_mode(i)->get_alpha(),modeDatabase.get_mode(i)->get_beta(),
                                modeDatabase.get_mode(i)->get_perturbationalLoss(),real(modeDatabase.get_mode(i)->get_Pzavg()));
      i++;
   }

   return result;
}

fem2D::~fem2D()
{
   if (fespace_ND) {delete fespace_ND; fespace_ND=nullptr;}
   if (fec_ND) {delete fec_ND; fec_ND=nullptr;}

   if (fespace_H1) {delete fespace_H1; fespace_H1=nullptr;}
   if (fec_H1) {delete fec_H1; fec_H1=nullptr;}

   if (fespace_L2) {delete fespace_L2; fespace_L2=nullptr;}
   if (fec_L2) {delete fec_L2; fec_L2=nullptr;}

}

// last used error is 122
