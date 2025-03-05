////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM3D - A fullwave 2D electromagnetic simulator.                  //
//    Copyright (C) 2025 Brian Young                                          //
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

#include "Hsolve.h"

PetscErrorCode printMatInfo (const char *mat_name, Mat *mat)
{
   PetscErrorCode ierr=0;
   MatInfo info;
   PetscInt m,n;
   MatType type;

   ierr=MatGetInfo(*mat,MAT_GLOBAL_SUM,&info); CHKERRQ(ierr);
   ierr=MatGetSize(*mat,&m,&n); CHKERRQ(ierr);
   ierr=MatGetType(*mat,&type); CHKERRQ(ierr);

   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s:\n",mat_name);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   rows=%ld, columns=%ld\n",m,n);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   type=%s\n",type);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   block_size=%g\n",info.block_size);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   nz_allocated=%g\n",info.nz_allocated);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   nz_used=%g\n",info.nz_used);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   nz_unneeded=%g\n",info.nz_unneeded);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   memory_allocated=%g\n",info.memory);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   number_of_assemblies=%g\n",info.assemblies);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"   mallocs=%g\n",info.mallocs);

   return ierr;
}

PetscErrorCode convertToIdentity (Mat *mat)
{
   PetscErrorCode ierr=0;
   PetscInt rows,cols;
   int i,j;

   ierr=MatGetSize(*mat,&rows,&cols); CHKERRQ(ierr);
   i=0;
   while (i < rows) {
      j=0;
      while (j < cols) {
         ierr=MatSetValue(*mat, i, j, 0.0, INSERT_VALUES); CHKERRQ(ierr);
         if (i == j) {
            ierr=MatSetValue(*mat, i, j, 1.0, INSERT_VALUES); CHKERRQ(ierr);
         }
         j++;
      }
      i++;
   }

   ierr=MatAssemblyBegin(*mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr=MatAssemblyEnd(*mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

   return ierr;
}

// insert A into B at location rowOffset, colOffset
PetscErrorCode InsertSubMatrix (Mat *A, Mat *B, PetscInt rowOffset, PetscInt colOffset, int show_stats, int show_data)
{
   PetscErrorCode ierr=0;
   PetscInt i,j,nc;
   const PetscInt *aj;
   const PetscScalar *aa;
   PetscInt Arows;
   PetscInt offsetRow[1];
   PetscInt *offsetCol;

   if (show_stats) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"InsertSubMatrix: rowOffset=%ld colOffset=%ld\n",rowOffset,colOffset);
      printMatInfo ("MatInfo on A:",A);
      if (show_data) {ierr=MatView(*A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);}
   }

   // rows
   ierr=MatGetSize(*A,&Arows,NULL); CHKERRQ(ierr);
   i=0;
   while (i < Arows) {
      ierr=MatGetRow(*A,i,&nc,&aj,&aa); CHKERRQ(ierr);

      // offset the columns
      ierr=PetscMalloc(nc*sizeof(PetscInt),&offsetCol); CHKERRQ(ierr);
      j=0;
      while (j < nc) {
         offsetCol[j]=aj[j]+colOffset;
         j++;
      }

      offsetRow[0]=i+rowOffset;
      ierr=MatSetValues(*B,1,offsetRow,nc,offsetCol,aa,INSERT_VALUES); CHKERRQ(ierr);
      ierr=MatRestoreRow(*A,i,&nc,&aj,&aa); CHKERRQ(ierr);
      ierr=PetscFree(offsetCol); CHKERRQ(ierr);
      i++;
   }

   if (show_stats) {
      printMatInfo ("MatInfo on B:",B);
      if (show_data) {ierr=MatView(*B,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);}
   }

   return (ierr);
}

// the basic matrices are only dependent on frequency
int Hsetup (struct projectData *projData, Mat *Mt, Mat *Cz, Mat *Zt, Mat *Mz, Mat *Ct, PetscMPIInt rank)
{
   PetscErrorCode ierr=0;
   PetscInt ioffset,joffset;
   int location,transpose,sign;
   int iMt,iCz,iZt,iMz,iCt;
   char *resultsDir;
   long int MtHeight,MtWidth,MtSparseWidth;
   long int CzHeight,CzWidth,CzSparseWidth;
   long int ZtHeight,ZtWidth,ZtSparseWidth;
   long int MzHeight,MzWidth,MzSparseWidth;
   long int CtHeight,CtWidth,CtSparseWidth;
   int fail=0;

   if (projData->output_show_postprocessing) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"            setting up for H field calculation ...\n");}

   // directory to hold results
   resultsDir=(char *) malloc ((5+strlen(projData->project_name)+1+1)*sizeof(char));
   sprintf (resultsDir,"temp_%s/",projData->project_name);  // align this with OpenParEM2D.cpp, eigensolve.c, Hsolve.c, and Fields::saveFields

   // run through the files to get counts for allocating matrices

   if (rank == 0) {
      if (loadDataFileStats("Mt_mat",resultsDir,projData->project_name,&MtHeight,&MtWidth,&MtSparseWidth) != 0) {
         prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2151: Failed to scan \"Mt_mat\".\n");
         fail=1;
      }

      if (loadDataFileStats("Cz_mat",resultsDir,projData->project_name,&CzHeight,&CzWidth,&CzSparseWidth) != 0) {
         prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2152: Failed to scan \"Cz_mat\".\n");
         fail=1;
      }

      if (loadDataFileStats("Zt_mat",resultsDir,projData->project_name,&ZtHeight,&ZtWidth,&ZtSparseWidth) != 0) {
         prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2153: Failed to scan \"Zt_mat\".\n");
         fail=1;
      }

      if (loadDataFileStats("Mz_mat",resultsDir,projData->project_name,&MzHeight,&MzWidth,&MzSparseWidth) != 0) {
         prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2154: Failed to scan \"Mz_mat\".\n");
         fail=1;
      }

      if (loadDataFileStats("Ct_mat",resultsDir,projData->project_name,&CtHeight,&CtWidth,&CtSparseWidth) != 0) {
         prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2155: Failed to scan \"Ct_mat\".\n");
         fail=1;
      }
   }

   if (MPI_Bcast(&MtHeight,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&MtWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&MtSparseWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;

   if (MPI_Bcast(&CzHeight,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&CzWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&CzSparseWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;

   if (MPI_Bcast(&ZtHeight,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&ZtWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&ZtSparseWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;

   if (MPI_Bcast(&MzHeight,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&MzWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&MzSparseWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;

   if (MPI_Bcast(&CtHeight,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&CtWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;
   if (MPI_Bcast(&CtSparseWidth,1,MPI_LONG,0,PETSC_COMM_WORLD)) fail=1;

   if (MPI_Barrier(PETSC_COMM_WORLD)) fail=1;

   if (fail) goto EXIT;
   //prefix(); PetscPrintf (PETSC_COMM_WORLD,"MtHeight=%zu  MtWidth=%zu  MtSparseWidth=%zu\n",MtHeight,MtWidth,MtSparseWidth);
   //prefix(); PetscPrintf (PETSC_COMM_WORLD,"CzHeight=%zu  CzWidth=%zu  CzSparseWidth=%zu\n",CzHeight,CzWidth,CzSparseWidth);
   //prefix(); PetscPrintf (PETSC_COMM_WORLD,"ZtHeight=%zu  ZtWidth=%zu  ZtSparseWidth=%zu\n",ZtHeight,ZtWidth,ZtSparseWidth);
   //prefix(); PetscPrintf (PETSC_COMM_WORLD,"MzHeight=%zu  MzWidth=%zu  MzSparseWidth=%zu\n",MzHeight,MzWidth,MzSparseWidth);
   //prefix(); PetscPrintf (PETSC_COMM_WORLD,"CtHeight=%zu  CtWidth=%zu  CtSparseWidth=%zu\n",CtHeight,CtWidth,CtSparseWidth);
   //prefix(); PetscPrintf (PETSC_COMM_WORLD,"matrix size: Ht=%zu, Hz=%zu\n",MtHeight,MzHeight);

   // Build up the matrics for the standard linear problem Ax=b.

   // Mt
   ierr=MatCreate(PETSC_COMM_WORLD,Mt); CHKERRQ(ierr);
   ierr=MatSetSizes(*Mt,PETSC_DECIDE,PETSC_DECIDE,MtHeight,MtWidth); CHKERRQ(ierr);
   ierr=MatSetFromOptions(*Mt); CHKERRQ(ierr);
   ierr=MatSeqAIJSetPreallocation(*Mt,MtSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatMPIAIJSetPreallocation(*Mt,MtSparseWidth,NULL,MtSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatZeroEntries(*Mt); CHKERRQ(ierr);

   ioffset=0; joffset=0; location=1; transpose=0; sign=0;
   iMt=loadDataFile ("Mt_mat",resultsDir,projData->project_name,Mt,ioffset,joffset,location,transpose,sign,rank);
   if (iMt) {prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2156: Failed to load \"Mt_mat\" data file.\n"); fail=1;}

   ierr=MatAssemblyBegin(*Mt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr=MatAssemblyEnd(*Mt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

   // Cz
   ierr=MatCreate(PETSC_COMM_WORLD,Cz); CHKERRQ(ierr);
   ierr=MatSetSizes(*Cz,PETSC_DECIDE,PETSC_DECIDE,MtHeight,CzWidth); CHKERRQ(ierr);
   ierr=MatSetFromOptions(*Cz); CHKERRQ(ierr);
   ierr=MatSeqAIJSetPreallocation(*Cz,CzSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatMPIAIJSetPreallocation(*Cz,CzSparseWidth,NULL,CzSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatZeroEntries(*Cz); CHKERRQ(ierr);

   ioffset=0; joffset=0; location=0; transpose=0; sign=0;
   iCz=loadDataFile ("Cz_mat",resultsDir,projData->project_name,Cz,ioffset,joffset,location,transpose,sign,rank);
   if (iCz) {prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2157: Failed to load \"Cz_mat\" data file.\n"); fail=1;}

   ierr=MatAssemblyBegin(*Cz,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr=MatAssemblyEnd(*Cz,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

   // Zt
   ierr=MatCreate(PETSC_COMM_WORLD,Zt); CHKERRQ(ierr);
   ierr=MatSetSizes(*Zt,PETSC_DECIDE,PETSC_DECIDE,ZtHeight,ZtWidth); CHKERRQ(ierr);
   ierr=MatSetFromOptions(*Zt); CHKERRQ(ierr);
   ierr=MatSeqAIJSetPreallocation(*Zt,ZtSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatMPIAIJSetPreallocation(*Zt,ZtSparseWidth,NULL,ZtSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatZeroEntries(*Zt); CHKERRQ(ierr);

   ioffset=0; joffset=0; location=0; transpose=0; sign=0;
   iZt=loadDataFile ("Zt_mat",resultsDir,projData->project_name,Zt,ioffset,joffset,location,transpose,sign,rank);
   if (iZt) {prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2158: Failed to load \"Zt_mat\" data file.\n"); fail=1;}

   ierr=MatAssemblyBegin(*Zt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr=MatAssemblyEnd(*Zt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

   // Mz
   ierr=MatCreate(PETSC_COMM_WORLD,Mz); CHKERRQ(ierr);
   ierr=MatSetSizes(*Mz,PETSC_DECIDE,PETSC_DECIDE,MzHeight,MzWidth); CHKERRQ(ierr);
   ierr=MatSetFromOptions(*Mz); CHKERRQ(ierr);
   ierr=MatSeqAIJSetPreallocation(*Mz,MzSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatMPIAIJSetPreallocation(*Mz,MzSparseWidth,NULL,MzSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatZeroEntries(*Mz); CHKERRQ(ierr);

   ioffset=0; joffset=0; location=1; transpose=0; sign=0;
   iMz=loadDataFile ("Mz_mat",resultsDir,projData->project_name,Mz,ioffset,joffset,location,transpose,sign,rank);
   if (iMz) {prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2159: Failed to load \"Mz_mat\" data file.\n"); fail=1;}

   ierr=MatAssemblyBegin(*Mz,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr=MatAssemblyEnd(*Mz,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

   // Ct
   ierr=MatCreate(PETSC_COMM_WORLD,Ct); CHKERRQ(ierr);
   ierr=MatSetSizes(*Ct,PETSC_DECIDE,PETSC_DECIDE,CtHeight,CtWidth); CHKERRQ(ierr);
   ierr=MatSetFromOptions(*Ct); CHKERRQ(ierr);
   ierr=MatSeqAIJSetPreallocation(*Ct,CtSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatMPIAIJSetPreallocation(*Ct,CtSparseWidth,NULL,CtSparseWidth,NULL); CHKERRQ(ierr);
   ierr=MatZeroEntries(*Ct); CHKERRQ(ierr);

   ioffset=0; joffset=0; location=0; transpose=0; sign=1;
   iCt=loadDataFile ("Ct_mat",resultsDir,projData->project_name,Ct,ioffset,joffset,location,transpose,sign,rank);
   if (iCt) {prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2160: Failed to load \"Ct_mat\" data file.\n"); fail=1;}
   ierr=MatAssemblyBegin(*Ct,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr=MatAssemblyEnd(*Ct,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

   EXIT:

   free(resultsDir);

   return fail;
}

void convergence (struct projectData *projData, KSPConvergedReason reason)
{
   if (projData->output_show_postprocessing) {
      if (reason < 0) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"                  NOT CONVERGED: ");}
      else {prefix(); PetscPrintf(PETSC_COMM_WORLD,"                  Converged: ");}

      if (reason == KSP_CONVERGED_ITERATING) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"CONVERGED_ITERATING");}
      if (reason == KSP_CONVERGED_RTOL_NORMAL) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"RTOL_NORMAL");}
      if (reason == KSP_CONVERGED_ATOL_NORMAL) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"ATOL_NORMAL");}
      if (reason == KSP_CONVERGED_RTOL) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"RTOL");}
      if (reason == KSP_CONVERGED_ATOL) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"ATOL");}
      if (reason == KSP_CONVERGED_ITS) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"ITS");}
      if (reason == KSP_CONVERGED_NEG_CURVE) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"NEG_CURVE");}
      if (reason == KSP_CONVERGED_STEP_LENGTH) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"STEP_LENGTH");}
      if (reason == KSP_CONVERGED_HAPPY_BREAKDOWN) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"HAPPY_BREAKDOWN");}
      if (reason == KSP_DIVERGED_NULL) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"NULL");}
      if (reason == KSP_DIVERGED_ITS) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"ITS");}
      if (reason == KSP_DIVERGED_DTOL) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"DTOL");}
      if (reason == KSP_DIVERGED_BREAKDOWN) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"BREAKDOWN - Generic breakdown during solution.");}
      if (reason == KSP_DIVERGED_BREAKDOWN_BICG) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"BREAKDOWN_BICG");}
      if (reason == KSP_DIVERGED_NONSYMMETRIC) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"NONSYMMETRIC");}
      if (reason == KSP_DIVERGED_INDEFINITE_PC) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"INDEFINITE_PC");}
      if (reason == KSP_DIVERGED_NANORINF) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"NANORINF");}
      if (reason == KSP_DIVERGED_INDEFINITE_MAT) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"INDEFINITE_MAT");}
      if (reason == KSP_DIVERGED_PC_FAILED) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"PC_FAILED - Could not build or use the requested preconditioner.");}
   } else {
      if (reason < 0) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"           Hfield NOT CONVERGED");}
   }
}

// split out a subvector of A from low to high-1 and return in B starting at index 0
// B is allocated and must be destroyed elsewhere
int VecSplit (Vec *A, PetscInt low, PetscInt high, Vec *B)
{
   PetscErrorCode ierr=0;
   PetscInt i;
   PetscInt *idx_from,low_from,high_from,count_from;
   PetscInt *idx_to,low_to,high_to,count_to;
   PetscScalar *data_from,*data_to;
   struct mpi_complex_int *sendData,*recvData;
   PetscMPIInt size,rank;

   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // set up a datatype for data transfer
   // use an int for an array location then two doubles for a complex value

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

   // space to hold transfer lengths and displacements
   int *counts_recv,*displacements_recv;
   counts_recv=(int *)malloc(size*sizeof(int));
   displacements_recv=(int *)malloc(size*sizeof(int));

   // data from A
   ierr=VecGetOwnershipRange(*A,&low_from,&high_from); CHKERRQ(ierr);
   ierr=PetscMalloc((high-low)*sizeof(PetscInt),&idx_from); CHKERRQ(ierr);
   ierr=PetscMalloc((high-low)*sizeof(PetscScalar),&data_from); CHKERRQ(ierr);
   count_from=0;
   i=low;
   while (i < high) {
      if (i >= low_from && i < high_from) {
         idx_from[count_from]=i;
         count_from++;
      }
      i++;
   }

   // gather the transfer lengths
   if (MPI_Gather(&count_from,1,MPI_INT,counts_recv,1,MPI_INT,0,PETSC_COMM_WORLD)) {
      prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2161: Failed to gather data.\n");
      return 1;
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

   // get the data
   ierr=VecGetValues(*A,count_from,idx_from,data_from); CHKERRQ(ierr);

   // assemble the data for transfer
   ierr=PetscMalloc(count_from*sizeof(struct mpi_complex_int),&sendData); CHKERRQ(ierr);
   i=0;
   while (i < count_from) {
      sendData[i].location=idx_from[i];
      sendData[i].real=PetscRealPart(data_from[i]);
      sendData[i].imag=PetscImaginaryPart(data_from[i]);
      i++;
   }
   ierr=PetscFree(data_from); CHKERRQ(ierr);
   ierr=PetscFree(idx_from); CHKERRQ(ierr);

   // space for the received data
   ierr=PetscMalloc((high-low)*sizeof(struct mpi_complex_int),&recvData); CHKERRQ(ierr);

   // send all to rank 0
   if (MPI_Gatherv(sendData,count_from,mpi_complex_int_type,recvData,counts_recv,displacements_recv,mpi_complex_int_type,0,PETSC_COMM_WORLD)) {
      prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2162: Failed to gather data.\n");
      return 1;
   }

   if (counts_recv) {free(counts_recv); counts_recv=NULL;}
   if (displacements_recv) {free(displacements_recv); displacements_recv=NULL;}
   ierr=PetscFree(sendData); CHKERRQ(ierr); 

   // broadcast recvData
   if (MPI_Bcast(recvData,high-low,mpi_complex_int_type,0,PETSC_COMM_WORLD)) {
      prefix(); PetscPrintf (PETSC_COMM_WORLD,"ERROR2163: Failed to broadcast data.\n");
      return 1;
   }

   // transfer data to B
   ierr=VecCreate(PETSC_COMM_WORLD,B); CHKERRQ(ierr);
   ierr=VecSetType(*B,VECSTANDARD); CHKERRQ(ierr);
   ierr=VecSetSizes(*B,PETSC_DECIDE,high-low); CHKERRQ(ierr);

   ierr=VecGetOwnershipRange(*B,&low_to,&high_to); CHKERRQ(ierr);
   ierr=PetscMalloc((high_to-low_to)*sizeof(PetscInt),&idx_to); CHKERRQ(ierr);
   ierr=PetscMalloc((high_to-low_to)*sizeof(PetscScalar),&data_to); CHKERRQ(ierr);

   // setup for transfer
   count_to=0;
   i=0;
   while (i < high-low) {
      if (recvData[i].location-low >= low_to && recvData[i].location-low < high_to) {
         idx_to[count_to]=recvData[i].location-low;
         data_to[count_to]=recvData[i].real+PETSC_i*recvData[i].imag;
         count_to++;
      }
      i++;
   }

   // transfer
   ierr=VecSetValues(*B,count_to,idx_to,data_to,INSERT_VALUES); CHKERRQ(ierr);
   ierr=VecAssemblyBegin(*B); CHKERRQ(ierr);
   ierr=VecAssemblyEnd(*B); CHKERRQ(ierr);

   // cleanup
   ierr=PetscFree(recvData); CHKERRQ(ierr);
   ierr=PetscFree(idx_to); CHKERRQ(ierr);
   ierr=PetscFree(data_to); CHKERRQ(ierr);
   MPI_Type_free(&mpi_complex_int_type);

   return 0;
}

// Solve Ax=b
int Hsolve (struct projectData *projData, Mat *Mt, Mat *Cz, Mat *Zt, Mat *Mz, Mat *Ct,
            PetscInt EtSize, PetscInt EzSize, Vec *Efield, PetscScalar *gamma, Vec *Hfield, PetscMPIInt rank)
{
   PetscErrorCode ierr=0;
   Mat Ztgamma;
   Vec Et,Ez;
   Vec Ht,Hz;
   Vec btt,btz,bz;
   KSP ksp;
   PetscInt its,maxits;
   PetscReal rtol,atol,dtol;
   KSPConvergedReason reason;
   //PC pc;
   int fail=0;

   if (projData->output_show_postprocessing) {prefix(); PetscPrintf(PETSC_COMM_WORLD,"               solving for H field ...\n");}

   // space for Ht and Hz

   ierr=VecCreate(PETSC_COMM_WORLD,&Ht); CHKERRQ(ierr);
   ierr=VecSetType(Ht,VECSTANDARD); CHKERRQ(ierr);
   ierr=VecSetSizes(Ht,PETSC_DECIDE,EtSize); CHKERRQ(ierr);

   ierr=VecCreate(PETSC_COMM_WORLD,&Hz); CHKERRQ(ierr);
   ierr=VecSetType(Hz,VECSTANDARD); CHKERRQ(ierr);
   ierr=VecSetSizes(Hz,PETSC_DECIDE,EzSize); CHKERRQ(ierr);

   // split out Et from Efield
   VecSplit(Efield,0,EtSize,&Et);

   // split out Ez from Efield
   VecSplit(Efield,EtSize,EtSize+EzSize,&Ez);

   //*******************************************************************************************
   // Ht
   //*******************************************************************************************

   // copy Zt and scale by gamma
   ierr=MatDuplicate(*Zt, MAT_COPY_VALUES, &Ztgamma); CHKERRQ(ierr);
   ierr=MatScale(Ztgamma, *gamma); CHKERRQ(ierr);

   // calculate btt

   ierr=VecCreate(PETSC_COMM_WORLD,&btt); CHKERRQ(ierr);
   ierr=VecSetType(btt,VECSTANDARD); CHKERRQ(ierr);
   ierr=VecSetSizes(btt,PETSC_DECIDE,EtSize); CHKERRQ(ierr);

   ierr=MatMult(Ztgamma, Et, btt);

   ierr=VecAssemblyBegin(btt); CHKERRQ(ierr);
   ierr=VecAssemblyEnd(btt); CHKERRQ(ierr);

   // calculate btz

   ierr=VecCreate(PETSC_COMM_WORLD,&btz); CHKERRQ(ierr);
   ierr=VecSetType(btz,VECSTANDARD); CHKERRQ(ierr);
   ierr=VecSetSizes(btz,PETSC_DECIDE,EtSize); CHKERRQ(ierr);

   ierr=MatMult(*Cz, Ez, btz); CHKERRQ(ierr);

   ierr=VecAssemblyBegin(btz); CHKERRQ(ierr);
   ierr=VecAssemblyEnd(btz); CHKERRQ(ierr);

   // get the sum
   VecAXPY(btt,1,btz);

   // solve Ax=b

   // setup
   ierr=KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
   //ierr=KSPSetType(ksp,KSPLSQR); CHKERRQ(ierr);
   ierr=KSPSetOperators(ksp, *Mt, *Mt); CHKERRQ(ierr);
   //ierr=KSPGetPC(ksp,&pc); CHKERRQ(ierr);
   //ierr=PCSetType(pc,PCNONE); CHKERRQ(ierr);
   ierr=KSPSetFromOptions(ksp); CHKERRQ(ierr);

   maxits=projData->solution_iteration_limit;       // PETSc default = 1e4
   rtol=projData->solution_tolerance;               // PETSc default = 1e-5
   //atol=projData->solution_tolerance;             // PETSC default = 1e-50
   atol=1e-50;                                      // PETSC default = 1e-50
   dtol=1e5;                                        // PETSC default = 1e5
   ierr=KSPSetTolerances(ksp, rtol, atol, dtol, maxits); CHKERRQ(ierr);

   // solve
   ierr=KSPSolve(ksp, btt, Ht); CHKERRQ(ierr);

   // get stats
   ierr=KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
   ierr=KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);

   convergence (projData,reason);
   if (projData->output_show_postprocessing) {prefix(); ierr=PetscPrintf (PETSC_COMM_WORLD,", number of Ht iterations: %ld\n",its); CHKERRQ(ierr);}

   // cleanup
   ierr=MatDestroy(&Ztgamma); CHKERRQ(ierr);
   ierr=VecDestroy(&btt); CHKERRQ(ierr);
   ierr=VecDestroy(&btz); CHKERRQ(ierr);
   ierr=VecDestroy(&Ez); CHKERRQ(ierr);
   ierr=KSPDestroy(&ksp); CHKERRQ(ierr);

   //*******************************************************************************************
   // Hz
   //*******************************************************************************************

   // calculate bz

   ierr=VecCreate(PETSC_COMM_WORLD,&bz); CHKERRQ(ierr);
   ierr=VecSetType(bz,VECSTANDARD); CHKERRQ(ierr);
   ierr=VecSetSizes(bz,PETSC_DECIDE,EzSize); CHKERRQ(ierr);

   ierr=MatMult(*Ct, Et, bz);

   ierr=VecAssemblyBegin(bz); CHKERRQ(ierr);
   ierr=VecAssemblyEnd(bz); CHKERRQ(ierr);

   // solve Ax=b

   // setup
   ierr=KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
   //ierr=KSPSetType(ksp,KSPLSQR); CHKERRQ(ierr);
   ierr=KSPSetOperators(ksp, *Mz, *Mz); CHKERRQ(ierr);
   //ierr=KSPGetPC(ksp,&pc); CHKERRQ(ierr);
   //ierr=PCSetType(pc,PCNONE); CHKERRQ(ierr);
   ierr=KSPSetFromOptions(ksp); CHKERRQ(ierr);

   maxits=projData->solution_iteration_limit;       // PETSc default = 1e4
   rtol=projData->solution_tolerance;               // PETSc default = 1e-5
   //atol=projData->solution_tolerance;             // PETSC default = 1e-50
   atol=1e-50;                                      // PETSC default = 1e-50
   dtol=1e5;                                        // PETSC default = 1e5
   ierr=KSPSetTolerances(ksp, rtol, atol, dtol, maxits); CHKERRQ(ierr);

   // solve
   ierr=KSPSolve(ksp, bz, Hz); CHKERRQ(ierr);

   // get stats
   ierr=KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
   ierr=KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);

   convergence (projData,reason);
   if (projData->output_show_postprocessing) {prefix(); ierr=PetscPrintf (PETSC_COMM_WORLD,", number of Hz iterations: %ld\n",its); CHKERRQ(ierr);}

   // cleanup
   ierr=VecDestroy(&bz); CHKERRQ(ierr);
   ierr=KSPDestroy(&ksp); CHKERRQ(ierr);

   //*******************************************************************************************
   // combine Ht and Hz into one vector
   //*******************************************************************************************

   Vec vecList[2];
   vecList[0]=Ht;
   vecList[1]=Hz;
   ierr=VecConcatenate(2,vecList,Hfield,NULL);  // allocates memory

   // cleanup
   ierr=VecDestroy(&Et); CHKERRQ(ierr);
   ierr=VecDestroy(&Ez); CHKERRQ(ierr);
   ierr=VecDestroy(&Ht); CHKERRQ(ierr);
   ierr=VecDestroy(&Hz); CHKERRQ(ierr);

   return fail;
}

