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

#include "Zsolve.h"

double tol=1e-14;

void printComplex (lapack_complex_double a)
{
   double re, im;

   re=lapack_complex_double_real(a);
   im=lapack_complex_double_imag(a);

   if (fabs(re) < 1e-15) re=0;
   if (fabs(im) < 1e-15) im=0;

   PetscPrintf(PETSC_COMM_WORLD,"%#12.6g",re);
   if (im < 0) PetscPrintf(PETSC_COMM_WORLD,"-j%#11.6g",-im);
   else PetscPrintf(PETSC_COMM_WORLD,"+j%#-11.6g",im);
}

void linearPrint (lapack_complex_double *A, lapack_int n)
{
   lapack_int m=0;
   while (m < n*n) {
      printComplex(A[m]);
      PetscPrintf(PETSC_COMM_WORLD,"\n");
      m++;
   }
}

void matrixPrint (lapack_complex_double *A, lapack_int n)
{
   lapack_int i,j;

   PetscPrintf(PETSC_COMM_WORLD,"                   ");
   j=0;
   while (j < n) {
      PetscPrintf(PETSC_COMM_WORLD,"             %3d           ",j+1);
      j++;
   }
   PetscPrintf(PETSC_COMM_WORLD,"\n");

   i=0;
   while (i < n) {  // rows
      PetscPrintf(PETSC_COMM_WORLD,"                 %3d: ",i+1);
      j=0;
      while (j < n) {  // columns
         printComplex(A[i+j*n]);
         PetscPrintf(PETSC_COMM_WORLD,"  ");
         j++;
      }
      PetscPrintf(PETSC_COMM_WORLD,"\n");
      i++;
   }
   PetscPrintf(PETSC_COMM_WORLD,"\n");
}

void matrixDiagonalPrint (lapack_complex_double *A, lapack_int n)
{
   lapack_int i,j;

   PetscPrintf(PETSC_COMM_WORLD,"\n");

   i=0;
   while (i < n) {  // rows
      PetscPrintf(PETSC_COMM_WORLD,"                 %3d: ",i+1);
      j=0;
      while (j < n) {  // columns
         if (i == j) {
            printComplex(A[i+j*n]);
         }
         j++;
      }
      PetscPrintf(PETSC_COMM_WORLD,"\n");
      i++;
   }
   PetscPrintf(PETSC_COMM_WORLD,"\n");
}

int isIdentityMatrix (lapack_complex_double *A, lapack_int n)
{
   lapack_int i,j;

   i=0;
   while (i < n) {  // rows
      j=0;
      while (j < n) {  // columns
         if (i == j) {
            if (cabs(A[i+j*n]-CMPLX(1,0)) > tol) return 1;
         } else {
            if (cabs(A[i+j*n]) > tol) return 1;
         }
         j++;
      }
      i++;
   }

   return 0;
}

int isEqualMatrix (lapack_complex_double *A, lapack_complex_double *B, lapack_int n)
{
   lapack_int i,j;

   i=0;
   while (i < n) {  // rows
      j=0;
      while (j < n) {  // columns
         if (cabs(A[i+j*n]-B[i+j*n]) > tol) return 1;
         j++;
      }
      i++;
   }

   return 0;
}

void matrixTranspose (lapack_complex_double *A, lapack_int n)
{
   lapack_int i,j;
   lapack_complex_double temp;

   i=0;
   while (i < n) {  // rows
      j=i+1;
      while (j < n) {  // columns
         temp=A[i+j*n];
         A[i+j*n]=A[i*n+j];
         A[i*n+j]=temp;
         j++;
      }
      i++;
   }
}

void matrixConjugate (lapack_complex_double *A, lapack_int n)
{
   lapack_int i,j;

   i=0;
   while (i < n) {
      j=0;
      while (j < n) {
         A[i+j*n]=CMPLX(lapack_complex_double_real(A[i+j*n]),
                        lapack_complex_double_imag(A[i+j*n]));
         j++;
      }
      i++;
   }
}

void matrixScale (lapack_complex_double *A, lapack_complex_double *scale, lapack_int n)
{
   lapack_int i,j;

   i=0;
   while (i < n) {
      j=0;
      while (j < n) {
         A[i+j*n]*=(*scale);
         j++;
      }
      i++;
   }
}

int matrixInverse (lapack_complex_double *A, lapack_int n)
{
   if (n == 1) {
      A[0]=1/A[0];
      return 0;
   }

   lapack_int info;
   lapack_int *ipiv = (lapack_int *) malloc (n*n*sizeof(lapack_int));

   info=LAPACKE_zgetrf(LAPACK_COL_MAJOR,n,n,A,n,ipiv);
   if (info > 0) return info;

   info=LAPACKE_zgetri(LAPACK_COL_MAJOR,n,A,n,ipiv);
   if (info > 0) return info;

   free (ipiv);

   return 0;
}

// A*B, result returned in B
void matrixMultiply (lapack_complex_double *A, lapack_complex_double *B, lapack_int N)
{
   lapack_int i,j,n;
   lapack_complex_double *C;

   // temp space
   C=(lapack_complex_double *) malloc(N*N*sizeof(lapack_complex_double));

   i=0;
   while (i < N) {
      j=0;
      while (j < N) {
         C[i+j*N]=CMPLX(0,0);
         j++;
      }
      i++;
   }

   // C=A*B
   i=0;
   while (i < N) {
      j=0;
      while (j < N) {
         n=0;
         while (n < N) {
            C[i+j*N]+=A[i+n*N]*B[n+j*N];
            n++;
         }
         j++;
      }
      i++;
   }

   // copy to B
   i=0;
   while (i < N) {
      j=0;
      while (j < N) {
         B[i+j*N]=C[i+j*N];
         j++;
      }
      i++;
   }

   free(C); 
}

void colMajorTest()
{
   int n=3;
   lapack_complex_double *A;
   A=(lapack_complex_double *) malloc(n*n*sizeof(lapack_complex_double));

   /*
   matrix:
       0 1 2
       3 4 5
       6 7 8
   linear for column major:
       0 3 6 1 4 7 2 5 8
   */

   A[0]=CMPLX(0,0);
   A[1]=CMPLX(3,3);
   A[2]=CMPLX(6,6);
   A[3]=CMPLX(1,1);
   A[4]=CMPLX(4,4);
   A[5]=CMPLX(7,7);
   A[6]=CMPLX(2,2);
   A[7]=CMPLX(5,5);
   A[8]=CMPLX(8,8);

   PetscPrintf(PETSC_COMM_WORLD,"linear:\n");
   linearPrint(A,n);

   PetscPrintf(PETSC_COMM_WORLD,"matrix:\n");
   matrixPrint(A,n);

   free(A);
}

void matrixTest()
{
   lapack_int i,j,n;
   lapack_complex_double *A,*B,*C;

   srand(time(0));
 
   n=30;
   A=(lapack_complex_double *) malloc(n*n*sizeof(lapack_complex_double));
   B=(lapack_complex_double *) malloc(n*n*sizeof(lapack_complex_double));
   C=(lapack_complex_double *) malloc(n*n*sizeof(lapack_complex_double));

   // generate a random test matrix
   // good for n up to about 40
   i=0;
   while (i < n) {
      j=0;
      while (j < n) {
         A[i+j*n]=CMPLX((double)rand()/(double)RAND_MAX-0.5,(double)rand()/(double)RAND_MAX-0.5);
         B[i+j*n]=A[i+j*n];
         C[i+j*n]=A[i+j*n];
         j++;
      }
      i++;
   }

   matrixTranspose(A,n);
   matrixTranspose(A,n);
   if (isEqualMatrix (A,B,n)) PetscPrintf(PETSC_COMM_WORLD,"testMatrix: (A^T)^T==A FAILED\n");
   else PetscPrintf(PETSC_COMM_WORLD,"testMatrix: (A^T)^T==A passed\n");

   if (matrixInverse(A,n)) {PetscPrintf(PETSC_COMM_WORLD,"ERROR2277: matrix inversion error.\n"); return;}
   matrixMultiply(A,B,n);
   if (isIdentityMatrix(B,n)) PetscPrintf(PETSC_COMM_WORLD,"testMatrix: A^-1*A==I FAILED\n");
   else PetscPrintf(PETSC_COMM_WORLD,"testMatrix: A^-1*A==I passed\n");

   if (matrixInverse(A,n)) {PetscPrintf(PETSC_COMM_WORLD,"ERROR2278: matrix inversion error.\n"); return;}
   if (isEqualMatrix (A,C,n)) PetscPrintf(PETSC_COMM_WORLD,"testMatrix: (A^-1)^-1==A FAILED\n");
   else PetscPrintf(PETSC_COMM_WORLD,"testMatrix: (A^-1)^-1==A passed\n");
 

   free(A);
   free(B);
   free(C);
}


