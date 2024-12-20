////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM2D - A fullwave 2D electromagnetic simulator.                  //
//    Copyright (C) 2024 Brian Young                                          //
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

#ifndef HSOLVE_H
#define HSOLVE_H

#include "eigensolve.h"

int Hsetup (struct projectData *, Mat *, Mat *, Mat *, Mat *, Mat *, PetscMPIInt);
int Hsolve (struct projectData *, Mat *, Mat *, Mat *, Mat *, Mat *, PetscInt, PetscInt, Vec *, PetscScalar *, Vec *, PetscMPIInt);
PetscErrorCode printMatInfo (const char *, Mat *);
void convergence (struct projectData *, KSPConvergedReason);

struct mpi_complex_int {
   double real;
   double imag;
   int location;
};

int VecSplit (Vec *, PetscInt, PetscInt, Vec *);

#endif

