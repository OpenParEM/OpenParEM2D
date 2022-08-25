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

#include "petscsys.h"
#include <iostream>

using namespace std;

void print_copyright_notice ()
{
   PetscPrintf(PETSC_COMM_WORLD,"OpenParEM2D Copyright (C) 2022 Brian Young\n\n");
   PetscPrintf(PETSC_COMM_WORLD,"This program comes with ABSOLUTELY NO WARRANTY.  This is free software,\n");
   PetscPrintf(PETSC_COMM_WORLD,"and you are welcome to distribute it under certain conditions. To see\n");
   PetscPrintf(PETSC_COMM_WORLD,"the full license, set\n");
   PetscPrintf(PETSC_COMM_WORLD,"  output.show.license true\n");
   PetscPrintf(PETSC_COMM_WORLD,"in the setup file.\n\n");
}

void print_license()
{
   cout << "Printing license text ..." << endl;
   cout << endl;

   #include "formatted_license.txt"
}
