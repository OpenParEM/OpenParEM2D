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

// Class for reading the $PhysicalNames block from Gmsh mesh files format 2.2

#ifndef MESH_H
#define MESH_H

#include "petscsys.h"
#include <slepceps.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

bool is_comment (string);
void split_on_space (vector<string> *, string);

class meshMaterialList {
    private:
       string GMSH_version_number="2.2";
       string regionsFile_version_number="1.0";
       int file_type;
       int data_size;
       vector<int> index;
       vector<string> list;
    public:
       int loadGMSH (const char *);
       int loadMFEM (const char *);
       void print ();
       int size();
       int get_index (long unsigned int);
       string get_name (long unsigned int);
};

#endif

