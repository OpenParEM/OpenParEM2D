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

#include "mesh.hpp"

bool is_comment (string a) {
   long unsigned int i;

   // blank line
   if (a.compare("") == 0) return false;

   // first /
   i=0;
   while (i < a.length()) {
      if ((a.data())[i] == ' ' || (a.data())[i] == '\t') {
         // do nothing
      } else {
         if ((a.data())[i] == '/') {i++; break;}
         else return false;
      }
      i++;
   }

   // second /
   if (i < a.length()) {
      if ((a.data())[i] == '/') return true;
   }

   return false;
}

void split_on_space (vector<string> *tokens, string a) {
   string buf;
   stringstream ss(a);

   tokens->clear();
   while (ss >> buf) tokens->push_back(buf);
}

// parse the msh file for the information in the $PhysicalNames block
int meshMaterialList::loadGMSH (const char *filename)
{
   long unsigned int materialCount;
   int lineCount=0;
   bool startedFormat=false,completedFormat=false,loadedFormat=false;
   bool startedNames=false,completedNames=false,loadedEntryCount=false;
   string line,version_number;
   vector<string> tokens;
   size_t pos1,pos2;

   // the mesh file in Gmsh 2.2 format 
   ifstream meshFile;
   meshFile.open(filename,ifstream::in);

   if (meshFile.is_open()) {

      while (getline(meshFile,line)) {
         lineCount++;

         // skip blank lines
         if (line.compare("") != 0) {

            // skip comment lines
            if (! is_comment(line)) {

               // chop off comments
               line=line.substr(0,line.find("//",0));

               if (startedFormat) {
                  if (line.compare("$EndMeshFormat") == 0) {
                     startedFormat=false; completedFormat=true;

                     if (GMSH_version_number.compare(version_number) != 0) {
                        PetscPrintf(PETSC_COMM_WORLD,"ERROR300: Incorrect mesh format of %s in file \"%s\".\n",version_number.c_str(),filename);
                        meshFile.close();
                        return 1;
                     }

                     if (file_type != 0) {
                        PetscPrintf(PETSC_COMM_WORLD,"ERROR301: Binary file form is not supported for file \"%s\".\n",filename);
                        meshFile.close();
                        return 1;
                     }
                  } else {
                     if (!loadedFormat) {
                        split_on_space (&tokens,line);
                        if (tokens.size() != 3) {
                           PetscPrintf(PETSC_COMM_WORLD,"ERROR302: Incorrect number of tokens in file \"%s\" at line %d.\n",filename,lineCount);
                           meshFile.close();
                           return 1;
                        }

                        version_number=tokens[0];
                        file_type=stoi(tokens[1]);
                        data_size=stoi(tokens[2]);

                        tokens.clear();
                        loadedFormat=true;
                     }
                  }
               }

               if (startedNames) {
                  if (line.compare("$EndPhysicalNames") == 0) {
                     startedNames=false; completedNames=true;

                     if (materialCount != list.size()) {
                        PetscPrintf(PETSC_COMM_WORLD,"ERROR303: Incorrect format in $PhysicalNames block in file \"%s\" at line %d\n",filename,lineCount);
                        meshFile.close();
                        return 1;
                     }
                  } else {

                     if (! loadedEntryCount) {
                        loadedEntryCount=true;
                     } else {
                        split_on_space (&tokens,line);
                        if (tokens.size() != 3) {
                           PetscPrintf(PETSC_COMM_WORLD,"ERROR304: Incorrect number of tokens in file \"%s\" at line %d.\n",filename,lineCount);
                           meshFile.close();
                           return 1;
                        }
                        int dim=stoi(tokens[0]);
                        if (dim != 2) {
                           PetscPrintf(PETSC_COMM_WORLD,"Warning: Dimension %d!=2 in file \"%s\" at line %d.\n",dim,filename,lineCount);
                        }

                        index.push_back(stoi(tokens[1])-1);

                        // strip off "
                        pos1=tokens[2].find("\"",0);
                        if (pos1 >= 0) {
                           pos2=tokens[2].rfind("\"",tokens[2].length());
                           tokens[2]=tokens[2].substr(pos1+1,pos2-pos1-1);
                        }

                        list.push_back(tokens[2]);
                        tokens.clear();
                        materialCount++;
                     }
                  }
               }

               if (line.compare("$MeshFormat") == 0) startedFormat=true;
               if (line.compare("$PhysicalNames") == 0) {startedNames=true; materialCount=0;}
            }
         }
      }

      if (meshFile.bad()) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR305: Error while reading file \"%s\".\n",filename);
         meshFile.close();
         return 1;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR306: File \"%s\" is not available for reading.\n",filename);
      return 1;
   }

   int retval=0;

   if (! completedFormat) {
      if (startedFormat) PetscPrintf(PETSC_COMM_WORLD,"ERROR307: $MeshFormat block is missing the $EndMeshFormat statement in file \"%s\".\n",filename);
      else PetscPrintf(PETSC_COMM_WORLD,"ERROR308: $EndMeshFormat block is missing in file \"%s\".\n",filename);
      retval=1;
   }

   if (! completedNames) {
      if (startedNames) PetscPrintf(PETSC_COMM_WORLD,"ERROR309: $PhysicalNames block is missing the $EndPhysicalNames statement in file \"%s\".\n",filename);
      else PetscPrintf(PETSC_COMM_WORLD,"ERROR310: $PhysicalNames block is missing in file \"%s\".\n",filename);
      retval=1;
   }

   meshFile.close();

   return retval;
}

int meshMaterialList::size()
{
   return list.size();
}

int meshMaterialList::get_index(long unsigned int m)
{
   if (m >= 0 && m < index.size()) return index[m];
   return -1;
}

string meshMaterialList::get_name(long unsigned int m)
{
   string a="ERROR311: out of bounds";
   if (m >= 0 && m < list.size()) return list[m];
   return a;
}

void meshMaterialList::print ()
{
   PetscPrintf(PETSC_COMM_WORLD,"meshMaterialList:\n");
   PetscPrintf(PETSC_COMM_WORLD,"GMSH_version_number=%s\n",GMSH_version_number.c_str());
   PetscPrintf(PETSC_COMM_WORLD,"regionsFile_version_number=%s\n",regionsFile_version_number.c_str());
   PetscPrintf(PETSC_COMM_WORLD,"file-type=%d\n",file_type);
   PetscPrintf(PETSC_COMM_WORLD,"data-size=%d\n",data_size);
   
   long unsigned int i=0;
   while (i < list.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d %s\n",index[i],list[i].c_str());
      i++;
   }
}


int meshMaterialList::loadMFEM (const char *filename)
{
   int lineCount=0;
   string line;
   vector<string> tokens;
   size_t pos1,pos2;

   ifstream regionsFile;
   regionsFile.open(filename,ifstream::in);

   if (regionsFile.is_open()) {

      stringstream ss;
      ss << "openEM Regions v" << regionsFile_version_number;

      // first line is the version
      getline(regionsFile,line);
      lineCount++;
      if (line.compare(ss.str()) != 0) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR312: Incorrect regions format of \"%s\" in file \"%s\" at line 1.\n",ss.str().c_str(),filename);
         regionsFile.close();
         return 1;
      }

      // parse the rest
      while (getline(regionsFile,line)) {
         lineCount++;

         // skip blank lines
         if (line.compare("") != 0) {

            // skip comment lines
            if (! is_comment(line)) {

               // chop off comments
               line=line.substr(0,line.find("//",0));

               split_on_space (&tokens,line);
               if (tokens.size() != 2) {
                  PetscPrintf(PETSC_COMM_WORLD,"ERROR313: Incorrect number of tokens in file \"%s\" at line %d.\n",filename,lineCount);
                  regionsFile.close();
                  return 1;
               }
               index.push_back(stoi(tokens[0])-1);

               // strip off "
               pos1=tokens[1].find("\"",0);
               if (pos1 >= 0) {
                  pos2=tokens[1].rfind("\"",tokens[1].length());
                  tokens[1]=tokens[1].substr(pos1+1,pos2-pos1-1);
               }

               list.push_back(tokens[1]);

               tokens.clear();
            }
         }
      }

      if (regionsFile.bad()) {
         PetscPrintf(PETSC_COMM_WORLD,"ERROR314: Error while reading file \"%s\".\n",filename);
         regionsFile.close();
         return 1;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR315: File \"%s\" is not available for reading.\n",filename);
      return 1;
   }

   regionsFile.close();

   return 0;
}

// last ERROR - shared with fieldPoints.cpp and frequencyPlan.cpp










