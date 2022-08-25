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

#ifndef MODES_H
#define MODES_H

#include "mfem.hpp"
#include "mesh.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include "OpenParEMmaterials.hpp"

using namespace std;
using namespace mfem;

class SourceFile2 {
   private:
      int startLine;
      int endLine;
      keywordPair name;
   public:
      SourceFile2(int,int);
      bool load(string *, inputFile *);
      bool inBlock(int);
      string get_name() {return name.get_value();}
      bool check(string *);
      void print();
};

class Path {
   private:
      int startLine;
      int endLine;
      keywordPair name;
      vector<keywordPair *> points;
      keywordPair closed;
      double tol=1e-12;
   public:
      Path (int, int);
      ~Path();
      bool load(string *, inputFile *);
      bool inBlock (int);
      bool check(string *);
      bool checkBoundingBox(Vector *, Vector *, string *, double);
      string get_name() {return name.get_value();}
      bool name_is_loaded() {return name.is_loaded();}
      int get_name_lineNumber() {return name.get_lineNumber();}
      bool get_closed() {return closed.get_bool_value();}
      int get_startLine() {return startLine;}
      int get_endLine() {return endLine;}
      long unsigned int get_points_size() {return points.size();}
      point get_point (long unsigned int i) {return points[i]->get_point_value();}
      double get_point_x (long unsigned int i) {return points[i]->get_point_value_x();}
      double get_point_y (long unsigned int i) {return points[i]->get_point_value_y();}
      bool compare (long unsigned int i, keywordPair test_point);
      bool is_closed() {return closed.get_bool_value();}
      keywordPair* get_startPoint() {return points[0];}
      keywordPair* get_endPoint() {if (closed.get_bool_value()) return points[0]; return points[points.size()-1];}
      long unsigned int is_segmentOnLine (double, double, double, double);
      void subdivide (Path *);
      void print();
};

// storage for the information stored in for an attribute tagged onto an edge element in the mesh
struct EdgeAttribute {
   long unsigned int boundary;   // index to a boundary
   long unsigned int path;       // index to a path on the boundary
   long unsigned int segment;    // index to a segment of the path
};

class Border
{
   private:
      // An edge element can have only one boundary condition, but it can be part of many modes.
      struct EdgeAttribute boundary;          // for a BoundaryDatabase for a boundary condition (Zs, PEC, PMC)
      vector<struct EdgeAttribute> mode;      // for a BoundaryDatabase for a mode, indexed by modeNumber
   public:
      Border();
      bool has_boundary();
      bool has_mode();
      bool has_mode(long unsigned int);
      void set(struct EdgeAttribute);
      void set(struct EdgeAttribute, long unsigned int);
      bool exists(struct EdgeAttribute);
      bool exists(struct EdgeAttribute, long unsigned int);
      bool exists_any_mode(struct EdgeAttribute);
      bool exists_test_boundary_only(struct EdgeAttribute);
      bool mode_exists(long unsigned int);
      struct EdgeAttribute get_boundary() {return boundary;}
      struct EdgeAttribute get_mode(long unsigned int);
      void print ();
};

class BorderDatabase
{
   private:
      vector<Border *> borderList;
   public:
      BorderDatabase();
      ~BorderDatabase();
      long unsigned int exists (struct EdgeAttribute);
      long unsigned int exists (struct EdgeAttribute, long unsigned int);
      long unsigned int exists_any_mode (struct EdgeAttribute);
      long unsigned int exists_test_boundary_only (struct EdgeAttribute);
      long unsigned int add(struct EdgeAttribute, int);
      long unsigned int add(struct EdgeAttribute, long unsigned int, int);
      long unsigned int get_size() {return borderList.size();}
      Border* get_border(long unsigned int i) {return borderList[i];}
      Array<int>* build_entry_list (bool, long unsigned int, int);
      void print ();
};

class Boundary
{
   private:
      int startLine;
      int endLine;
      keywordPair name;                    // for boundaries only
      keywordPair mode;                    // for modes only (numerical, 0-based)
      keywordPair mode_block_type;         // for modes only (modal, line)
      keywordPair type;                    // boundaries: (surface impedanec, perfect electric conductor, perfect magnetic conductor); modes/lines: (voltage, current)
      keywordPair material;                // surface impedance boundary only
      keywordPair attribute;
      vector<keywordPair *> pathNameList;
      vector<long unsigned int> pathIndexList;
      vector<bool> reverseList;
   public:
      Boundary (int,int,string);
      ~Boundary();
      bool load(string *, inputFile *, bool, int);
      bool inBlock (int);
      bool check(string *);
      bool checkBoundingBox(Vector *, Vector *, string *, double, vector<Path *> *);
      bool check_current_paths (string *, vector<Path *> *, bool);
      int get_startLine() {return startLine;}
      int get_endLine() {return endLine;}
      string get_name() {return name.get_value();}
      int get_mode() {return mode.get_int_value();}
      bool name_is_loaded() {return name.is_loaded();}
      bool mode_is_loaded() {return mode.is_loaded();}
      int get_name_lineNumber() {return name.get_lineNumber();}
      int get_mode_lineNumber() {return mode.get_lineNumber();}
      int get_attribute() {return attribute.get_int_value();}
      string get_type() {return type.get_value();}
      string get_material() {return material.get_value();}
      long unsigned int get_pathNameList_size() {return pathNameList.size();}
      string get_pathName(long unsigned int i) {return pathNameList[i]->get_value();}
      int get_pathName_lineNumber(long unsigned int i) {return pathNameList[i]->get_lineNumber();}
      bool get_reverse(long unsigned int i) {return reverseList[i];}
      long unsigned int get_path_size() {return pathIndexList.size();}
      long unsigned int get_path(long unsigned int i) {return pathIndexList[i];}
      void push(long unsigned int a) {pathIndexList.push_back(a);}
      struct EdgeAttribute is_segmentOnPath (double, double, double, double, vector<Path *> *);
      void set_attribute(int i) {attribute.set_int_value(i);}
      void markMeshBoundaries (Mesh *, BorderDatabase *, vector<Path *> *);
      bool is_surface_impedance();
      bool is_perfect_electric_conductor();
      bool is_perfect_magnetic_conductor();
      bool is_mode_voltage();
      bool is_mode_current();
      bool is_boundary();
      bool is_mode();
      bool is_modal();
      bool is_line();
      string get_block_type();
      string get_mode_name();
      void print();
};

class BoundaryDatabase
{
   private:
      inputFile inputs;
      vector<SourceFile2 *> sourceFileList;
      vector<Path *> pathList;
      vector<Boundary *> boundaryList;
      double tol=1e-14;
      string indent="   ";
      string version_name="#OpenParEMmodes";
      string version_value="1.0";
   public:
      ~BoundaryDatabase();
      bool inBlocks(int);
      bool check();
      bool check_scale (Mesh *, int);
      void subdivide_paths ();
      bool findSourceFileBlocks();
      bool findPathBlocks();
      bool findBoundaryBlocks();
      bool findModeBlocks();
      bool findLineBlocks();
      bool load(const char *, bool);
      void mark_boundaries (Mesh *, BorderDatabase *);
      long unsigned int get_boundary_size() {return boundaryList.size();}
      Boundary* get_boundary (long unsigned int i) {return boundaryList[i];}
      vector<Path *>* get_pathList() {return &pathList;}
      void print();
};

#endif

