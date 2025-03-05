////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    builder - A geometry builder for OpenParEM2D                            //
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

#ifndef BUILDER_H
#define BUILDER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include "petscsys.h"
#include <vector>
#include "keywordPair.hpp"

using namespace std;

class RectangularWaveguide;
class Strip;
class CoupledStrip;

class Point
{
   private:
      double x;
      double y;
      double z;
      double grid;
      int index;
      bool active;
   public:
      Point (int index_, double x_, double y_, double grid_, bool active_) {
         x=x_;
         y=y_;
         z=0;
         index=index_;
         grid=grid_;
         active=active_;
      }
      Point (int index_, double x_, double y_, double z_, double grid_, bool active_) {
         x=x_;
         y=y_;
         z=z_;
         index=index_;
         grid=grid_;
         active=active_;
      }
      double get_x() {return x;}
      double get_y() {return y;}
      double get_z() {return z;}
      double get_grid() {return grid;}
      int get_index() {return index;}
      bool is_active() {return active;}
      void set_x(double x_) {x=x_;}
      void set_y(double y_) {y=y_;}
      void set_z(double z_) {z=z_;}
      void set_grid(double grid_) {grid=grid_;}
      void set_index(int index_) {index=index_;}
      void set_active() {active=true;}
      void set_inactive() {active=false;}
      void shift (double shift) {x+=shift;}
      void offset_index (int n) {index+=n;}
      bool is_equal (Point *);
      double distance (Point *a) {return sqrt(pow(x-a->x,2)+pow(y-a->y,2)+pow(z-a->z,2));}
      void write(ofstream *);
      Point* clone ();
      void print();
};

class Line
{
   private:
      Point *first;
      Point *second;
      int index;
      bool active;
   public:
      Line (int index_, Point *first_, Point *second_, bool active_) {
         first=first_;
         second=second_;
         index=index_;
         active=active_;
      }
      Point* get_first () {return first;}
      Point* get_second () {return second;}
      bool is_active () {return active;}
      void set_first (Point *first_) {first=first_;}
      void set_second (Point *second_) {second=second_;}
      void set_index (int index_) {index=index_;} 
      void set_active () {active=true;}
      void set_inactive () {active=false;}
      int get_index () {return index;}
      bool get_active () {return active;}
      void offset_index (int n) {index+=n;}
      bool is_overlap (Line *);
      bool is_touching (Line *);
      bool is_touching_first (Point *);
      bool is_touching_second (Point *);
      bool reverse (Line *, bool);
      void write (ofstream *);
};

class Curve
{
   private:
      string name;
      vector<Line *> lineList;
      vector<bool> reverseList;
      int index;
      bool active;
   public:
      Curve (int index_, string name_, bool active_) {
         name=name_;
         index=index_;
         active=active_;
      }
      void push_line (Line *line, bool reverse) {
         lineList.push_back(line);
         reverseList.push_back(reverse);
      }
      string get_name() {return name;}
      int get_index() {return index;}
      bool is_active() {return active;}
      void set_active() {active=true;}
      void set_inactive() {active=false;}
      void set_index (int index_) {index=index_;}
      void offset_index (int n) {index+=n;}
      long unsigned int get_lineList_size() {return lineList.size();}
      Line* get_line (int i) {return lineList[i];}
      bool get_reverse (int i) {return reverseList[i];}
      bool is_touching (Line *);
      void merge (Curve *);
      void write (ofstream *);
      void print ();
};

class Surface
{
   private:
      string name;
      vector<Curve *> curveList;
      int index;
      bool active;
   public:
      Surface (int index_, string name_, bool active_) {
         name=name_;
         index=index_;
         active=active_;
      }
      void push_curve (Curve *curve) {
         curveList.push_back(curve);
      }
      string get_name() {return name;}
      int get_index() {return index;}
      bool is_active() {return active;}
      void set_active() {active=true;}
      void set_inactive() {active=false;}
      //void offset_index (int n) {index+=n;}
      //long unsigned int get_curveList_size() {return curveList.size();}
      //Curve* get_curve (int i) {return curveList[i];}
      void write (ofstream *);
};

class Geo
{
   private:
      vector<Point *> pointList;
      vector<Line *> lineList;
      vector<Curve *> curveList;
      vector<Surface *> surfaceList;
   public:
      void push_point (Point *point) {pointList.push_back(point);}
      void push_line (Line *line) {lineList.push_back(line);}
      void push_curve (Curve *curve) {curveList.push_back(curve);}
      void push_surface (Surface *surface) {surfaceList.push_back(surface);}
      long unsigned int get_pointList_size () {return pointList.size();}
      long unsigned int get_lineList_size () {return lineList.size();}
      long unsigned int get_curveList_size () {return curveList.size();}
      long unsigned int get_surfaceList_size () {return surfaceList.size();}
      Point* get_point (int i) {return pointList[i];}
      Line* get_line (int i) {return lineList[i];}
      Curve* get_curve (int i) {return curveList[i];}
      Surface* get_surface (int i) {return surfaceList[i];}
      Point* find_point (int);
      Line* find_line (int);
      void shift (double);
      void merge (Geo *);
      void offset_indices (int);
      long unsigned int find_point (Point *);
      long unsigned int find_line (Line *);
      long unsigned int find_curve (Curve *);
      int get_next_point_index ();
      int get_next_line_index ();
      int get_next_curve_index ();
      int get_next_surface_index ();
      void extrude (double);
      void write (ofstream *);
      void audit ();
};


class Control
{
   private:
      int startLine;                     // inclusive of "Control"
      int endLine;                       // inclusive of "EndControl"
      keywordPair build;
      keywordPair checkLimits;
   public:
      Control (int,int,bool);
      Control(){}
      bool load (string *, inputFile *, bool);
      bool inBlock(int);
      keywordPair get_build() {return build;}
      keywordPair get_checkLimits() {return checkLimits;}
      int get_startLine() {return startLine;}
      void print(string);
      void print_commented(ofstream *);
      bool check(string);
};


// included_type: 0 for RectangularWaveguide
//                1 for Strip
//                2 for CoupledStrip

class RectangularWaveguide
{
   private:
      int startLine;                     // inclusive of "RectangularWaveguide"
      int endLine;                       // inclusive of "EndRectangularWaveguide"
      keywordPair name;
      keywordPair include;
      keywordPair width;
      keywordPair height;
      keywordPair material;
      keywordPair default_conductor_material;
      keywordPair conductor_material_top;
      keywordPair conductor_material_bottom;
      keywordPair conductor_material_left;
      keywordPair conductor_material_right;
      keywordPair length;
      int included_type;
      void *included; 
   public:
      RectangularWaveguide (int,int,bool);
      RectangularWaveguide(){}
      bool load (string *, inputFile *, bool);
      bool inBlock(int);
      keywordPair get_name() {return name;}
      void set_included_type(int included_type_) {included_type=included_type_;}
      int get_included_type() {return included_type;}
      keywordPair get_include() {return include;}
      keywordPair get_width() {return width;}
      keywordPair get_height() {return height;}
      keywordPair get_material() {return material;}
      keywordPair get_default_conductor_material() {return default_conductor_material;}
      keywordPair get_conductor_material_top () {if (conductor_material_top.is_loaded()) {return conductor_material_top;} else {return default_conductor_material;}}
      keywordPair get_conductor_material_bottom () {if (conductor_material_bottom.is_loaded()) {return conductor_material_bottom;} else {return default_conductor_material;}}
      keywordPair get_conductor_material_left () {if (conductor_material_left.is_loaded()) {return conductor_material_left;} else {return default_conductor_material;}}
      keywordPair get_conductor_material_right () {if (conductor_material_right.is_loaded()) {return conductor_material_right;} else {return default_conductor_material;}}
      void set_included (void *included_) {included=included_;}
      void* get_included() {return included;}
      int get_startLine() {return startLine;}
      void apply_include ();
      void print(string);
      void print_commented(ofstream *);
      bool check(string);
      bool checkInclude(string);
      bool checkLimits();
      void build_geo (Geo *);
      bool write_geo(string, Control *, Geo *);
      bool write_modes_and_boundaries (string, Control *);
      bool write_ports_and_boundaries (string, Control *);
      bool write_OpenParEM2D_proj(string);
      bool write_OpenParEM3D_proj(string);
};

class Strip
{
   private:
      int startLine;                          // inclusive of "Strip"
      int endLine;                            // inclusive of "EndStrip"
      keywordPair name;
      keywordPair include;
      keywordPair use_symmetry;
      keywordPair default_conductor_material;
      keywordPair left_side_gap;              // measured from trace_width/2
      keywordPair left_side_material;
      keywordPair right_side_gap;             // measured from trace_width/2
      keywordPair right_side_material;
      keywordPair upper_thickness;            // measured from the top of the lower thickness
      keywordPair upper_material;
      keywordPair soldermask_thickness;
      keywordPair soldermask_material;
      keywordPair lower_thickness;
      keywordPair lower_material;
      keywordPair trace_thickness;
      keywordPair trace_width;
      keywordPair trace_etch_angle;
      keywordPair trace_material_top;
      keywordPair trace_material_bottom;
      keywordPair trace_material_sides;
      keywordPair upper_groundplane_material;
      keywordPair lower_groundplane_material;
      keywordPair length;
      int included_type;
      void *included;
   public:
      Strip (int,int,bool);
      Strip(){}
      bool load (string *, inputFile *, bool);
      bool inBlock(int);
      keywordPair get_name() {return name;}
      int get_included_type() {return included_type;}
      void set_included_type(int included_type_) {included_type=included_type_;}
      keywordPair get_include() {return include;}
      keywordPair get_use_symmetry() {return use_symmetry;}
      keywordPair get_upper_material() {return upper_material;}
      keywordPair get_upper_thickness() {return upper_thickness;}
      keywordPair get_soldermask_thickness() {return soldermask_thickness;}
      keywordPair get_soldermask_material() {return soldermask_material;}
      keywordPair get_lower_thickness() {return lower_thickness;}
      keywordPair get_lower_material() {return lower_material;}
      keywordPair get_trace_thickness() {return trace_thickness;}
      keywordPair get_trace_width() {return trace_width;}
      keywordPair get_trace_etch_angle() {return trace_etch_angle;}
      keywordPair get_default_conductor_material() {return default_conductor_material;}
      keywordPair get_trace_material_top () {if (trace_material_top.is_loaded()) {return trace_material_top;} else {return default_conductor_material;}}
      keywordPair get_trace_material_bottom () {if (trace_material_bottom.is_loaded()) {return trace_material_bottom;} else {return default_conductor_material;}}
      keywordPair get_trace_material_sides () {if (trace_material_sides.is_loaded()) {return trace_material_sides;} else {return default_conductor_material;}}
      keywordPair get_upper_groundplane_material () {if (upper_groundplane_material.is_loaded()) {return upper_groundplane_material;} else {return default_conductor_material;}}
      keywordPair get_lower_groundplane_material () {if (lower_groundplane_material.is_loaded()) {return lower_groundplane_material;} else {return default_conductor_material;}}
      keywordPair get_left_side_material () {if (left_side_material.is_loaded()) {return left_side_material;} else {return default_conductor_material;}}
      keywordPair get_left_side_gap () {return left_side_gap;}
      keywordPair get_right_side_material () {if (right_side_material.is_loaded()) {return right_side_material;} else {return default_conductor_material;}}
      keywordPair get_right_side_gap () {return right_side_gap;}
      void set_use_symmetry (bool use_symmetry_) {use_symmetry.set_bool_value(use_symmetry_);}
      void set_default_conductor_material (keywordPair default_conductor_material_) {default_conductor_material.copy(default_conductor_material_);}
      void set_left_side_gap (keywordPair left_side_gap_) {left_side_gap.copy(left_side_gap_);}
      void set_left_side_material (keywordPair left_side_material_) {left_side_material.copy(left_side_material_);}
      void set_right_side_gap (keywordPair right_side_gap_) {right_side_gap.copy(right_side_gap_);}
      void set_right_side_material (keywordPair right_side_material_) {right_side_material.copy(right_side_material_);}
      void set_upper_thickness (keywordPair upper_thickness_) {upper_thickness.copy(upper_thickness_);}
      void set_upper_material (keywordPair upper_material_) {upper_material.copy(upper_material_);}
      void set_soldermask_thickness (keywordPair soldermask_thickness_) {soldermask_thickness.copy(soldermask_thickness_);}
      void set_soldermask_material (keywordPair soldermask_material_) {soldermask_material.copy(soldermask_material_);}
      void set_lower_thickness (keywordPair lower_thickness_) {lower_thickness.copy(lower_thickness_);}
      void set_lower_material (keywordPair lower_material_) {lower_material.copy(lower_material_);}
      void set_trace_thickness (keywordPair trace_thickness_) {trace_thickness.copy(trace_thickness_);}
      void set_trace_width (keywordPair trace_width_) {trace_width.copy(trace_width_);}
      void set_trace_etch_angle (keywordPair trace_etch_angle_) {trace_etch_angle.copy(trace_etch_angle_);}
      void set_trace_material_top (keywordPair trace_material_top_) {trace_material_top.copy(trace_material_top_);}
      void set_trace_material_bottom (keywordPair trace_material_bottom_) {trace_material_bottom.copy(trace_material_bottom_);}
      void set_trace_material_sides (keywordPair trace_material_sides_) {trace_material_sides.copy(trace_material_sides_);}
      void set_upper_groundplane_material (keywordPair upper_groundplane_material_) {upper_groundplane_material.copy(upper_groundplane_material_);}
      void set_lower_groundplane_material (keywordPair lower_groundplane_material_) {lower_groundplane_material.copy(lower_groundplane_material_);}
      void set_included (void *included_) {included=included_;}
      void* get_included() {return included;}
      int get_startLine() {return startLine;}
      void apply_include ();
      void print(string);
      void print_commented(ofstream *);
      bool check(string);
      bool checkInclude(string);
      bool checkLimits();
      bool check_geo(Geo *);
      void split_right_side_gap() {right_side_gap.set_dbl_value(right_side_gap.get_dbl_value()/2);}
      void split_left_side_gap() {left_side_gap.set_dbl_value(left_side_gap.get_dbl_value()/2);}
      void build_geo (Geo *, bool);
      bool write_geo(string, Control *, Geo*);
      bool write_modes_and_boundaries (string, Control *);
      bool write_ports_and_boundaries (string, Control *);
      bool write_OpenParEM2D_proj(string);
      bool write_OpenParEM3D_proj(string);
};

class CoupledStrip
{
   private:
      int startLine;                          // inclusive of "CoupledStrip"
      int endLine;                            // inclusive of "EndEndCoupledStrip"
      keywordPair name;
      keywordPair include;
      keywordPair solution_impedance_calculation;
      keywordPair default_conductor_material;
      keywordPair left_side_gap;              // measured from trace_width/2
      keywordPair left_side_material;
      keywordPair right_side_gap;             // measured from trace_width/2
      keywordPair right_side_material;
      keywordPair upper_thickness;            // measured from the top of the lower thickness
      keywordPair upper_material;
      keywordPair soldermask_thickness;
      keywordPair soldermask_material;
      keywordPair lower_thickness;
      keywordPair lower_material;
      keywordPair trace_left_width;
      keywordPair trace_right_width;
      keywordPair trace_thickness;
      keywordPair trace_air_gap;
      keywordPair trace_etch_angle;
      keywordPair trace_material_top;
      keywordPair trace_material_bottom;
      keywordPair trace_material_sides;
      keywordPair upper_groundplane_material;
      keywordPair lower_groundplane_material;
      keywordPair length;
      int included_type;
      void *included;
   public:
      CoupledStrip (int,int,bool);
      CoupledStrip(){}
      bool load (string *, inputFile *, bool);
      bool inBlock(int);
      keywordPair get_name() {return name;}
      int get_included_type() {return included_type;}
      void set_included_type(int included_type_) {included_type=included_type_;}
      keywordPair get_include() {return include;}
      keywordPair get_upper_material() {return upper_material;}
      keywordPair get_upper_thickness() {return upper_thickness;}
      keywordPair get_soldermask_thickness() {return soldermask_thickness;}
      keywordPair get_soldermask_material() {return soldermask_material;}
      keywordPair get_lower_thickness() {return lower_thickness;}
      keywordPair get_lower_material() {return lower_material;}
      keywordPair get_trace_left_width() {return trace_left_width;}
      keywordPair get_trace_right_width() {return trace_right_width;}
      keywordPair get_trace_thickness() {return trace_thickness;}
      keywordPair get_trace_air_gap() {return trace_air_gap;}
      keywordPair get_trace_etch_angle() {return trace_etch_angle;}
      keywordPair get_default_conductor_material() {return default_conductor_material;}
      keywordPair get_trace_material_top () {if (trace_material_top.is_loaded()) {return trace_material_top;} else {return default_conductor_material;}}
      keywordPair get_trace_material_bottom () {if (trace_material_bottom.is_loaded()) {return trace_material_bottom;} else {return default_conductor_material;}}
      keywordPair get_trace_material_sides () {if (trace_material_sides.is_loaded()) {return trace_material_sides;} else {return default_conductor_material;}}
      keywordPair get_upper_groundplane_material () {if (upper_groundplane_material.is_loaded()) {return upper_groundplane_material;} else {return default_conductor_material;}}
      keywordPair get_lower_groundplane_material () {if (lower_groundplane_material.is_loaded()) {return lower_groundplane_material;} else {return default_conductor_material;}}
      keywordPair get_left_side_material () {if (left_side_material.is_loaded()) {return left_side_material;} else {return default_conductor_material;}}
      keywordPair get_left_side_gap () {return left_side_gap;}
      keywordPair get_right_side_material () {if (right_side_material.is_loaded()) {return right_side_material;} else {return default_conductor_material;}}
      keywordPair get_right_side_gap () {return right_side_gap;}
      void set_included (void *included_) {included=included_;}
      void* get_included() {return included;}
      int get_startLine() {return startLine;}
      void apply_include ();
      void print(string);
      void print_commented(ofstream *);
      bool check(string);
      bool checkInclude(string);
      bool checkLimits();
      void build_geo (Geo *);
      bool write_geo(string, Control *, Geo*);
      bool write_modes_and_boundaries (string, Control *);
      bool write_ports_and_boundaries (string, Control *);
      bool write_OpenParEM2D_proj (string);
      bool write_OpenParEM3D_proj (string);
};

class StructureDatabase
{
   private:
      inputFile inputs;
      vector<Point *> pointList;
      vector<Line *> lineList;
      vector<Control *> controlList;
      vector<RectangularWaveguide *> rectangularWaveguideList;
      vector<Strip *> stripList;
      vector<CoupledStrip *> coupledStripList;
      double tol=1e-12;     // tolerance for floating point matches
      string indent="";  // for error messages
      string version_name="#builder";
      string version_value="1.0";
      Geo geo;
   public:
      ~StructureDatabase();
      bool load (const char *, const char *, bool);
      bool set_include();
      void apply_includes();
      void print();
      bool findBlocks(bool);
      bool inBlocks(int);
      bool check();
      bool checkInclude();
      bool checkLimits();
      RectangularWaveguide* get_rectangularWaveguide(string);
      double get_tol() {return tol;}
      string get_indent() {return indent;}
      bool build_geo();
      bool write_geo();
      bool write_modes_and_boundaries();
      bool write_proj();
      void audit () {geo.audit();}
};

#endif
