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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "mesh.hpp"
#include <limits>
#include <cfloat>

using namespace std;

#ifndef MATERIALS_H
#define MATERIALS_H

void get_token_pair (string *, string *, string *, int *, string);
bool is_double (string *);
bool is_int (string *);
bool double_compare (double, double, double);

struct point {
   double x;
   double y;
};

class keywordPair
{
   private:
      vector<string> aliases;  // token name plus aliases ex. "frequency", "freq", "f"
      string keyword;          // the text value from the input file
      string value;            // the text value from the input file
      int lineNumber;

      int int_value;
      double dbl_value;
      bool bool_value;
      struct point point_value;

      bool loaded;

      double lowerLimit;
      double upperLimit;
      bool positive_required;
      bool non_negative_required;

      string indent="   ";  // for error messages
      double dbl_tolerance=1e-14;
      bool checkLimits=true;
   public:
      void push_alias (string a) {aliases.push_back(a);}
      bool match_alias (string *);
      void set_keyword (string a) {keyword=a;}
      void set_value (string a) {value=a;}
      void set_int_value (int i) {int_value=i;}
      void set_point_value (double x, double y) {point_value.x=x; point_value.y=y;}
      void set_dbl_value (double a) {dbl_value=a;}
      void set_bool_value (bool a) {bool_value=a;}
      void set_lineNumber (int a) {lineNumber=a;}
      void set_loaded (bool a) {loaded=a;}
      void set_positive_required (bool a) {positive_required=a;}
      void set_non_negative_required (bool a) {non_negative_required=a;}
      void set_lowerLimit (double a) {lowerLimit=a;}
      void set_upperLimit (double a) {upperLimit=a;}
      void set_checkLimits (bool a) {checkLimits=a;}

      bool is_loaded () {return loaded;}
      string get_keyword () {return keyword;} 
      string get_value () {return value;}
      int get_lineNumber () {return lineNumber;}
      int get_int_value() {return int_value;}
      struct point get_point_value() {return point_value;}
      double get_point_value_x() {return point_value.x;}
      double get_point_value_y() {return point_value.y;}
      double get_point_distance (keywordPair *);
      bool is_close_point (keywordPair *);
      double get_dbl_value () {return dbl_value;}
      bool get_bool_value () {return bool_value;}
      bool dbl_compare(keywordPair *a);
      bool value_compare(keywordPair *a);
      bool point_compare(keywordPair *a);

      bool loadBool (string *, string *, int);
      bool loadInt (string *, string *, int);
      bool loadDouble (string *, string *, int);
      bool loadPoint (string *, string *, int);
      bool int_limit_checks(string *, int);
      bool dbl_limit_checks(string *, int);
      bool point_limit_checks(string *, int);
      bool limit_check (string);

      void copy (keywordPair a);
      keywordPair* clone ();

      bool is_any() {
         if (value.compare("any") == 0) return true;
         return false;
      }
      void print();
};


class inputFile
{
   private:
      vector<string> lineTextList;
      vector<int> lineNumberList;
      vector<int> crossReferenceList;   // for easy lookup
      string indent="   ";
   public:
      bool load(const char *);
      void createCrossReference();
      bool checkVersion(string, string);
      void print(); 
      bool findBlock(int, int, int *, int *, string, string, bool);
      unsigned long int get_size() {return lineNumberList.size();}
      int get_lineNumber(int i) {return lineNumberList[i];}
      int get_first_lineNumber();
      int get_last_lineNumber();
      int get_previous_lineNumber (int);
      int get_next_lineNumber (int);
      string get_line (int);
};

class Frequency
{
   private:
      int startLine;                     // inclusive of "Frequency"
      int endLine;                       // inclusive of "EndFrequency"
      keywordPair frequency;
      keywordPair relative_permittivity;
      keywordPair relative_permeability;
      keywordPair loss;                  // loss tangent or conductivity
      keywordPair Rz;                    // surface roughness
   public:
      Frequency (int,int,bool);
      Frequency (){}
      bool load (string *, inputFile *);
      bool inFrequencyBlock(int);
      keywordPair* get_frequency() {return &frequency;}
      keywordPair* get_relative_permittivity() {return &relative_permittivity;}
      keywordPair* get_relative_permeability() {return &relative_permeability;}
      keywordPair* get_loss() {return &loss;}
      keywordPair* get_Rz() {return &Rz;}
      int get_startLine() {return startLine;}
      void print(string);
      bool check(string);
};


class Temperature
{
   private:
      int startLine;                     // inclusive of "Temperature"
      int endLine;                       // inclusive of "EndTemperature"
      vector<Frequency *> frequencyList;
      keywordPair temperature;
      keywordPair er_infinity;           // for Debye model - no frequency blocks with Debye and vice versa
      keywordPair delta_er;
      keywordPair m1;
      keywordPair m2;
      keywordPair relative_permeability;
      keywordPair loss;                  // loss tangent or conductivity
   public:
      Temperature (int,int,bool);
      Temperature (){}
      ~Temperature();
      bool findFrequencyBlocks(inputFile *, bool);
      bool inFrequencyBlocks(int);
      bool inTemperatureBlock (int);
      keywordPair* get_temperature() {return &temperature;}
      Frequency* get_frequency(int i) {return frequencyList[i];}
      int get_startLine() {return startLine;}
      complex<double> get_eps (double, double, string);
      double get_mu (double, double, string);
      double get_Rs (double, double, string);
      bool load(string *, inputFile *, bool);
      void print(string);
      bool check(string);
};

class Source
{
   private:
      int startLine;  // inclusive of "Source"
      int endLine;    // inclusive of "EndSource"
      vector<int> lineNumberList;
      vector<string> lineList;
   public:
      Source (int,int);  // startLine,endLine
      bool inSourceBlock (int);
      bool load(inputFile *);
      void print(string);
};


class Material
{
   private:
      int startLine;  // inclusive of "Material"
      int endLine;    // inclusive of "EndMaterial"
      vector<Temperature *> temperatureList;
      vector<Source *> sourceList;
      keywordPair name;
      bool merged=false;
   public:
      Material (int,int);  // startLine,endLine
      Material (){}
      ~Material();
      bool load (string *, inputFile *, bool);
      bool get_merged () {return merged;}
      void set_merged (bool a) {merged=a;}
      bool findTemperatureBlocks(inputFile *, bool);
      bool findSourceBlocks(inputFile *);
      bool inTemperatureBlocks(int);
      bool inSourceBlocks(int);
      keywordPair* get_name() {return &name;}
      int get_startLine() {return startLine;}
      Temperature* get_temperature(double, double, string);
      complex<double> get_eps(double, double, double, string);
      double get_mu(double, double, double, string);
      double get_Rs(double, double, double, string);
      void print(string);
      bool check(string);
};

class MaterialDatabase
{
   private:
      inputFile inputs;
      vector<Material *> materialList;
      double tol=1e-12;     // tolerance for floating point matches
      string indent="   ";  // for error messages
      string version_name="#OpenParEMmaterials";
      string version_value="1.0";
      double isTransferred=false;
   public:
      ~MaterialDatabase();
      bool load (const char *, const char *, bool);
      bool merge(MaterialDatabase *, string);
      void push(Material *a) {materialList.push_back(a);}
      void print(string);
      bool findMaterialBlocks();
      bool check();
      Material* get(string);
      double get_tol() {return tol;}
      string get_indent() {return indent;}
};


#endif

