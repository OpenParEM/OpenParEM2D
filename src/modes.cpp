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

#include "modes.hpp"

// with common point (x0,y0)
double angle_between_two_lines (double x0, double y0, double x1, double y1, double x2, double y2)
{
   double theta=atan2(y1-y0,x1-x0)-atan2(y2-y0,x2-x0);

   return theta;
}

bool are_parallel (double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
   double tol=1e-12;
   double theta=angle_between_two_lines(0,0,x1-x0,y1-y0,x3-x2,y3-y2);
   if (abs(theta) < tol) return true;
   return false;
}

bool compare_xy (double x1, double y1, double x2, double y2)
{
   double tol=1e-12;
   if (double_compare(x1,x2,tol) && double_compare(y1,y2,tol)) return true;
   return false;
}

// checks to see if test point (xt,yt) falls on the line given by (x1,y1) to (x2,y2)
bool is_point_on_line (double xt, double yt, double x1, double y1, double x2, double y2, double tolerance)
{
   double beta,betat,betar;
   double length,lengtht;
   double xtr,ytr;
   double tol;

   // angle and length of the line
   beta=atan2(y2-y1,x2-x1);
   length=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

   // angle and length of the test point
   betat=atan2(yt-y1,xt-x1);
   lengtht=sqrt((xt-x1)*(xt-x1)+(yt-y1)*(yt-y1));

   // rotate the line and the test point so that the line aligns with the x-axis: rotate by -beta
   betar=betat-beta;

   // rotate the test point
   xtr=lengtht*cos(betar);
   ytr=lengtht*sin(betar);

   // scale the tolerance
   tol=tolerance*length;

   // check if the rotated test point falls on the rotated line
   if (xtr >= -tolerance && xtr <= length+tolerance && fabs(ytr) <= tol) return true;

   return false;
}

bool is_point_on_line_not_ends (double xt, double yt, double x1, double y1, double x2, double y2, double tolerance)
{
   if (compare_xy(xt,yt,x1,y1)) return false;
   if (compare_xy(xt,yt,x2,y2)) return false;
   if (is_point_on_line (xt,yt,x1,y1,x2,y2,tolerance)) return true;
   return false;
}

void test_is_point_on_line ()
{
  int i=1;
  double x1,y1,x2,y2,tolerance;

  x1=1; y1=1; x2=10,y2=10; tolerance=1e-8/9;
  if (is_point_on_line (2, 2, x1, y1, x2, y2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, 1-1e-9, x1, y1, x2, y2, tolerance))    PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, 10+1e-9, x1, y1, x2, y2, tolerance))  PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 2+1e-10, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++; 
  if (is_point_on_line (2, 2-1e-10, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++; 
  if (!is_point_on_line (-2, -2, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++; 
  if (!is_point_on_line (11, 11, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 3, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++; 
  if (!is_point_on_line (3, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++; 
  if (!is_point_on_line (2+1e-7, 2, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2-1e-7, 2, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=1; y1=-1; x2=10,y2=-10;
  if (is_point_on_line (2, -2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, -1-1e-9, x1, y1, x2, y2, tolerance))   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, -10+1e-9, x1, y1, x2, y2, tolerance)) PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -2+1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -2-1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, -11, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (3, -2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2+1e-7, -2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2-1e-7, -2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=1; x2=-10,y2=10;
  if (is_point_on_line (-2, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1-1e-10, 1-1e-9, x1, y1, x2, y2, tolerance))   PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10+1e-10, 10+1e-9, x1, y1, x2, y2, tolerance)) PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 2+1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 2-1e-10, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, 11, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-3, 2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2+1e-7, 2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2-1e-7, 2, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=-1; x2=-10,y2=-10;
  if (is_point_on_line (-2, -2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1-1e-10, -1-1e-9, x1, y1, x2, y2, tolerance))  PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10+1e-10, -10+1e-9, x1, y1, x2, y2, tolerance))PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -2+1e-10, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -2-1e-10, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, -11, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -3, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-3, -2, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2+1e-7, -2, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2-1e-7, -2, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=-1; y1=0; x2=-10,y2=0;
  if (is_point_on_line (-2, 0, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-1+1e-10, 0, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-10-1e-10, 0, x1, y1, x2, y2, tolerance))       PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, 1e-10, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (-2, -1e-10, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 2, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-11, -11, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, -3, x1, y1, x2, y2, tolerance))            PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (-2, 1e-7, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++; 
  if (!is_point_on_line (-2, -1e7, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;

  x1=1; y1=0; x2=10,y2=0;
  if (is_point_on_line (2, 0, x1, y1, x2, y2, tolerance))               PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (1-1e-10, 0, x1, y1, x2, y2, tolerance))         PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (10+1e-10, 0, x1, y1, x2, y2, tolerance))        PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, 1e-10, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (is_point_on_line (2, -1e-10, x1, y1, x2, y2, tolerance))          PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 2, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (11, 0, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -3, x1, y1, x2, y2, tolerance))             PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 3, x1, y1, x2, y2, tolerance))              PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, 1e-7, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
  if (!is_point_on_line (2, -1e7, x1, y1, x2, y2, tolerance))           PetscPrintf(PETSC_COMM_WORLD,"%d PASS\n",i); else PetscPrintf(PETSC_COMM_WORLD,"%d FAIL\n",i); i++;
}

bool operator==(struct EdgeAttribute a, struct EdgeAttribute b)
{
   if (a.boundary != b.boundary) return false;
   if (a.path != b.path) return false;
   if (a.segment != b.segment) return false;
   return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Path
///////////////////////////////////////////////////////////////////////////////////////////

Path::Path(int startLine_, int endLine_)
{
   startLine=startLine_;
   endLine=endLine_;

   // name
   name.push_alias("name");
   name.set_loaded(false);
   name.set_positive_required(false);
   name.set_non_negative_required(false);
   name.set_lowerLimit(0);
   name.set_upperLimit(0);
   name.set_checkLimits(false);

   // closed
   closed.push_alias("closed");
   closed.set_loaded(false);
   closed.set_positive_required(false);
   closed.set_non_negative_required(false);
   closed.set_lowerLimit(0);
   closed.set_upperLimit(0);
   closed.set_checkLimits(false);
}


// return true if the point is close
bool Path::compare (long unsigned int i, keywordPair test_point)
{
   if (points[i]->get_point_value_x() == 0) {
      if (fabs(test_point.get_point_value_x()) > tol) return false;
   }
   if (fabs((test_point.get_point_value_x()-points[i]->get_point_value_x())/points[i]->get_point_value_x()) > tol) return false;

   if (points[i]->get_point_value_y() == 0) {
      if (fabs(test_point.get_point_value_y()) > tol) return false;
   }
   if (fabs((test_point.get_point_value_y()-points[i]->get_point_value_y())/points[i]->get_point_value_y()) > tol) return false;
   return true;
}

void Path::print ()
{
   PetscPrintf(PETSC_COMM_WORLD,"Path\n");
   PetscPrintf(PETSC_COMM_WORLD,"   name=%s\n",get_name().c_str());
   long unsigned int i=0;
   while (i < points.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"   point=(%g,%g)\n",get_point_x(i),get_point_y(i));
      i++;
   }
   if (closed.get_bool_value()) PetscPrintf(PETSC_COMM_WORLD,"   closed=true\n");
   else PetscPrintf(PETSC_COMM_WORLD,"   closed=false\n");
   PetscPrintf(PETSC_COMM_WORLD,"EndPath\n");
}

bool Path::load(string *indent, inputFile *inputs)
{
   bool fail=false;

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      string token,value,line;
      line=inputs->get_line(lineNumber);
      get_token_pair(&line,&token,&value,&lineNumber,*indent);

      int recognized=0;

      if (name.match_alias(&token)) {
         if (name.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR200: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber,name.get_lineNumber());
            fail=true;
         } else {
            name.set_keyword(token);
            name.set_value(value);
            name.set_lineNumber(lineNumber);
            name.set_loaded(true);
         }
         recognized++;
      }

      if (token.compare("point") == 0) {
         keywordPair *point=new keywordPair;
         point->push_alias("point");
         point->set_keyword(token);
         point->set_value(value);
         point->set_lineNumber(lineNumber);
         point->set_positive_required(false);
         point->set_non_negative_required(false);
         point->set_lowerLimit(-100);
         point->set_upperLimit(100);
         point->set_checkLimits(true);
         point->set_loaded(false);

         if (point->loadPoint(&token,&value,lineNumber)) delete point;
         else points.push_back(point);

         recognized++;
      }

      if (closed.match_alias(&token)) {
         recognized++;
         if (closed.loadBool(&token, &value, lineNumber)) fail=true;
      }

      // should recognize one keyword
      if (recognized != 1) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR201: Unrecognized keyword at line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
         fail=true;
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   return fail;
}

bool Path::checkBoundingBox(Vector *lowerLeft, Vector *upperRight, string *indent, double tol)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < points.size()) {
      bool point_fail=false;

      if (points[i]->get_point_value_x() < lowerLeft->Elem(0)-tol) point_fail=true;
      if (points[i]->get_point_value_x() > upperRight->Elem(0)+tol) point_fail=true;

      if (points[i]->get_point_value_y() < lowerLeft->Elem(1)-tol) point_fail=true;
      if (points[i]->get_point_value_y() > upperRight->Elem(1)+tol) point_fail=true;

      if (point_fail) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR202: Path block at line %d has point (%g,%g) outside of the mesh bounding box.\n",
                                      indent->c_str(),indent->c_str(),startLine,points[i]->get_point_value_x(),points[i]->get_point_value_y());
         fail=true;
      }

      i++;
   }

   return fail;
}

bool Path::inBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Path::check(string *indent)
{
   bool fail=false;

   if (!name.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR203: Path block at line %d must specify a name.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (!closed.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR204: Path block at line %d must specify \"closed\".\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (points.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR205: Path block at line %d must specify points.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (points.size() == 1) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR206: Path block at line %d must specify more than one point.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   if (points.size() == 2 && closed.get_bool_value()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR207: Path block at line %d cannot be closed with just two points.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   return fail;
}

// returns the segment index on which the segment falls
long unsigned int Path::is_segmentOnLine (double x1, double y1, double x2, double y2)
{
   long unsigned int max=-1;

   long unsigned int i=0;
   while (points.size() > 0 && i < points.size()-1) {
      if (is_point_on_line (x1,y1,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                  points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),1e-8) &&
          is_point_on_line (x2,y2,points[i]->get_point_value_x(),points[i]->get_point_value_y(),
                                  points[i+1]->get_point_value_x(),points[i+1]->get_point_value_y(),1e-8)) return i;
      i++;
   }

   if (is_closed()) {
      if (is_point_on_line (x1,y1,points[points.size()-1]->get_point_value_x(),points[points.size()-1]->get_point_value_y(),
                                  points[0]->get_point_value_x(),points[0]->get_point_value_y(),1e-8) &&
          is_point_on_line (x2,y2,points[points.size()-1]->get_point_value_x(),points[points.size()-1]->get_point_value_y(),
                                  points[0]->get_point_value_x(),points[0]->get_point_value_y(),1e-8)) return points.size()-1;
   }

   return max;
}

// eliminate partial overlaps of paths by subdividing
// crossing paths are ok
void Path::subdivide(Path *test)
{
   bool modified=false;
   vector<keywordPair *> newPoints;

   if (points.size() == 0) return;

   newPoints.push_back(test->points[0]->clone());

   long unsigned int i=0;
   while (i < points.size()-1) {

      double x1=points[i]->get_point_value_x();
      double y1=points[i]->get_point_value_y();
      double x2=points[i+1]->get_point_value_x();
      double y2=points[i+1]->get_point_value_y();

      long unsigned int j=0;
      while (j < test->points.size()-1) {

         double xt1=test->points[j]->get_point_value_x();
         double yt1=test->points[j]->get_point_value_y();
         double xt2=test->points[j+1]->get_point_value_x();
         double yt2=test->points[j+1]->get_point_value_y();

         bool break_on_1=false;
         bool break_on_2=false;

         bool parallel=are_parallel (xt1,yt1,xt2,yt2,x1,y1,x2,y2);

         if (parallel && is_point_on_line_not_ends(xt1,yt1,x1,y1,x2,y2,1e-8)) break_on_1=true;
         if (parallel && is_point_on_line_not_ends(xt2,yt2,x1,y1,x2,y2,1e-8)) break_on_2=true;

         if (break_on_1) {
            if (break_on_2) {
               // segment is fully enclosed

               // maintain ordering along the line
               if (points[i]->get_point_distance(test->points[j]) < points[i]->get_point_distance(test->points[j+1])) {
                  if (! points[points.size()-1]->is_close_point (test->points[j])) {
                     newPoints.push_back(test->points[j]->clone());
                  }
                  newPoints.push_back(test->points[j+1]->clone());
               } else {
                  if (! points[points.size()-1]->is_close_point (test->points[j+1])) {
                     newPoints.push_back(test->points[j+1]->clone());
                  }
                  newPoints.push_back(test->points[j]->clone());
               }
               modified=true;
            } else {
               // partial overlap - break the segment at test point 1
               newPoints.push_back(test->points[j]->clone());
               modified=true;
            }
         } else {
            if (break_on_2) {
               // partial overlap - break the segment at test point 2
               newPoints.push_back(test->points[j+1]->clone());
               modified=true;
            } else {
               // nothing to do
            }
         }

         j++;
      }

      // finish the segment
      newPoints.push_back(points[i+1]->clone());

      i++;
   }

   if (modified) {
      long unsigned int k=0;
      while (k < points.size()) {
         delete points[k];
         k++;
      }
      points.clear();

      k=0;
      while (k < newPoints.size()) {
         points.push_back(newPoints[k]);
         k++;
      }
   } else {
      long unsigned int k=0;
      while (k < newPoints.size()) {
         delete newPoints[k];
         k++;
      }
   }
}

Path::~Path ()
{
   long unsigned int i=0;
   while (i < points.size()) {
      delete points[i];
      i++;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Boundary
///////////////////////////////////////////////////////////////////////////////////////////

Boundary::Boundary(int startLine_, int endLine_, string mode_block_type_)
{
   startLine=startLine_;
   endLine=endLine_;

   // name
   name.push_alias("name");
   name.set_loaded(false);
   name.set_positive_required(false);
   name.set_non_negative_required(false);
   name.set_lowerLimit(0);
   name.set_upperLimit(0);
   name.set_checkLimits(false);

   // mode/line number (conductor number)
   mode.push_alias("mode");
   mode.push_alias("line");
   mode.set_loaded(false);
   mode.set_positive_required(true);
   mode.set_non_negative_required(false);
   mode.set_lowerLimit(1);
   mode.set_upperLimit(10000);
   mode.set_checkLimits(true);

   // mode_block_type
   mode_block_type.push_alias("mode_block_type");
   mode_block_type.set_value(mode_block_type_);
   mode_block_type.set_loaded(true);
   mode_block_type.set_positive_required(false);
   mode_block_type.set_non_negative_required(false);
   mode_block_type.set_lowerLimit(0);
   mode_block_type.set_upperLimit(0);
   mode_block_type.set_checkLimits(false);

   // type
   type.push_alias("type");
   type.set_loaded(false);
   type.set_positive_required(false);
   type.set_non_negative_required(false);
   type.set_lowerLimit(0);
   type.set_upperLimit(0);
   type.set_checkLimits(false);

   // material
   material.push_alias("material");
   material.set_loaded(false);
   material.set_positive_required(false);
   material.set_non_negative_required(false);
   material.set_lowerLimit(0);
   material.set_upperLimit(0);
   material.set_checkLimits(false);

   // attribute
   attribute.push_alias("attribute");
   attribute.set_loaded(false);
   attribute.set_positive_required(false);
   attribute.set_non_negative_required(true);
   attribute.set_lowerLimit(0);
   attribute.set_upperLimit(1000000);
   attribute.set_checkLimits(true);
}

// is_boundary=false => load as mode
// is_boundary=true => load as boundary
bool Boundary::load(string *indent, inputFile *inputs, bool is_boundary, int attribute_)
{
   bool fail=false;
   bool found_first_path=false;

   attribute.set_int_value(attribute_);

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      string token,value,line;
      line=inputs->get_line(lineNumber);
      get_token_pair(&line,&token,&value,&lineNumber,*indent);

      int recognized=0;

      if (is_boundary && name.match_alias(&token)) {
         if (name.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR208: Duplicate entry at line %d for previous entry at line %d.\n",indent->c_str(),indent->c_str(),lineNumber,name.get_lineNumber());
            fail=true;
         } else {
            name.set_keyword(token);
            name.set_value(value);
            name.set_lineNumber(lineNumber);
            name.set_loaded(true);
         }
         recognized++;
      }

      if (!is_boundary && mode.match_alias(&token)) {
         if (token.compare("mode") == 0 && is_line()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR249: \"mode\" defined at line %d should be \"line\".\n",indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         if (token.compare("line") == 0 && is_modal()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR250: \"line\" defined at line %d should be \"mode\".\n",indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         if (!fail && mode.loadInt(&token, &value, lineNumber)) fail=true;
         recognized++;
      }

      if (type.match_alias(&token)) {
         if (type.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR209: Duplicate entry at line %d for previous entry at line %d.\n",indent->c_str(),indent->c_str(),lineNumber,type.get_lineNumber());
            fail=true;
         } else {
            type.set_keyword(token);
            type.set_value(value);
            type.set_lineNumber(lineNumber);
            type.set_loaded(true);
         }
         recognized++;
      }

      if (is_boundary && material.match_alias(&token)) {
         if (material.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR210: Duplicate entry at line %d for previous entry at line %d.\n",indent->c_str(),indent->c_str(),lineNumber,material.get_lineNumber());
            fail=true;
         } else {
            material.set_keyword(token);
            material.set_value(value);
            material.set_lineNumber(lineNumber);
            material.set_loaded(true);
         }
         recognized++;
      }

      if (token.compare("path") == 0) {
         if (found_first_path) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR211: Extraneous path= statement at line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         } else {
            bool reverse=false;
            if (value.substr(0,1).compare("+") == 0) value=value.substr(1);
            else if (value.substr(0,1).compare("-") == 0) {value=value.substr(1); reverse=true;}

            keywordPair *path=new keywordPair();
            path->push_alias("path");
            path->set_keyword(token);
            path->set_value(value);
            path->set_lineNumber(lineNumber);
            path->set_positive_required(false);
            path->set_non_negative_required(false);
            path->set_lowerLimit(0);
            path->set_upperLimit(0);
            path->set_checkLimits(false);
            path->set_loaded(true);

            pathNameList.push_back(path);
            reverseList.push_back(reverse);
         }

         found_first_path=true;
         recognized++;
      }

      if (token.compare("path+") == 0) {
         if (found_first_path) {
            if (value.substr(0,1).compare("+")  == 0 || value.substr(0,1).compare("-") == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR212: Misformatted path at line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            } else {

               keywordPair *path=new keywordPair();
               path->push_alias("path");
               path->set_keyword(token);
               path->set_value(value);
               path->set_lineNumber(lineNumber);
               path->set_positive_required(false);
               path->set_non_negative_required(false);
               path->set_lowerLimit(0);
               path->set_upperLimit(0);
               path->set_checkLimits(false);
               path->set_loaded(true);

               pathNameList.push_back(path);
               reverseList.push_back(false);
            }
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR213: Missing path= statement before line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         recognized++;
      }

      if (token.compare("path-") == 0) {
         if (found_first_path) {
            if (value.substr(0,1).compare("+") == 0 || value.substr(0,1).compare("-") == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR214: Misformatted path at line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
               fail=true;
            } else {

               keywordPair *path=new keywordPair();
               path->push_alias("path");
               path->set_keyword(token);
               path->set_value(value);
               path->set_lineNumber(lineNumber);
               path->set_positive_required(false);
               path->set_non_negative_required(false);
               path->set_lowerLimit(0);
               path->set_upperLimit(0);
               path->set_checkLimits(false);
               path->set_loaded(true);

               pathNameList.push_back(path);
               reverseList.push_back(true);
            }
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR215: Missing path= statement before line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
            fail=true;
         }
         recognized++;
      }

      // should recognize one keyword
      if (recognized != 1) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR216: Unrecognized keyword at line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
         fail=true;
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   return fail;
}

struct EdgeAttribute Boundary::is_segmentOnPath (double x1, double y1, double x2, double y2, vector<Path *> *pathList)
{
   long unsigned int max=-1;
   long unsigned int i=0;
   struct EdgeAttribute attribute;

   attribute.boundary=-1;
   attribute.path=-1;
   attribute.segment=-1;

   while (i < pathIndexList.size()) {
      long unsigned int segment_index=(*pathList)[pathIndexList[i]]->is_segmentOnLine(x1,y1,x2,y2);
      if (segment_index != max) {
         attribute.boundary=get_attribute();
         attribute.path=i;
         attribute.segment=segment_index;
         return attribute;
      }
      i++;
   }
   return attribute;
}

bool Boundary::is_surface_impedance()
{
   if (type.get_value().compare("surface_impedance") == 0) return true;
   return false;
}

bool Boundary::is_perfect_electric_conductor()
{
   if (type.get_value().compare("perfect_electric_conductor") == 0) return true;
   return false;
}

bool Boundary::is_perfect_magnetic_conductor()
{
   if (type.get_value().compare("perfect_magnetic_conductor") == 0) return true;
   return false;
}

bool Boundary::is_mode_voltage()
{
   if (type.get_value().compare("voltage") == 0) return true;
   return false;
}

bool Boundary::is_mode_current()
{
   if (type.get_value().compare("current") == 0) return true;
   return false;
}

bool Boundary::is_boundary()
{
   if (is_surface_impedance()) return true;
   if (is_perfect_electric_conductor()) return true;
   if (is_perfect_magnetic_conductor()) return true;
   return false;
}

bool Boundary::is_modal()
{
   if (mode_block_type.get_value().compare("modal") == 0) return true;
   return false;
}

bool Boundary::is_line()
{
   if (mode_block_type.get_value().compare("line") == 0) return true;
   return false;
}

bool Boundary::is_mode()
{
   if (is_modal()) return true;
   if (is_line()) return true;
   return false;
}

string Boundary::get_block_type()
{
   if (is_modal()) return "Mode";
   if (is_line()) return "Line";
   return "null";
}

string Boundary::get_mode_name()
{
   if (is_modal()) return "mode";
   if (is_line()) return "line";
   return "null";
}

void Boundary::print()
{
   if (is_boundary()) {
      PetscPrintf(PETSC_COMM_WORLD,"Boundary\n");
      PetscPrintf(PETSC_COMM_WORLD,"   name=%s\n",get_name().c_str());
      PetscPrintf(PETSC_COMM_WORLD,"   type=%s\n",get_type().c_str());
      if (is_surface_impedance() && type.is_loaded()) PetscPrintf(PETSC_COMM_WORLD,"   material=%s\n",get_material().c_str());
      long unsigned int i=0;
      while (i < pathNameList.size()) {
         if (i == 0) {
            if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path=-%s\n",pathNameList[i]->get_value().c_str());
            else PetscPrintf(PETSC_COMM_WORLD,"   path=%s\n",pathNameList[i]->get_value().c_str());
         } else {
            if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path-=%s\n",pathNameList[i]->get_value().c_str());
            else PetscPrintf(PETSC_COMM_WORLD,"   path+=%s\n",pathNameList[i]->get_value().c_str());
         }
         i++;
      }
      PetscPrintf(PETSC_COMM_WORLD,"   [attribute=%d]\n",attribute.get_int_value());
      PetscPrintf(PETSC_COMM_WORLD,"EndBoundary\n");
      return;
   }

   if (is_mode()) {
      if (is_modal()) {
         PetscPrintf(PETSC_COMM_WORLD,"Mode\n");
         PetscPrintf(PETSC_COMM_WORLD,"   mode=%d\n",get_mode());
      }
      if (is_line()) {
         PetscPrintf(PETSC_COMM_WORLD,"Line\n");
         PetscPrintf(PETSC_COMM_WORLD,"   line=%d\n",get_mode());
      }
      PetscPrintf(PETSC_COMM_WORLD,"   type=%s\n",get_type().c_str());
      long unsigned int i=0;
      while (i < pathNameList.size()) {
         if (i == 0) {
            if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path=-%s\n",pathNameList[i]->get_value().c_str());
            else PetscPrintf(PETSC_COMM_WORLD,"   path=%s\n",pathNameList[i]->get_value().c_str());
         } else {
            if (reverseList[i]) PetscPrintf(PETSC_COMM_WORLD,"   path-=%s\n",pathNameList[i]->get_value().c_str());
            else PetscPrintf(PETSC_COMM_WORLD,"   path+=%s\n",pathNameList[i]->get_value().c_str());
         }
         i++;
      }
      PetscPrintf(PETSC_COMM_WORLD,"   [attribute=%d]\n",attribute.get_int_value());
      if (is_modal()) PetscPrintf(PETSC_COMM_WORLD,"EndMode\n");
      if (is_line()) PetscPrintf(PETSC_COMM_WORLD,"EndLine\n");
      return;
   }

   PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Boundary::print() did not find valid type.\n");
   return;
}

bool Boundary::inBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Boundary::check(string *indent)
{
   bool fail=false;

   if (! type.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR217: Block at line %d must specify a type.\n",indent->c_str(),indent->c_str(),startLine);
      return true;
   }

   if (is_boundary()) {
      if (!name.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR218: Boundary block at line %d must specify a name.\n",indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }

      if (mode.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR219: Boundary block at line %d must not specify a mode or line.\n",indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }

      if (pathNameList.size() == 0) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR220: Boundary block at line %d must specify at least one path.\n",indent->c_str(),indent->c_str(),startLine);
         fail=true;
      }

      if (is_surface_impedance()) {
         if (!material.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR221: Boundary block at line %d must specify a material.\n",indent->c_str(),indent->c_str(),startLine);
            fail=true;
         }
      } else if (is_perfect_electric_conductor() || is_perfect_magnetic_conductor()) {
         if (material.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR222: Boundary block at line %d must not specify a material.\n",indent->c_str(),indent->c_str(),startLine);
            fail=true;
         }
      }
   } else if (is_mode()) {

      if (name.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR223: %s block at line %d not must specify a name.\n",indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine);
         fail=true;
      }

      if (! mode.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR224: %s block at line %d must specify a %s.\n",indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine,get_mode_name().c_str());
         fail=true;
      }

      if (material.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR225: %s block at line %d not must specify a material.\n",indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine);
         fail=true;
      }

      if (pathNameList.size() == 0) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR226: %s block at line %d must specify at least one path.\n",indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine);
         fail=true;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR227: type at line %d is invalid.\n",indent->c_str(),indent->c_str(),type.get_lineNumber());
      fail=true;
   }

   return fail;
}

// for currents, paths must form 1 or more closed loops
bool Boundary::check_current_paths (string *indent, vector<Path *> *pathList, bool check_closed_loop)
{
   bool fail=false;

   if (is_mode_current()) {
      vector<bool> closed;
      vector<bool> connectedStart;
      vector<bool> connectedEnd;

      // to keep track of what has been looked at
      long unsigned int i=0;
      while (i < pathIndexList.size()) {
         if ((*pathList)[pathIndexList[i]]->is_closed()) {
            closed.push_back(true);
            connectedStart.push_back(true);
            connectedEnd.push_back(true);
         } else {
            closed.push_back(false);
            connectedStart.push_back(false);
            connectedEnd.push_back(false);
         }
         i++;
      }

      // line up the ends of the open sections
      i=0;
      while (pathIndexList.size() > 0 && i < pathIndexList.size()-1) {
         if (! closed[i]) {

            long unsigned int j=i+1;
            while (j < pathIndexList.size()) {
               if (! closed[j]) {

                  // start to start
                  if ((*pathList)[pathIndexList[i]]->get_startPoint()->point_compare((*pathList)[pathIndexList[j]]->get_startPoint())) {
                     if (! connectedStart[j]) {
                        connectedStart[i]=true;
                        connectedStart[j]=true;
                     } else {
                        PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR228: %s block at line %d topology error at (%g,%g).\n",
                                                     indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine,
                                                     (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_x(),
                                                     (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_y());
                        fail=true;
                     }
                  }

                  // start to end
                  if ((*pathList)[pathIndexList[i]]->get_startPoint()->point_compare((*pathList)[pathIndexList[j]]->get_endPoint())) {
                     if (! connectedEnd[j]) {
                        connectedStart[i]=true;
                        connectedEnd[j]=true;
                     } else {
                        PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR229: %s block at line %d topology error at (%g,%g).\n",
                                                     indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine,
                                                     (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_x(),
                                                     (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_y());
                        fail=true;
                     }
                  }

                  // end to start
                  if ((*pathList)[pathIndexList[i]]->get_endPoint()->point_compare((*pathList)[pathIndexList[j]]->get_startPoint())) {
                     if (! connectedStart[j]) {
                        connectedEnd[i]=true;
                        connectedStart[j]=true;
                     } else {
                        PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR230: %s block at line %d topology error at (%g,%g).\n",
                                                     indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine,
                                                     (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_x(),
                                                     (*pathList)[pathIndexList[j]]->get_startPoint()->get_point_value_y());
                        fail=true;
                     }
                  }

                  // end to end
                  if ((*pathList)[pathIndexList[i]]->get_endPoint()->point_compare((*pathList)[pathIndexList[j]]->get_endPoint())) {
                     if (! connectedEnd[j]) {
                        connectedEnd[i]=true;
                        connectedEnd[j]=true;
                     } else {
                        PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR231: %s block at line %d topology error at (%g,%g).\n",
                                                     indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine,
                                                     (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_x(),
                                                     (*pathList)[pathIndexList[j]]->get_endPoint()->get_point_value_y());
                        fail=true;
                     }
                  }

               }
               j++;
            }
         }
         i++;
      }

      // check for dangling ends
      if (check_closed_loop) {
         i=0;
         while (i < pathIndexList.size()) {
            if (! closed[i]) {
               if (! connectedStart[i]) {
                  PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR232: %s block at line %d topology error with dangling point at (%g,%g).\n",
                                               indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine,
                                               (*pathList)[pathIndexList[i]]->get_startPoint()->get_point_value_x(),
                                               (*pathList)[pathIndexList[i]]->get_startPoint()->get_point_value_y());
                  fail=true;
               }
               if (! connectedEnd[i]) {
                  PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR233: %s block at line %d topology error with dangling point at (%g,%g).\n",
                                               indent->c_str(),indent->c_str(),get_block_type().c_str(),startLine,
                                               (*pathList)[pathIndexList[i]]->get_endPoint()->get_point_value_x(),
                                               (*pathList)[pathIndexList[i]]->get_endPoint()->get_point_value_y());
                  fail=true;
               }
            }
            i++;
         }
      }
   }

   return fail;
}

bool Boundary::checkBoundingBox (Vector *lowerLeft, Vector *upperRight, string *indent, double tol, vector<Path *> *pathList)
{
   bool fail=false;

   long unsigned int i=0;
   while (i < pathIndexList.size()) {
      if ((*pathList)[pathIndexList[i]]->checkBoundingBox(lowerLeft, upperRight, indent, tol)) fail=true;
      i++;
   }
   return fail;
}

// Mark mesh borders with an attribute indexing into the border database.
// The border database includes indexes into the boundary database for both boundaries and modes.
void Boundary::markMeshBoundaries (Mesh *mesh, BorderDatabase *borderDatabase, vector<Path *> *pathList)
{
   long unsigned int max=-1;
   struct EdgeAttribute test_attribute;

   // loop through the mesh border elements
   int i=0;
   while (i < mesh->GetNBE()) {

      // get the border vertices
      Array<int> vertices;
      mesh->GetBdrElementVertices(i,vertices);

      // loop through the vertices and see if they fall on the path of the boundary
      if (vertices.Size() == 2) {
         double *vertex0=mesh->GetVertex(vertices[0]);
         double *vertex1=mesh->GetVertex(vertices[1]);
         test_attribute=is_segmentOnPath(vertex0[0],vertex0[1],vertex1[0],vertex1[1],pathList);

         if (test_attribute.boundary != max) {
            int current_bdrAttribute=mesh->GetBdrAttribute(i);
            int new_bdrAttribute;
            if (is_boundary()) new_bdrAttribute=borderDatabase->add(test_attribute,current_bdrAttribute);
            else new_bdrAttribute=borderDatabase->add(test_attribute,get_mode(),current_bdrAttribute);
            mesh->SetBdrAttribute(i,new_bdrAttribute);
         }
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"ASSERT: Boundary::markMeshBoundaries found boundary edge without two points.\n");
      }
      i++;
   }
   return;
}

Boundary::~Boundary()
{
   long unsigned int i=0;
   while (i < pathNameList.size()) {
      delete pathNameList[i];
      i++;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// SourceFile2
///////////////////////////////////////////////////////////////////////////////////////////

SourceFile2::SourceFile2(int startLine_, int endLine_)
{
   startLine=startLine_;
   endLine=endLine_;

   // name
   name.push_alias("name");
   name.set_loaded(false);
   name.set_positive_required(false);
   name.set_non_negative_required(false);
   name.set_lowerLimit(0);
   name.set_upperLimit(0);
   name.set_checkLimits(false);
}

bool SourceFile2::load(string *indent, inputFile *inputs)
{
   bool fail=false;

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      string token,value,line;
      line=inputs->get_line(lineNumber);
      get_token_pair(&line,&token,&value,&lineNumber,*indent);

      int recognized=0;

      if (name.match_alias(&token)) {
         if (name.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR234: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),indent->c_str(),lineNumber,name.get_lineNumber());
            fail=true;
         } else {
            name.set_keyword(token);
            name.set_value(value);
            name.set_lineNumber(lineNumber);
            name.set_loaded(true);
         }
         recognized++;
      }

      // should recognize one keyword
      if (recognized != 1) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR235: Unrecognized keyword at line %d.\n",indent->c_str(),indent->c_str(),lineNumber);
         fail=true;
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   return fail;
}

bool SourceFile2::inBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool SourceFile2::check(string *indent)
{
   bool fail=false;

   if (!name.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR236: File block at line %d must specify a name.\n",indent->c_str(),indent->c_str(),startLine);
      fail=true;
   }

   return fail;
}

void SourceFile2::print() {
   PetscPrintf(PETSC_COMM_WORLD,"File\n");
   PetscPrintf(PETSC_COMM_WORLD,"   name=%s\n",name.get_value().c_str());
   PetscPrintf(PETSC_COMM_WORLD,"EndFile\n");
   return;
}

///////////////////////////////////////////////////////////////////////////////////////////
// BoundaryDatabase
///////////////////////////////////////////////////////////////////////////////////////////

void BoundaryDatabase::mark_boundaries (Mesh *mesh, BorderDatabase *borderDatabase)
{

   // pre-set the mesh attributes to 1 - assumed PEC if there are BorderDatabase elements
   //                                    Note the MFEM does not allow 0 as an attribute, so starting with 1 instead.
   // The mesh attribute is an index into a BorderDatabase that in turn holds indices into a BoundaryDatabase
   int j=0;
   while (j < mesh->GetNBE()) {
      mesh->SetBdrAttribute(j,1);
      j++;
   }

   // set the attributes to index into a BorderDatabase
   // the BorderDatabase is set up in the process
   long unsigned int i=0;
   while (i < boundaryList.size()) {
      boundaryList[i]->markMeshBoundaries(mesh,borderDatabase,&pathList);
      i++;
   }

   return;
}

bool BoundaryDatabase::findSourceFileBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "File", "EndFile", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            SourceFile2 *newSourceFile2=new SourceFile2(block_start,block_stop);
            sourceFileList.push_back(newSourceFile2);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::findPathBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Path", "EndPath", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Path *newPath=new Path(block_start,block_stop);
            pathList.push_back(newPath);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::findBoundaryBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Boundary", "EndBoundary", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Boundary *newBoundary=new Boundary(block_start,block_stop,"");
            boundaryList.push_back(newBoundary);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::findModeBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Mode", "EndMode", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Boundary *newBoundary=new Boundary(block_start,block_stop,"modal");
            boundaryList.push_back(newBoundary);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::findLineBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Line", "EndLine", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Boundary *newBoundary=new Boundary(block_start,block_stop,"line");
            boundaryList.push_back(newBoundary);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

bool BoundaryDatabase::load(const char *filename, bool check_closed_loop) {
   bool fail=false;
   string line;

   if (strcmp(filename,"") == 0) return false;  // this file is optional

   PetscPrintf(PETSC_COMM_WORLD,"   loading mode definition file \"%s\"\n",filename);

   if (inputs.load(filename)) return true;
   inputs.createCrossReference();

   if (inputs.checkVersion(version_name, version_value)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR237: Version mismatch.  Expecting the first line to be: %s %s\n",
                                   indent.c_str(),indent.c_str(),version_name.c_str(),version_value.c_str());
      return true;
   }

   // Source
   if (findSourceFileBlocks()) {fail=true;}

   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      if (sourceFileList[i]->load(&indent, &inputs)) fail=true;
      i++;
   }

   // Path
   if (findPathBlocks()) fail=true;

   i=0;
   while (i < pathList.size()) {
      if (pathList[i]->load(&indent, &inputs)) fail=true;
      i++;
   }

   // Boundary

   if (findBoundaryBlocks()) fail=true;
   long unsigned int boundaryCount=boundaryList.size();

   i=0;
   while (i < boundaryCount) {
      if (boundaryList[i]->load(&indent, &inputs, true, i)) fail=true;
      i++;
   }

   // Mode
   if (findModeBlocks()) fail=true;

   // Line
   if (findLineBlocks()) fail=true;

   i=boundaryCount;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->load(&indent, &inputs, false, i)) fail=true;
      i++;
   }

   // database checks
   if (check()) fail=true;

   // fill out with path pointers for later convenience
   i=0;
   while (i < boundaryList.size()) {
      long unsigned int j=0;
      while (j < boundaryList[i]->get_pathNameList_size()) {
         bool found=false;
         long unsigned int k=0;
         while (k < pathList.size()) {
            if (pathList[k]->get_name().compare(boundaryList[i]->get_pathName(j)) == 0) {
               boundaryList[i]->push(k);
               found=true;
               break;
            }
            k++;
         }
         if (! found) {
            // errors previously reported
            fail=true;
         }
         j++;
      }
      i++;
   }

   // additional Mode check
   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->check_current_paths(&indent,&pathList,check_closed_loop)) fail=true;
      i++;
   }

   // subdivide the paths to eliminate partial overlaps
   subdivide_paths();

   if (fail) PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR238: Failed to load mode definitions.\n",indent.c_str(),indent.c_str());

   return fail;
}

void BoundaryDatabase::print()
{
   PetscPrintf(PETSC_COMM_WORLD,"%s %s\n",version_name.c_str(),version_value.c_str());

   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      sourceFileList[i]->print();
      i++;
   }

   i=0;
   while (i < pathList.size()) {
      pathList[i]->print();
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      boundaryList[i]->print();
      i++;
   }
}

bool BoundaryDatabase::inBlocks(int lineNumber)
{
   long unsigned int i=0;
   while (i < pathList.size()) {
      if (pathList[i]->inBlock(lineNumber)) return true;
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->inBlock(lineNumber)) return true;
      i++;
   }

   i=0;
   while (i < sourceFileList.size()) {
      if (sourceFileList[i]->inBlock(lineNumber)) return true;
      i++;
   }

   return false;
}


bool BoundaryDatabase::check()
{
   bool fail=false;

   // Source
   if (sourceFileList.size() > 1) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR239: Only one File block is allowed.\n",indent.c_str(),indent.c_str());
      fail=true;
   }

   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      if (sourceFileList[i]->check (&indent)) fail=true;
      i++;
   }

   // Path
   i=0;
   while (i < pathList.size()) {

      // individual block checks
      if (pathList[i]->check(&indent)) fail=true;

      // cross block checks

      // duplicated names
      long unsigned int j=i+1;
      while (pathList.size() > 0 && i < pathList.size()-1 && j < pathList.size()) {
         if (pathList[i]->name_is_loaded() && pathList[j]->name_is_loaded()) {
            if (pathList[i]->get_name().compare(pathList[j]->get_name()) == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR240: name at line %d duplicates the name at line %d.\n",
                                            indent.c_str(),indent.c_str(),pathList[j]->get_name_lineNumber(),pathList[i]->get_name_lineNumber());
               fail=true;
            }
         }
         j++;
      }

      i++;
   }

   // Boundary
   i=0;
   while (i < boundaryList.size()) {

      // individual block checks
      if (boundaryList[i]->check(&indent)) fail=true;

      // cross block checks

      if (boundaryList[i]->is_boundary()) {
         // duplicated names
         long unsigned int j=i+1;
         while (boundaryList.size() > 0 && i < boundaryList.size()-1 && j < boundaryList.size()) {
            if (boundaryList[i]->name_is_loaded() && boundaryList[j]->name_is_loaded()) {
               if (boundaryList[i]->get_name().compare(boundaryList[j]->get_name()) == 0) {
                  PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR241: name at line %d duplicates the name at line %d.\n",
                                               indent.c_str(),indent.c_str(),boundaryList[j]->get_name_lineNumber(),boundaryList[i]->get_name_lineNumber());
                  fail=true;
               }
            }
            j++;
         }
      }

      // paths exist
      long unsigned int j=0;
      while (j < boundaryList[i]->get_pathNameList_size()) {
         bool found=false;
         long unsigned int k=0;
         while (k < pathList.size()) {
            if (boundaryList[i]->get_pathName(j).compare(pathList[k]->get_name()) == 0) {found=true; break;}
            k++;
         }
         if (! found) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR242: path at line %d is not defined by a Path block.\n",
                                         indent.c_str(),indent.c_str(),boundaryList[i]->get_pathName_lineNumber(j));
            fail=true;
         }

         j++;
      }

      // paths are not duplicated
      j=0;
      while (boundaryList[i]->get_pathNameList_size() < 0 && j < boundaryList[i]->get_pathNameList_size()-1) {
         long unsigned int k=j+1;
         while (k < boundaryList[i]->get_pathNameList_size()) {
            if (boundaryList[i]->get_pathName(j).compare(boundaryList[i]->get_pathName(k)) == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR243: path at line %d duplicates the path at line %d.\n",
                                            indent.c_str(),indent.c_str(),boundaryList[i]->get_pathName_lineNumber(k),boundaryList[i]->get_pathName_lineNumber(j));
               fail=true;
            }
            k++;
         }
         j++; 
      }

      i++;
   }

   // only one voltage and one current per mode
   int highestModeNumber=0;
   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->is_mode() && boundaryList[i]->get_mode() > highestModeNumber) highestModeNumber=boundaryList[i]->get_mode();
      i++;
   }

   int modeNumber=1;
   while (modeNumber <= highestModeNumber) {
      bool foundVoltage=false;
      bool foundCurrent=false;
      long unsigned int i=0;
      while (i < boundaryList.size()) {
         if (boundaryList[i]->is_mode() && boundaryList[i]->get_mode() == modeNumber) {
            if (boundaryList[i]->is_mode_voltage()) {
               if (!foundVoltage) foundVoltage=true;
               else {
                  PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR244: %s block at line %d makes an extraneous voltage definition for mode %d.\n",
                                               indent.c_str(),indent.c_str(),boundaryList[i]->get_block_type().c_str(),boundaryList[i]->get_startLine(),boundaryList[i]->get_mode());
                  fail=true;
               }
            }

            if (boundaryList[i]->is_mode_current()) {
               if (!foundCurrent) foundCurrent=true;
               else {
                  PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR245: %s block at line %d makes an extraneous current definition for mode %d.\n",
                                               indent.c_str(),indent.c_str(),boundaryList[i]->get_block_type().c_str(),boundaryList[i]->get_startLine(),boundaryList[i]->get_mode());
                  fail=true;
               }
            }
         }
         i++;
      }
      modeNumber++;
   }

   // mode numbering
   // check voltage and current separately
   vector<int> modeNumbers_voltage;
   vector<int> modeNumbers_current;
   i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->is_mode_voltage()) modeNumbers_voltage.push_back(0);
      if (boundaryList[i]->is_mode_current()) modeNumbers_current.push_back(0);
      i++;    
   }
   modeNumbers_voltage.push_back(0);  // one more since the mode numbers are not zero based
   modeNumbers_current.push_back(0);  // one more since the mode numbers are not zero based

   i=0;
   while (i < boundaryList.size()) {

      // voltage
      if (boundaryList[i]->is_mode_voltage() && boundaryList[i]->mode_is_loaded()) {
         if (boundaryList[i]->get_mode() < (int)modeNumbers_voltage.size()) {
            modeNumbers_voltage[boundaryList[i]->get_mode()]++;
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR246: %s number at line number %d is too large.\n",
                                         indent.c_str(),indent.c_str(),boundaryList[i]->get_block_type().c_str(),boundaryList[i]->get_mode_lineNumber());
            fail=true;
         }
      }

      // current
      if (boundaryList[i]->is_mode_current() && boundaryList[i]->mode_is_loaded()) {
         if (boundaryList[i]->get_mode() < (int)modeNumbers_current.size()) {
            modeNumbers_current[boundaryList[i]->get_mode()]++;
         } else {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR246: %s number at line number %d is too large.\n",
                                         indent.c_str(),indent.c_str(),boundaryList[i]->get_block_type().c_str(),boundaryList[i]->get_mode_lineNumber());
            fail=true;
         }
      }

      i++;
   }

   // voltage
   long unsigned int lastNonZero_voltage=modeNumbers_voltage.size()-1;
   while (modeNumbers_voltage.size() > 0 && lastNonZero_voltage > 0) {
      if (modeNumbers_voltage[lastNonZero_voltage] > 0) break;
      lastNonZero_voltage--;
   }

   // current
   long unsigned int lastNonZero_current=modeNumbers_current.size()-1;
   while (modeNumbers_current.size() > 0 && lastNonZero_current > 0) {
      if (modeNumbers_current[lastNonZero_current] > 0) break;
      lastNonZero_current--;
   }   

   // voltage
   i=1;
   while (i <= lastNonZero_voltage) {
      if (modeNumbers_voltage[i] == 0) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR247: Modes/Lines are not sequentially numbered starting at 1.\n",indent.c_str(),indent.c_str());
         fail=true;
         break;
      }
      i++;
   }

   // current
   i=1;
   while (i <= lastNonZero_current) {
      if (modeNumbers_current[i] == 0) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR247: Modes/Lines are not sequentially numbered starting at 1.\n",indent.c_str(),indent.c_str());
         fail=true;
         break;
      }
      i++;
   }

   // check for extraneous text
   i=1;  // skip the first line, which is the version information
   while (i < inputs.get_size()) {
      if (! inBlocks(inputs.get_lineNumber(i))) {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR248: Invalid input at line %d.\n",indent.c_str(),indent.c_str(),inputs.get_lineNumber(i));
         fail=true;
      }
      i++;
   }

   return fail;
}

// Make sure the paths fall within the bounding box
// to catch scaling errors.
bool BoundaryDatabase::check_scale (Mesh *mesh, int order)
{
   bool fail=false;

   Vector lowerLeft,upperRight;
   mesh->GetBoundingBox(lowerLeft,upperRight,max(order,1));

   long unsigned int i=0;
   while (i < boundaryList.size()) {
      if (boundaryList[i]->checkBoundingBox(&lowerLeft,&upperRight,&indent,tol,&pathList)) fail=true;
      i++;
   }

   return fail;
}

// remove overlaps in paths
void BoundaryDatabase::subdivide_paths ()
{
   long unsigned int i=0;
   while (i < pathList.size()) {
      long unsigned int j=0;
      while (j < pathList.size()) {
         if (i != j) pathList[i]->subdivide(pathList[j]);
         j++;
      }
      i++;
   }
}

BoundaryDatabase::~BoundaryDatabase ()
{
   long unsigned int i=0;
   while (i < sourceFileList.size()) {
      delete sourceFileList[i];
      i++;
   }

   i=0;
   while (i < boundaryList.size()) {
      delete boundaryList[i];
      i++;
   }

   i=0;
   while (i < pathList.size()) {
      delete pathList[i];
      i++;    
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Border
///////////////////////////////////////////////////////////////////////////////////////////

Border::Border()
{
   boundary.boundary=-1;
   boundary.path=-1;
   boundary.segment=-1;
}

bool Border::has_boundary()
{
   long unsigned int max=-1;
   if (boundary.boundary == max) return false;
   return true;
}

// any mode
bool Border::has_mode()
{
   if (mode.size() > 0) return true;
   return false;
}

struct EdgeAttribute Border::get_mode(long unsigned int i)
{
   struct EdgeAttribute temp;
   temp.boundary=-1;
   temp.path=-1;
   temp.segment=-1;

   if (i >= mode.size()) return temp;
   return mode[i];
}

// specified modeNumber
bool Border::has_mode(long unsigned int modeNumber)
{
   long unsigned int max=-1;

   if (mode.size() < modeNumber+1) return false;
   if (mode[modeNumber].boundary == max) return false;
   return true;
}

// for boundary
void Border::set(struct EdgeAttribute attribute)
{
   boundary=attribute;
}

// for mode 
void Border::set(struct EdgeAttribute attribute, long unsigned int modeNumber)
{
   struct EdgeAttribute temp;
   temp.boundary=-1;
   temp.path=-1;
   temp.segment=-1;

   // make sure there is enough space
   long unsigned int i=mode.size();
   while (i <= modeNumber) {
      mode.push_back(temp);
      i++;
   }

   mode[modeNumber]=attribute;
}

// for boundaries
bool Border::exists(struct EdgeAttribute test_attribute)
{
   if (boundary == test_attribute) return true;
   return false;
}

// for modes
bool Border::exists(struct EdgeAttribute test_attribute, long unsigned int modeNumber)
{
   if (! (modeNumber < mode.size())) return false;
   if (mode[modeNumber] == test_attribute) return true;
   return false;
}

// for modes
bool Border::exists_any_mode(struct EdgeAttribute test_attribute)
{
   long unsigned int i=0;
   while (i < mode.size()) {
      if (mode[i] == test_attribute) return true;
      i++;
   }
   return false;
}

// for modes
bool Border::exists_test_boundary_only(struct EdgeAttribute test_attribute)
{
   long unsigned int i=0;
   while (i < mode.size()) {
      if (mode[i].boundary == test_attribute.boundary) return true;
      i++;
   }
   return false;
}

void Border::print()
{
   cout << "Border=" << this << endl;

   cout << "   boundary.boundary=" << boundary.boundary 
        << "   boundary.path=" << boundary.path
        << "   boundary.segment=" << boundary.segment << endl;

   long unsigned int i=0;
   while (i < mode.size()) {
      cout << "   mode[" << i << "].boundary=" << mode[i].boundary
           << "   mode[" << i << "].path=" << mode[i].path
           << "   mode[" << i << "].segment=" << mode[i].segment << endl;
      i++;
   }

}

///////////////////////////////////////////////////////////////////////////////////////////
// BorderDatabase
///////////////////////////////////////////////////////////////////////////////////////////

BorderDatabase::BorderDatabase()
{
   // MFEM does not support border attributes of 0 on the mesh, so create a dummy entry for the 0 slot.
   // This dummy border is never accessed.
   Border *border0=new Border;
   borderList.push_back(border0);

   // The attribute of 1 for the mesh border attributes indicates the default PEC boundary condition,
   // so create another dummy entry for the 1 slot the points to nothing.  That will yield the PEC boundary.
   Border *border1=new Border;
   borderList.push_back(border1);
}

// for boundaries
long unsigned int BorderDatabase::exists(struct EdgeAttribute test_attribute)
{
   long unsigned int i=0;
   while (i < borderList.size()) {
      if (borderList[i]->exists(test_attribute)) return i;
      i++;
   }
   return -1;
}

// for modes
long unsigned int BorderDatabase::exists(struct EdgeAttribute test_attribute, long unsigned int modeNumber)
{
   long unsigned int i=0;
   while (i < borderList.size()) {
      if (borderList[i]->exists(test_attribute,modeNumber)) return i;
      i++;
   }
   return -1;
}

// for modes
long unsigned int BorderDatabase::exists_any_mode(struct EdgeAttribute test_attribute)
{
   long unsigned int i=0;
   while (i < borderList.size()) {
      if (borderList[i]->exists_any_mode(test_attribute)) return i;
      i++;
   }
   return -1;
}

// for modes
long unsigned int BorderDatabase::exists_test_boundary_only(struct EdgeAttribute test_attribute)
{
   long unsigned int i=0;
   while (i < borderList.size()) {
      if (borderList[i]->exists_test_boundary_only(test_attribute)) return i;
      i++;
   }
   return -1;
}

// for boundaries
long unsigned int BorderDatabase::add(struct EdgeAttribute edgeAttribute, int existing_bdr_attribute)
{
   long unsigned int max=-1;
   long unsigned int border_index=-1;

   if (existing_bdr_attribute == 1) {
      // see if it exists already
      border_index=exists(edgeAttribute);
      if (border_index == max) {
         // does not exist
         // create a new slot
         Border *border=new Border;
         border->set(edgeAttribute);
         borderList.push_back(border);
         border_index=borderList.size()-1;
      } else {
         // do nothing - already here
      }
   } else {
      // attribute already assigned
      // save here
      border_index=existing_bdr_attribute;
      borderList[existing_bdr_attribute]->set(edgeAttribute);
   }
   return border_index;
}

// for modes
long unsigned int BorderDatabase::add(struct EdgeAttribute edgeAttribute, long unsigned int modeNumber, int existing_bdr_attribute)
{
   long unsigned int max=-1;
   long unsigned int border_index=-1;

   if (existing_bdr_attribute == 1) {
      // see if this case exists already
      border_index=exists(edgeAttribute,modeNumber);
      if (border_index == max) {
         // does not exist
         // see if it exists for other mode numbers
         border_index=exists_any_mode(edgeAttribute);
         if (border_index == max) {
            // does not exist
            // create a new slot
            Border *border=new Border;
            border->set(edgeAttribute,modeNumber);
            borderList.push_back(border);
            border_index=borderList.size()-1;
         } else {
            // exists with other mode numbers
            // save here with new mode number
            border_index=existing_bdr_attribute;
            borderList[existing_bdr_attribute]->set(edgeAttribute,modeNumber);
         }
      } else {
         // exists - do nothing
      }
   } else {
      // attribute already assigned

      // see if the mode space is empty
      if (borderList[existing_bdr_attribute]->has_mode()) {
         // something is here
         // see if this case already exists
         if (borderList[existing_bdr_attribute]->exists(edgeAttribute,modeNumber)) {
            // exists - do nothing
            border_index=existing_bdr_attribute;
         } else {
            // save here with new mode
            border_index=existing_bdr_attribute;
            borderList[existing_bdr_attribute]->set(edgeAttribute,modeNumber);
         }
      } else {
         // empty
         // save here
         border_index=existing_bdr_attribute;
         borderList[existing_bdr_attribute]->set(edgeAttribute,modeNumber);
      }
   }
   return border_index;
}

// return an Array that includes the borderDatabase entries that reference a given boundary
Array<int>* BorderDatabase::build_entry_list (bool is_boundary, long unsigned int test_index, int modeNumber)
{
   Array<int> *entry_list=new Array<int>;

   long unsigned int i=0;
   while (i < borderList.size()) {

      if (is_boundary && borderList[i]->has_boundary() && borderList[i]->get_boundary().boundary == test_index) {
         entry_list->Append(i);
      }

      if (! is_boundary && borderList[i]->has_mode() && borderList[i]->get_mode(modeNumber).boundary == test_index) {
         entry_list->Append(i);
      }
      
      i++;
   }

   return entry_list;
}

void BorderDatabase::print()
{
   long unsigned int i=0;
   while (i < borderList.size()) {
      borderList[i]->print();
      i++;
   }
}

BorderDatabase::~BorderDatabase()
{
   long unsigned int i=0;
   while (i < borderList.size()) {
      delete borderList[i];
      i++;
   }
}


// last used ERROR is 250
