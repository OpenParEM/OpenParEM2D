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

#include "OpenParEMmaterials.hpp"

bool double_compare (double a, double b, double tol)
{
   if (a == b) return true;
   if (a == 0 && fabs(b) < tol) return true;
   if (b == 0 && fabs(a) < tol) return true;
   if (fabs((a-b)/a) < tol) return true;
   return false;
}

bool is_double (const char *b) {
   bool foundPeriod=false;
   bool foundE=false;
   bool foundChar=false;
   size_t length;

   // skip trailing white space
   size_t i=strlen(b)-1;
   while (i >= 0) {
      if (b[i] != ' ' && b[i] != '\t') break;
      i--;
   }
   length=i+1;

   // skip leading whitespace
   i=0;
   while (i < length) {
      if (b[i] != ' ' && b[i] != '\t') break;
      i++;
   }

   // skip leading sign
   if (b[i] == '+' || b[i] == '-') i++;

   while (i < length) {
      if (b[i] == '.') {
         if (foundPeriod) return false;
         foundPeriod=true;
      }
      else if (b[i] == ' ' || b[i] == '\t') return false;
      else if (i > 0 && !foundE && (b[i] == '+' || b[i] == '-')) return false;
      else if (isalpha(b[i])) {
         foundChar=true;
         if (b[i] == 'e' || b[i] == 'E') {
            if (foundE) return false;
            foundE=true;

            if (i < length-1 && isdigit(b[i+1])) {foundChar=false; i++;}
            if (i < length-2 && b[i+1] == '-' && isdigit(b[i+2])) {foundChar=false; i+=2;}
            if (i < length-2 && b[i+1] == '+' && isdigit(b[i+2])) {foundChar=false; i+=2;}
         }
      }
      if (foundChar) {return false;}
      i++;
   }
   return true;
}

bool is_double (string *a) {
   return is_double((*a).c_str());
}

bool is_int (string *a) {
   const char *b=(*a).c_str();

   // skip leading whitespace
   size_t i=0;
   while (i < strlen(b)) {
      if (b[i] != ' ' && b[i] != '\t') break;
      i++;
   }

   while (i < strlen(b)) {
      if (! isdigit(b[i])) return false;
      i++;
   }
   return true;
}

bool is_point (string *a)
{
   const char *b=(*a).c_str();

   // count (, ), and comma
   size_t countOpen=0;
   size_t countClose=0;
   size_t countComma=0;
   size_t i=0;
   while (i < strlen(b)) {
      if (b[i] == '(') countOpen++;
      if (b[i] == ')') countClose++;
      if (b[i] == ',') countComma++;
      i++;
   }

   if (countOpen != 1) return false;
   if (countClose != 1) return false;
   if (countComma != 1) return false;

   // check for leading characters
   i=0;
   while (i < strlen(b)) {
      if (b[i] != ' ' && b[i] != '\t') {
         if (b[i] == '(') break;
         return false;
      }
      i++;
   }

   // find the numbers
   size_t allocSize=256;
   char x[allocSize],y[allocSize];
   size_t size_x=0, size_y=0;
   bool foundOpen=false;
   bool foundComma=false;
   while (i < strlen(b)) {
      if (b[i] == '(' || b[i] == ',' || b[i] == ')') {
         if (b[i] == '(') foundOpen=true;
         if (b[i] == ',') foundComma=true;
         if (b[i] == ')') break;
      } else {
         if (foundOpen && ! foundComma) {
            x[size_x]=b[i];
            size_x++;
            if (size_x == allocSize) return false;
         }
         if (foundOpen && foundComma) {
            y[size_y]=b[i];
            size_y++;
            if (size_y == allocSize) return false;
         }
      }
      i++;
   }

   // complete the strings
   x[size_x]='\0';
   y[size_y]='\0';

   // check for trailing characters
   i++;
   while (i < strlen(b)) {
      if (b[i] != ' ' && b[i] != '\t') return false;
      i++;
   }

   if (! is_double(x)) return false;
   if (! is_double(y)) return false;

   return true;
}

bool point_get_xy (string *a, double *x_value, double *y_value, string indent, int lineNumber_)
{
   const char *b=(*a).c_str();

   size_t allocSize=256;
   char x[allocSize],y[allocSize];
   size_t size_x=0, size_y=0;
   bool foundOpen=false;
   bool foundComma=false;
   size_t i=0;
   while (i < strlen(b)) {
      if (b[i] == '(' || b[i] == ',' || b[i] == ')') {
         if (b[i] == '(') foundOpen=true;
         if (b[i] == ',') foundComma=true;
         if (b[i] == ')') break;
      } else {
         if (foundOpen && ! foundComma) {
            x[size_x]=b[i];
            size_x++;
            if (size_x == allocSize) return false;
         }
         if (foundOpen && foundComma) {
            y[size_y]=b[i];
            size_y++;
            if (size_y == allocSize) return false;
         }
      }
      i++;
   }

   // complete the strings
   x[size_x]='\0';
   y[size_y]='\0';

   // convert to strings to use stod
   string x_str=x;
   string y_str=y;

   // get the values
   try {*x_value=stod(x_str);}
   catch (const std::invalid_argument& ia) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR500: %s value at line %d is invalid.\n",indent.c_str(),indent.c_str(),a->c_str(),lineNumber_);
      return true;
   }

   // get the values
   try {*y_value=stod(y_str);}
   catch (const std::invalid_argument& ia) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR501: %s value at line %d is invalid.\n",indent.c_str(),indent.c_str(),a->c_str(),lineNumber_);
      return true;
   }

   return false;
}

void get_token_pair (string *line, string *token, string *value, int *lineNumber, string indent) {
   string test;
   int count=0;

   stringstream ss(*line);
   while (getline(ss,test,'=')) {
      if (count == 0) *token=test;
      else if (count == 1) *value=test;
      else {
         PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR502: Incorrectly formatted entry at line %d.\n",indent.c_str(),indent.c_str(),*lineNumber);
      }
      count++;
   }

   if (count < 2) PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR503: Missing value at line %d.\n",indent.c_str(),indent.c_str(),*lineNumber);

   // chop off white space from front and back

   size_t first=(*token).find_first_not_of(" \t");
   size_t last=(*token).find_last_not_of(" \t");
   *token=(*token).substr(first,last-first+1);

   if (value->length() > 0) {
      first=(*value).find_first_not_of(" \t");
      last=(*value).find_last_not_of(" \t");
      *value=(*value).substr(first,last-first+1);
   }

   return;
}

///////////////////////////////////////////////////////////////////////////////////////////
// keywordPair
///////////////////////////////////////////////////////////////////////////////////////////

bool keywordPair::int_limit_checks (string *keyword, int lineNumber)
{
   if (positive_required && int_value <= 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR504: %s at line %d is required to be positive.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (non_negative_required && int_value < 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR505: %s at line %d is required to be non-negative.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (int_value < lowerLimit) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR506: %s at line %d is required to be >= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,lowerLimit);
      return false;
   }

   if (int_value > upperLimit) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR507: %s at line %d is required to be <= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,upperLimit);
      return false;
   }

   return true;
}

bool keywordPair::dbl_limit_checks (string *keyword, int lineNumber)
{
   if (positive_required && dbl_value <= 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR508: %s at line %d is required to be positive.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false; 
   }

   if (non_negative_required && dbl_value < 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR509: %s at line %d is required to be non-negative.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (dbl_value < lowerLimit*(1-dbl_tolerance)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR510: %s at line %d is required to be >= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,lowerLimit);
      return false;
   }

   if (dbl_value > upperLimit*(1+dbl_tolerance)) { 
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR511: %s at line %d is required to be <= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,upperLimit);
      return false;
   }

   return true;
}

bool keywordPair::point_limit_checks (string *keyword, int lineNumber)
{
   if (positive_required && (point_value.x <= 0 || point_value.y <= 0)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR512: %s at line %d is required to be positive.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (non_negative_required && (point_value.x < 0 || point_value.y < 0)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR513: %s at line %d is required to be non-negative.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber);
      return false;
   }

   if (point_value.x < lowerLimit*(1-dbl_tolerance) || point_value.y < lowerLimit*(1-dbl_tolerance)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR514: %s at line %d s required to be >= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,lowerLimit);
      return false;
   }

   if (point_value.x > upperLimit*(1+dbl_tolerance) || point_value.y > upperLimit*(1+dbl_tolerance)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR515: %s at line %d is required to be <= %g.\n",
                                   indent.c_str(),indent.c_str(),keyword->c_str(),lineNumber,upperLimit);
      return false;
   }

   return true;
}

bool keywordPair::limit_check (string type)
{
   bool fail=false;

   if (type.compare("int") == 0) {
      if (! int_limit_checks (&keyword, lineNumber)) fail=true;
   } else if (type.compare("double") == 0) {
      if (! dbl_limit_checks (&keyword, lineNumber)) fail=true;
   } else if (type.compare("point") == 0) {
      if (! point_limit_checks (&keyword, lineNumber)) fail=true;
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ASSERT: bad selection in keywordPair::limit_check\n");
   }

   return fail;
}


bool keywordPair::match_alias (string *token)
{
   long unsigned int i=0;
   while (i < aliases.size()) {
      if (aliases[i].compare(*token) == 0) return true;
      i++;
   }
   return false;
}

bool keywordPair::loadBool (string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR516: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check for a boolean
   if (value_->compare("true") != 0 && value_->compare("false") != 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR517: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get the value
   if (value_->compare("true") == 0) bool_value=true;
   else bool_value=false;

   // save it
   keyword=*token;
   value=*value_;
   lineNumber=lineNumber_;
   loaded=true;

   return false;
}

bool keywordPair::loadInt (string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR518: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check for a pure number
   if (!is_int(value_)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR519: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get the value
   try {int_value=stoi(*value_);}
   catch (const std::invalid_argument& ia) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR520: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // check the limits
   if (checkLimits && ! int_limit_checks (token, lineNumber_)) return true;

   // save it
   keyword=*token;
   value=*value_;
   lineNumber=lineNumber_;
   loaded=true;

   return false;
}

bool keywordPair::loadDouble (string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR521: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check for a pure number
   if (!is_double(value_)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR522: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get the value
   try {dbl_value=stod(*value_);}
   catch (const std::invalid_argument& ia) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR523: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // check the limits
   if (checkLimits && ! dbl_limit_checks (token, lineNumber_)) return true;

   // save it
   keyword=*token;
   value=*value_;
   lineNumber=lineNumber_;
   loaded=true;

   return false;
}

bool keywordPair::loadPoint (string *token, string *value_, int lineNumber_)
{
   // check for duplicate
   if (loaded) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR524: Duplicate entry at line %d for previous entry at line %d.\n",
                                   indent.c_str(),indent.c_str(),lineNumber_,lineNumber);
      return true;
   }

   // check
   if (!is_point(value_)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR525: %s value at line %d is invalid.\n",
                                   indent.c_str(),indent.c_str(),token->c_str(),lineNumber_);
      return true;
   }

   // get values
   if (point_get_xy (value_, &point_value.x, &point_value.y, indent, lineNumber_)) return true;

   // check the limits
   if (checkLimits && ! point_limit_checks (token, lineNumber_)) return true;

   loaded=true;

   return false;
}

bool keywordPair::dbl_compare (keywordPair *test)
{
   if (dbl_value == test->dbl_value) return true;
   if (dbl_value == 0 && fabs(test->dbl_value) < dbl_tolerance) return true;
   if (test->dbl_value == 0 && fabs(dbl_value) < dbl_tolerance) return true;
   if (fabs((dbl_value-test->dbl_value)/dbl_value) < dbl_tolerance) return true;
   return false;
}

bool keywordPair::value_compare (keywordPair *test)
{
   if (value.compare(test->value) == 0) return true;
   return false;
}

bool keywordPair::point_compare (keywordPair *a)
{
   if (! double_compare(point_value.x,a->point_value.x,dbl_tolerance)) return false;
   if (! double_compare(point_value.y,a->point_value.y,dbl_tolerance)) return false;
   return true;
}

double keywordPair::get_point_distance (keywordPair *a)
{
   return sqrt(pow(get_point_value_x()-a->get_point_value_x(),2)+pow(get_point_value_y()-a->get_point_value_y(),2));
}

bool keywordPair::is_close_point (keywordPair *a)
{
   if (! double_compare(get_point_value_x(),a->get_point_value_x(),1e-12)) return false;
   if (! double_compare(get_point_value_y(),a->get_point_value_y(),1e-12)) return false;
   return true;
}

void keywordPair::copy (keywordPair a)
{
   aliases.clear();
   long unsigned int i=0;
   while (i < a.aliases.size()) {
      aliases.push_back(a.aliases[i]);
      i++;
   }

   keyword=a.keyword;
   value=a.value;
   lineNumber=a.lineNumber;
   int_value=a.int_value;
   dbl_value=a.dbl_value;
   bool_value=a.bool_value;
   point_value=a.point_value;
   loaded=a.loaded;
   lowerLimit=a.lowerLimit;
   upperLimit=a.upperLimit;
   positive_required=a.positive_required;
   non_negative_required=a.non_negative_required;
   indent=a.indent;
   dbl_tolerance=a.dbl_tolerance;
   checkLimits=a.checkLimits;
}

keywordPair* keywordPair::clone ()
{
   keywordPair *b=new keywordPair();
   b->copy(*this);
   return b;
}

void keywordPair::print()
{
  long unsigned int i=0;
  while (i < aliases.size()) {
     PetscPrintf(PETSC_COMM_WORLD,"alias: %s\n",aliases[i].c_str());
     i++;
  }
  PetscPrintf(PETSC_COMM_WORLD,"keyword: %s\n",keyword.c_str());
  PetscPrintf(PETSC_COMM_WORLD,"value: %s\n",value.c_str());
  PetscPrintf(PETSC_COMM_WORLD,"lineNumber: %d\n",lineNumber);
  PetscPrintf(PETSC_COMM_WORLD,"int_value: %d\n",int_value);
  PetscPrintf(PETSC_COMM_WORLD,"dbl_value: %g\n",dbl_value);
  PetscPrintf(PETSC_COMM_WORLD,"loaded: %d\n",loaded);
  PetscPrintf(PETSC_COMM_WORLD,"lowerLimit: %g\n",lowerLimit);
  PetscPrintf(PETSC_COMM_WORLD,"upperLimit: %g\n",upperLimit);
  PetscPrintf(PETSC_COMM_WORLD,"postive_required: %d\n",positive_required);
  PetscPrintf(PETSC_COMM_WORLD,"non_negative_required: %d\n",non_negative_required);
  PetscPrintf(PETSC_COMM_WORLD,"indent: [%s]\n",indent.c_str());
  PetscPrintf(PETSC_COMM_WORLD,"dbl_tolerance: %g\n",dbl_tolerance);
  PetscPrintf(PETSC_COMM_WORLD,"checkLimits: %d\n",checkLimits);
}

///////////////////////////////////////////////////////////////////////////////////////////
// inputFile
///////////////////////////////////////////////////////////////////////////////////////////

// return true on fail
bool inputFile::load(const char *filename)
{
   int lineNumber=0;
   string line;

   if (strcmp(filename,"") == 0) return true;

   ifstream materialFile;
   materialFile.open(filename,ifstream::in);

   if (materialFile.is_open()) {

      while (getline(materialFile,line)) {
         lineNumber++;

         // skip blank lines
         if (line.compare("") != 0) {

            // skip comment lines
            if (! is_comment(line)) {

               // chop off comments
               line=line.substr(0,line.find("//",0));

               // chop off white space from front and back
               size_t first=line.find_first_not_of(" \t");
               size_t last=line.find_last_not_of(" \t");

               if (first <= last) {
                  line=line.substr(first,last-first+1);

                  // save the line and the line number
                  lineTextList.push_back(line);
                  lineNumberList.push_back(lineNumber);
               }
            }
         }
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR526: File \"%s\" could not be opened for reading.\n",filename);
      return true;
   }
   return false;
}

void inputFile::createCrossReference()
{
   // initialize
   long unsigned int max=lineNumberList[lineNumberList.size()-1];
   long unsigned int i=0;
   while (i <= max) {
      crossReferenceList.push_back(-1);
      i++;
   }

   // index
   i=0;
   while (i < lineNumberList.size()) {
      crossReferenceList[lineNumberList[i]]=i;
      i++;
   }
}

// return true on fail
bool inputFile::checkVersion(string name, string version)
{
   string line,token,value,test;

   if (lineTextList.size() == 0) return true;

   line=lineTextList[0];

   // chop off comments 
   line=line.substr(0,line.find("//",0));

   // chop off white space from front and back
   size_t first=line.find_first_not_of(" \t");
   size_t last=line.find_last_not_of(" \t");
   line=line.substr(first,last-first+1);

   // name
   first=line.find_first_not_of(" \t");
   last=line.find_first_of(" \t");
   token=line.substr(first,last-first);

   if (token.compare(name) != 0) return true;

   // version
   first=line.find_first_of(" \t");
   last=line.find_last_not_of(" \t");
   value=line.substr(first+1,last-first+1);

   first=value.find_first_not_of(" \t");
   last=value.find_last_not_of(" \t");
   value=value.substr(first,last-first+1);

   if (value.compare(version) != 0) return true;

   return false;
}



int inputFile::get_first_lineNumber(){
   if (lineNumberList.size() > 0) return lineNumberList[0];
   return -1;
}

int inputFile::get_last_lineNumber() {
   if (lineNumberList.size() > 0) return lineNumberList[lineNumberList.size()-1];
   return -1;
}

int inputFile::get_previous_lineNumber (int lineNumber) {
   long unsigned int i=lineNumberList.size()-1;
   while (i > 0) {
      if (lineNumberList[i] == lineNumber) return lineNumberList[i-1];
      i--;
   }
   return lineNumberList[i];
}

int inputFile::get_next_lineNumber (int lineNumber) {
   long unsigned int i=0;
   while (i < lineNumberList.size()-1) {
      if (lineNumberList[i] == lineNumber) return lineNumberList[i+1];
      i++;
   }
   return lineNumberList[i];
}

string inputFile::get_line (int lineNumber) {
   return lineTextList[crossReferenceList[lineNumber]];
}

// find the starting and stopping line numbers for a block inclusive of the keywords
bool inputFile::findBlock(int search_startLine, int search_stopLine, 
                         int *block_startLine, int *block_stopLine,
                         string initiator, string terminator, bool reportUnmatchedText)
{
   bool fail=false;
   *block_startLine=-1;
   *block_stopLine=-1;

   // skip to the search range
   long unsigned int i=0;
   while (i < lineNumberList.size()) {
      if (lineNumberList[i] >= search_startLine) break;
      i++;
   }

   // find the start of the block
   while (i < lineNumberList.size()) {
      if (lineNumberList[i] <= search_stopLine) {
         if (lineTextList[i].compare(initiator) == 0) {
            *block_startLine=lineNumberList[i];
            break;
         } else if (lineTextList[i].compare(terminator) == 0) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR527: \"%s\" found at line %d is missing an opening \"%s\" keyword.\n",
                                         indent.c_str(),indent.c_str(),terminator.c_str(),lineNumberList[i],initiator.c_str());
            *block_stopLine=lineNumberList[i];
            return true;
         } else {
            if (reportUnmatchedText) {
               PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR528: Invalid entry at line %d.\n",indent.c_str(),indent.c_str(),lineNumberList[i] );
               fail=true;
            }
         }
      }
      i++;
   }
   if (fail) return true;

   // no block found
   if (*block_startLine < 0) {
      *block_stopLine=search_stopLine;
      return false;
   }

   // find the end of the block
   i++;
   while (i < lineNumberList.size()) {
      if (lineNumberList[i] <= search_stopLine) {
         if (lineTextList[i].compare(terminator) == 0) {
            *block_stopLine=lineNumberList[i];
            break;
         }

         // check for missing terminator
         if (lineTextList[i].compare(initiator) == 0) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR529: \"%s\" block at line %d is incorrectly terminated at line %d.\n",
                                         indent.c_str(),indent.c_str(),initiator.c_str(),*block_startLine,lineNumberList[i]);
            *block_stopLine=get_previous_lineNumber(lineNumberList[i]);
            return true;
         }
      }
      i++;
   }

   // missing block terminator
   if (*block_stopLine < 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR530: \"%s\" block at line %d is missing its terminator \"%s\".\n",
                                   indent.c_str(),indent.c_str(),initiator.c_str(),*block_startLine,terminator.c_str());
      *block_stopLine=search_stopLine;
      return true;
   }

   return false;
}

void inputFile::print()
{
   long unsigned int i=0;
   while (i < lineTextList.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s\n",lineNumberList[i],lineTextList[i].c_str());
      i++;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Frequency
///////////////////////////////////////////////////////////////////////////////////////////

Frequency::Frequency (int startLine_, int endLine_, bool checkLimits_)
{
   startLine=startLine_;
   endLine=endLine_;

   // frequency

   frequency.push_alias("frequency");
   frequency.push_alias("freq");
   frequency.push_alias("f");
   frequency.set_loaded(false);
   frequency.set_positive_required(true);
   frequency.set_non_negative_required(false);
   frequency.set_lowerLimit(0);
   frequency.set_upperLimit(1e12);
   frequency.set_checkLimits(checkLimits_);

   // relative permittivity

   relative_permittivity.push_alias("relative_permittivity");
   relative_permittivity.push_alias("er");
   relative_permittivity.push_alias("epsr");
   relative_permittivity.set_loaded(false);
   relative_permittivity.set_positive_required(true);
   relative_permittivity.set_non_negative_required(false);
   relative_permittivity.set_lowerLimit(1);
   relative_permittivity.set_upperLimit(1e6);
   relative_permittivity.set_checkLimits(checkLimits_);

   // relative permeability

   relative_permeability.push_alias("relative_permeability");
   relative_permeability.push_alias("mur");
   relative_permeability.set_loaded(false);
   relative_permeability.set_positive_required(true);
   relative_permeability.set_non_negative_required(false);
   relative_permeability.set_lowerLimit(1);
   relative_permeability.set_upperLimit(1e6);
   relative_permeability.set_checkLimits(checkLimits_);

   // loss

   loss.push_alias("loss_tangent");
   loss.push_alias("tand");
   loss.push_alias("tandel");
   loss.push_alias("conductivity");
   loss.push_alias("sigma");
   loss.set_loaded(false);
   loss.set_positive_required(false);
   loss.set_non_negative_required(true);
   loss.set_lowerLimit(0);
   loss.set_upperLimit(1e8);  // much too high for loss tangent, unavoidable
   loss.set_checkLimits(checkLimits_);

   // Rz

   Rz.push_alias("Rz");
   Rz.set_loaded(false);
   Rz.set_positive_required(false);
   Rz.set_non_negative_required(true);
   Rz.set_lowerLimit(0);
   Rz.set_upperLimit(0.0001);
   Rz.set_checkLimits(checkLimits_);
}

void Frequency::print(string indent)
{
   PetscPrintf(PETSC_COMM_WORLD,"%d: %s%sFrequency\n",startLine,indent.c_str(),indent.c_str());

   if (frequency.is_loaded()) {
      if (frequency.is_any()) {
         PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s%s=any\n",frequency.get_lineNumber(),indent.c_str(),indent.c_str(),indent.c_str(),frequency.get_keyword().c_str());
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s%s=%g\n",
                                     frequency.get_lineNumber(),indent.c_str(),indent.c_str(),indent.c_str(),frequency.get_keyword().c_str(),frequency.get_dbl_value());
      }
   }

   if (relative_permittivity.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s%s=%g\n",
                  relative_permittivity.get_lineNumber(),indent.c_str(),indent.c_str(),indent.c_str(),relative_permittivity.get_keyword().c_str(),relative_permittivity.get_dbl_value());
   }

   if (relative_permeability.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s%s=%g\n",
                  relative_permeability.get_lineNumber(),indent.c_str(),indent.c_str(),indent.c_str(),relative_permeability.get_keyword().c_str(),relative_permeability.get_dbl_value());
   }

   if (loss.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s%s=%g\n",
                  loss.get_lineNumber(),indent.c_str(),indent.c_str(),indent.c_str(),loss.get_keyword().c_str(),loss.get_dbl_value());
   }

   if (Rz.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s%s=%g\n",
                  Rz.get_lineNumber(),indent.c_str(),indent.c_str(),indent.c_str(),Rz.get_keyword().c_str(),Rz.get_dbl_value());
   }

   PetscPrintf(PETSC_COMM_WORLD,"%d: %s%sEndFrequency\n",endLine,indent.c_str(),indent.c_str());
}

bool Frequency::load (string *indent, inputFile *inputs)
{
   bool fail=false;

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {
      string token,value,line;
      line=inputs->get_line(lineNumber);
      get_token_pair(&line,&token,&value,&lineNumber,indent->c_str());

      int recognized=0;
      if (relative_permittivity.match_alias(&token)) {
         recognized++;
         if (relative_permittivity.loadDouble(&token, &value, lineNumber)) fail=true;
      }

      if (relative_permeability.match_alias(&token)) {
         recognized++;
         if (relative_permeability.loadDouble(&token, &value, lineNumber)) fail=true;
      }

      if (loss.match_alias(&token)) {
         recognized++;
         if (loss.loadDouble(&token, &value, lineNumber)) fail=true;
      }

      if (Rz.match_alias(&token)) {
         recognized++;
         if (Rz.loadDouble(&token, &value, lineNumber)) fail=true;
      }

      if (frequency.match_alias(&token)) {
         if (frequency.is_loaded()) {
            PetscPrintf(PETSC_COMM_WORLD,"%sERROR531: Duplicate entry at line %d for previous entry at line %d.\n",
                                         indent->c_str(),lineNumber,frequency.get_lineNumber());
            fail=true;
         } else {
            if (value.compare("any") == 0) {
               frequency.set_keyword(token);
               frequency.set_value(value);
               frequency.set_lineNumber(lineNumber);
               frequency.set_loaded(true);
            } else {
               if (frequency.loadDouble(&token, &value, lineNumber)) fail=true;
            }
         }
         recognized++;
      }

      // should recognize one keyword
      if (recognized != 1) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR532: Unrecognized keyword at line %d.\n",indent->c_str(),lineNumber);
         fail=true;
      }

      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }
   return fail;
}

bool Frequency::inFrequencyBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Frequency::check (string indent)
{
   bool fail=false;

   if (!frequency.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR533: Frequency block at line %d must specify a frequency.\n",indent.c_str(),startLine);
      fail=true;
   }

   if (!relative_permittivity.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR534: Frequency block at line %d must specify a relative permitivitty.\n",indent.c_str(),startLine);
      fail=true;
   }

   if (!relative_permeability.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR535: Frequency block at line %d must specify a relative permeability.\n",indent.c_str(),startLine);
      fail=true;
   }

   if (!loss.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR536: Frequency block at line %d must specify a loss tangent or a conductivity.\n",indent.c_str(),startLine);
      fail=true;
   }

   if (!Rz.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR536: Frequency block at line %d must specify Rz.\n",indent.c_str(),startLine);
      fail=true;
   }

   return fail;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Temperature
///////////////////////////////////////////////////////////////////////////////////////////

Temperature::Temperature (int startLine_, int endLine_, bool checkLimits_)
{
   startLine=startLine_;
   endLine=endLine_;

   // temperature

   temperature.push_alias("temperature");
   temperature.push_alias("temp");
   temperature.push_alias("t");
   temperature.set_loaded(false);
   temperature.set_positive_required(false);
   temperature.set_non_negative_required(false);
   temperature.set_lowerLimit(-273.15);
   temperature.set_upperLimit(1e5);
   temperature.set_checkLimits(checkLimits_);

   // er_infinity

   er_infinity.push_alias("er_infinity");
   er_infinity.push_alias("epsr_infinity");
   er_infinity.set_loaded(false);
   er_infinity.set_positive_required(true);
   er_infinity.set_non_negative_required(false);
   er_infinity.set_lowerLimit(1);
   er_infinity.set_upperLimit(1e6);
   er_infinity.set_checkLimits(checkLimits_);

   // delta_er

   delta_er.push_alias("delta_er");
   delta_er.push_alias("delta_epsr");
   delta_er.set_loaded(false);
   delta_er.set_positive_required(true);
   delta_er.set_non_negative_required(false);
   delta_er.set_lowerLimit(0);
   delta_er.set_upperLimit(1e6);
   delta_er.set_checkLimits(checkLimits_);

   // m1

   m1.push_alias("m1");
   m1.set_loaded(false);
   m1.set_positive_required(false);
   m1.set_non_negative_required(true);
   m1.set_lowerLimit(0);
   m1.set_upperLimit(100);
   m1.set_checkLimits(checkLimits_);

   // m2

   m2.push_alias("m2");
   m2.set_loaded(false);
   m2.set_positive_required(false);
   m2.set_non_negative_required(true);
   m2.set_lowerLimit(0);
   m2.set_upperLimit(100);
   m2.set_checkLimits(checkLimits_);

   // relative permeability

   relative_permeability.push_alias("relative_permeability");
   relative_permeability.push_alias("mur");
   relative_permeability.set_loaded(false);
   relative_permeability.set_positive_required(true);
   relative_permeability.set_non_negative_required(false);
   relative_permeability.set_lowerLimit(1);
   relative_permeability.set_upperLimit(1e6);
   relative_permeability.set_checkLimits(checkLimits_);

   // loss

   loss.push_alias("loss_tangent");
   loss.push_alias("tand");
   loss.push_alias("tandel");
   loss.push_alias("conductivity");
   loss.push_alias("sigma");
   loss.set_loaded(false);
   loss.set_positive_required(false);
   loss.set_non_negative_required(true);
   loss.set_lowerLimit(0);
   loss.set_upperLimit(1e8);  // much too high for loss tangent, unavoidable
   loss.set_checkLimits(checkLimits_);
}

void Temperature::print(string indent)
{
   PetscPrintf(PETSC_COMM_WORLD,"%d: %sTemperature\n",startLine,indent.c_str());

   if (temperature.is_loaded()) {
      if (temperature.is_any()) {
         PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=any\n",
                     temperature.get_lineNumber(),indent.c_str(),indent.c_str(),temperature.get_keyword().c_str());
      } else {
         PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=%g\n",
                     temperature.get_lineNumber(),indent.c_str(),indent.c_str(),temperature.get_keyword().c_str(),temperature.get_dbl_value());
      }
   }

   if (er_infinity.is_loaded()) {
     PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=%g\n",
                 er_infinity.get_lineNumber(),indent.c_str(),indent.c_str(),er_infinity.get_keyword().c_str(),er_infinity.get_dbl_value());
   }

   if (delta_er.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=%g\n",
                  delta_er.get_lineNumber(),indent.c_str(),indent.c_str(),delta_er.get_keyword().c_str(),delta_er.get_dbl_value());
   }

   if (m1.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=%g\n",
                  m1.get_lineNumber(),indent.c_str(),indent.c_str(),m1.get_keyword().c_str(),m1.get_dbl_value());
   }

   if (m2.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=%g\n",
                  m2.get_lineNumber(),indent.c_str(),indent.c_str(),m2.get_keyword().c_str(),m2.get_dbl_value());
   }

   if (relative_permeability.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=%g\n",
                  relative_permeability.get_lineNumber(),indent.c_str(),indent.c_str(),relative_permeability.get_keyword().c_str(),relative_permeability.get_dbl_value());
   }

   if (loss.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s=%d\n",
                  loss.get_lineNumber(),indent.c_str(),indent.c_str(),loss.get_keyword().c_str(),loss.get_lineNumber());
   }

   long unsigned int i=0;
   while (i < frequencyList.size()) {
      frequencyList[i]->print(indent);
      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"%d: %sEndTemperature\n",endLine,indent.c_str());
}

bool Temperature::findFrequencyBlocks(inputFile *inputs, bool checkLimits)
{
   bool fail=false;
   int start_lineNumber=startLine;
   int stop_lineNumber=endLine;
   int block_start,block_stop;

   while (start_lineNumber < stop_lineNumber) {
      if (inputs->findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                 "Frequency","EndFrequency", false)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Frequency *newFrequency=new Frequency(block_start,block_stop,checkLimits);
            frequencyList.push_back(newFrequency);
         }
      }
      start_lineNumber=inputs->get_next_lineNumber(block_stop);
   }
   return fail;
}

bool Temperature::inFrequencyBlocks(int lineNumber)
{
   long unsigned int i=0;
   while (i < frequencyList.size()) {
      if (frequencyList[i]->inFrequencyBlock(lineNumber)) return true;
      i++;
   }
   return false;
}

bool Temperature::load(string *indent, inputFile *inputs, bool checkInputs)
{
   bool fail=false;

   // Frequency blocks - must find before keywords

   if (findFrequencyBlocks(inputs, checkInputs)) fail=true;

   long unsigned int i=0;
   while (i < frequencyList.size()) {
      if (frequencyList[i]->load(indent,inputs)) fail=true;
      i++;
   }

   // now the keyword pairs

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      if (!inFrequencyBlocks(lineNumber)) {
         string token,value,line; 
         line=inputs->get_line(lineNumber);
         get_token_pair(&line,&token,&value,&lineNumber,*indent);

         int recognized=0;
         if (frequencyList.size() == 0) {

            if (er_infinity.match_alias(&token)) {
               recognized++;
               if (er_infinity.loadDouble(&token, &value, lineNumber)) fail=true;
            }

            if (delta_er.match_alias(&token)) {
               recognized++;
               if (delta_er.loadDouble(&token, &value, lineNumber)) fail=true;
            }

            if (m1.match_alias(&token)) {
               recognized++;
               if (m1.loadDouble(&token, &value, lineNumber)) fail=true;
            }

            if (m2.match_alias(&token)) {
               recognized++;
               if (m2.loadDouble(&token, &value, lineNumber)) fail=true;
            }

            if (relative_permeability.match_alias(&token)) {
               recognized++;
               if (relative_permeability.loadDouble(&token, &value, lineNumber)) fail=true;
            }

            if (loss.match_alias(&token)) {
               recognized++;
               if (loss.loadDouble(&token, &value, lineNumber)) fail=true;
            }

         } else {
            if (er_infinity.match_alias(&token) || delta_er.match_alias(&token) ||
                m1.match_alias(&token) || m2.match_alias(&token) ||
                relative_permeability.match_alias(&token) || loss.match_alias(&token)) {
                PetscPrintf(PETSC_COMM_WORLD,"%sERROR537: Debye variable at line %d is not allowed with frequency blocks defined.\n",
                                             indent->c_str(),lineNumber);
                fail=true;
            }
         }

         if (temperature.match_alias(&token)) {
            if (temperature.is_loaded()) {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR538: Duplicate entry at line %d for previous entry at line %d.\n",
                                            indent->c_str(),lineNumber,temperature.get_lineNumber());
               fail=true;
            } else {
               if (value.compare("any") == 0) {
                  temperature.set_keyword(token);
                  temperature.set_value(value);
                  temperature.set_lineNumber(lineNumber);
                  temperature.set_loaded(true);
               } else {
                  if (temperature.loadDouble(&token, &value, lineNumber)) fail=true;
               }
            }
            recognized++;
         }

         // should recognize one keyword
         if (recognized != 1) {
            PetscPrintf(PETSC_COMM_WORLD,"%sERROR539: Unrecognized keyword at line %d.\n",indent->c_str(),lineNumber);
            fail=true;
         }
      }

      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   return fail;
}

bool Temperature::inTemperatureBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

bool Temperature::check(string indent)
{
   bool fail=false;

   if (!temperature.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR540: Temperature block at line %d must specify a temperature.\n",
                                   indent.c_str(),startLine);
      fail=true;
   }

   if (frequencyList.size() == 0) {
      if (!er_infinity.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR541: Temperature block at line %d must specify an er_infinity.\n",
                                      indent.c_str(),startLine);
         fail=true;
      }

      if (!delta_er.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR542: Temperature block at line %d must specify a delta_er.\n",
                                      indent.c_str(),startLine);
         fail=true;
      }

      if (!m1.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR543: Temperature block at line %d must specify an m1.\n",
                                      indent.c_str(),startLine);
         fail=true;
      }

      if (!m2.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR544: Temperature block at line %d must specify an m2.\n",
                                     indent.c_str(),startLine);
         fail=true;
      }

      if (!relative_permeability.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR545: Temperature block at line %d must specify a relative permeability.\n",
                                      indent.c_str(),startLine);
         fail=true;
      }

      if (!loss.is_loaded()) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR546: Temperature block at line %d must specify a loss tangent or a conductivity.\n",
                                      indent.c_str(),startLine);
         fail=true;
      }

   }

   // frequency block checks
   long unsigned int i=0;
   while (i < frequencyList.size()) {

      // single block checks
      if (frequencyList[i]->check(indent)) fail=true;

      // cross-block checks
      Frequency *test=get_frequency(i);
      long unsigned int j=i+1; 
      while (i < frequencyList.size()-1 && j < frequencyList.size()) {
         Frequency *check=get_frequency(j);

         if (test->get_frequency()->is_any()) {
            if (check->get_frequency()->is_any()) {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR547: Temperature block at line %d incorrectly specifies another frequency block with frequency=any at line %d.\n",
                                            indent.c_str(),startLine,check->get_startLine());
               fail=true;
            } else {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR548: Temperature block at line %d incorrectly specifies a frequency block at line %d after specifying \"any\" at line %d.\n",
                                            indent.c_str(),startLine,check->get_startLine(),test->get_frequency()->get_lineNumber());
               fail=true;
            }
         } else {
            if (check->get_frequency()->is_any()) {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR549: Temperature block at line %d incorrectly specifies a frequency block at line %d after specifying \"any\" at line %d.\n",
                                            indent.c_str(),startLine,test->get_startLine(),check->get_frequency()->get_lineNumber());
               fail=true;
            } else {
               if (test->get_frequency()->dbl_compare(check->get_frequency())) {
                  PetscPrintf(PETSC_COMM_WORLD,"%sERROR550: Temperature block at line %d has frequency blocks with duplicated frequencies at lines %d and %d.\n",
                                               indent.c_str(),startLine,test->get_frequency()->get_lineNumber(),check->get_frequency()->get_lineNumber());
                  fail=true;
               }
            }
         }
         j++;
      }
      i++;
   }

   return fail;
}

complex<double> Temperature::get_eps(double frequency_, double tolerance, string indent)
{
   double eps;
   double loss_;
   complex<double> complex_eps=complex<double>(-DBL_MAX,0);
   double eps0=8.8541878176e-12;

   double freq_low=0;         // for linear interpolation - extrapolation is not supported
   double eps_low;            // ToDo - replace with spline interpolation or
   double loss_low;           //        with a Hilbert Transform algorithm for improved causality
   double freq_high=DBL_MAX;  //
   double eps_high;           // 
   double loss_high;          //
   double freq_test;          //
   bool found_low=false;      //
   bool found_high=false;     //

   if (frequencyList.size() > 0) {

      // frequency list
      bool found=false;
      long unsigned int i=0;
      while (i < frequencyList.size()) {

         // any
         if (frequencyList[i]->get_frequency()->get_value().compare("any") == 0) {
            eps=frequencyList[i]->get_relative_permittivity()->get_dbl_value();
            loss_=frequencyList[i]->get_loss()->get_dbl_value();
            found=true;
            break;
         }

         // exact match
         if (double_compare (frequency_, frequencyList[i]->get_frequency()->get_dbl_value(), tolerance)) {
            eps=frequencyList[i]->get_relative_permittivity()->get_dbl_value();
            loss_=frequencyList[i]->get_loss()->get_dbl_value();
            found=true;
            break;
         }

         // linear interoplation, if needed
         freq_test=frequencyList[i]->get_frequency()->get_dbl_value();

         if (freq_test < frequency_ && freq_test > freq_low) {
            freq_low=freq_test;
            eps_low=frequencyList[i]->get_relative_permittivity()->get_dbl_value();
            loss_low=frequencyList[i]->get_loss()->get_dbl_value();
            found_low=true;
         }

         if (freq_test > frequency_ && freq_test < freq_high) {
            freq_high=freq_test;
            eps_high=frequencyList[i]->get_relative_permittivity()->get_dbl_value();
            loss_high=frequencyList[i]->get_loss()->get_dbl_value();
            found_high=true;
         }

         i++;
      }

      // linear interpolation
      if (! found && found_low && found_high) {
         found=true;
         eps=eps_low+(frequency_-freq_low)/(freq_high-freq_low)*(eps_high-eps_low);
         loss_=loss_low+(frequency_-freq_low)/(freq_high-freq_low)*(loss_high-loss_low);
      }

      // wrap up
      if (found) {
         if (frequencyList[i]->get_loss()->get_keyword().compare("loss_tangent") == 0 ||
             frequencyList[i]->get_loss()->get_keyword().compare("tand") == 0 || 
             frequencyList[i]->get_loss()->get_keyword().compare("tandel") == 0) {
            complex_eps=complex<double>(eps*eps0,-loss_*eps*eps0);            // loss tangent
         } else {
            complex_eps=complex<double>(eps*eps0,-loss_/(2*M_PI*frequency_)); // conductivity
         }
      }

   } else {

      // Debye model

      // do the calculation in conductivity using the sigma variable
      double sigma;
      eps=relative_permeability.get_dbl_value();
      loss_=loss.get_dbl_value();

      if (loss.get_keyword().compare("loss_tangent") == 0 ||
          loss.get_keyword().compare("tand") == 0 ||
          loss.get_keyword().compare("tandel") == 0) {
         sigma=2*M_PI*frequency_*eps*eps0*loss_;   // loss tangent
      } else {
         sigma=loss_;                              // conductivity
      }

      complex<double> t1=complex<double>(pow(10,m1.get_dbl_value()),2*M_PI*frequency_);
      complex<double> t2=complex<double>(pow(10,m2.get_dbl_value()),2*M_PI*frequency_);
      complex<double> denom=complex<double>(log(10),0);
      complex<double> conductivity_term=complex<double>(0,-sigma/(2*M_PI*frequency_));
      complex<double> infinity_term=complex<double>(er_infinity.get_dbl_value()*eps0,0);
      complex<double> delta_term=complex<double>(delta_er.get_dbl_value()*eps0/(m2.get_dbl_value()-m1.get_dbl_value()),0);

      complex_eps=infinity_term+delta_term*log(t2/t1)/denom+conductivity_term;
   }

   return complex_eps;
}

double Temperature::get_mu(double frequency_, double tolerance, string indent)
{
   double mu=-DBL_MAX;

   double freq_low=0;         // for linear interpolation - extrapolation is not supported
   double mu_low;             // ToDo - replace with spline interpolation or
   double freq_high=DBL_MAX;  //        with a Hilbert Transform algorithm for improved causality
   double mu_high;            //
   double freq_test;          //
   bool found_low=false;      //
   bool found_high=false;     //

   if (frequencyList.size() > 0) {

      // frequency list
      bool found=false;
      long unsigned int i=0;
      while (i < frequencyList.size()) {

         // any
         if (frequencyList[i]->get_frequency()->get_value().compare("any") == 0) {
            mu=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            found=true;
            break;
         }

         // exact match
         if (double_compare (frequency_, frequencyList[i]->get_frequency()->get_dbl_value(), tolerance)) {
            mu=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            found=true;
            break;
         }

         // linear interoplation, if needed
         freq_test=frequencyList[i]->get_frequency()->get_dbl_value();

         if (freq_test < frequency_ && freq_test > freq_low) {
            freq_low=freq_test;
            mu_low=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            found_low=true;
         }

         if (freq_test > frequency_ && freq_test < freq_high) {
            freq_high=freq_test;
            mu_high=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            found_high=true;
         }

         i++;
      }

      // linear interpolation
      if (! found && found_low && found_high) {
         found=true;
         mu=mu_low+(frequency_-freq_low)/(freq_high-freq_low)*(mu_high-mu_low);
      }

      // wrap up
      if (found) mu=4e-7*M_PI*mu;

   } else {

      // Debye model
      mu=4e-7*M_PI*relative_permeability.get_dbl_value();
   }

   return mu;
}

double Temperature::get_Rs(double frequency_, double tolerance, string indent)
{
   double Rs=-DBL_MAX;
   double loss_;
   double mur_;
   double Rz_;
   string loss_tangent="loss_tangent";

   double freq_low=0;         // for linear interpolation - extrapolation is not supported
   double loss_low;           // ToDo - replace with spline interpolation or
   double mur_low;            //        with a Hilbert Transform algorithm for improved causality
   double Rz_low;             //
   double freq_high=DBL_MAX;  //
   double loss_high;          //
   double mur_high;           //
   double Rz_high;            //
   double freq_test;          //
   bool found_low=false;      //
   bool found_high=false;     //

   if (frequencyList.size() > 0) {

      // frequency list
      bool found=false;
      long unsigned int i=0;
      while (i < frequencyList.size()) {

         // any
         if (frequencyList[i]->get_frequency()->get_value().compare("any") == 0) {
            loss_=frequencyList[i]->get_loss()->get_dbl_value();
            mur_=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            Rz_=frequencyList[i]->get_Rz()->get_dbl_value();
            found=true;
            break;
         }

         // exact match
         if (double_compare (frequency_, frequencyList[i]->get_frequency()->get_dbl_value(), tolerance)) {
            loss_=frequencyList[i]->get_loss()->get_dbl_value();
            mur_=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            Rz_=frequencyList[i]->get_Rz()->get_dbl_value();
            found=true;
            break;
         }

         // linear interoplation, if needed
         freq_test=frequencyList[i]->get_frequency()->get_dbl_value();

         if (freq_test < frequency_ && freq_test > freq_low) {
            freq_low=freq_test;
            loss_low=frequencyList[i]->get_loss()->get_dbl_value();
            mur_low=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            Rz_low=frequencyList[i]->get_Rz()->get_dbl_value();
            found_low=true;
         }

         if (freq_test > frequency_ && freq_test < freq_high) {
            freq_high=freq_test;
            loss_high=frequencyList[i]->get_loss()->get_dbl_value();
            mur_high=frequencyList[i]->get_relative_permeability()->get_dbl_value();
            Rz_high=frequencyList[i]->get_Rz()->get_dbl_value();
            found_high=true;
         }

         i++;
      }

      // linear interpolation
      if (! found && found_low && found_high) {
         found=true;
         loss_=loss_low+(frequency_-freq_low)/(freq_high-freq_low)*(loss_high-loss_low);
         mur_=mur_low+(frequency_-freq_low)/(freq_high-freq_low)*(mur_high-mur_low);
         Rz_=Rz_low+(frequency_-freq_low)/(freq_high-freq_low)*(Rz_high-Rz_low);
      }

      // wrap up
      if (found) {
         if (frequencyList[i]->get_loss()->get_keyword().compare("loss_tangent") == 0 ||
             frequencyList[i]->get_loss()->get_keyword().compare("tand") == 0 ||
             frequencyList[i]->get_loss()->get_keyword().compare("tandel") == 0) {
            PetscPrintf(PETSC_COMM_WORLD,"%s%s%sERROR551: Attempt to use a dielectric model for an Rs calculation.\n",indent.c_str(),indent.c_str(),indent.c_str());
         } else {

            Rs=sqrt(M_PI*frequency_*4e-7*M_PI*mur_/loss_);

            // apply a correction factor for surface roughness using equation (5) from:
            //     Vladimir Dmitriev-Zdorov and Lambert Simonovich, "Causal Version of Conductor Roughness Models and
            //     its Effect on Characteristics of Transmission Lines", 2017 IEEE 26th Conference on Electrical Performance
            //     of Electronic Packaging and Systems (EPEPS), 2017.

            double r=Rz_/(2*sqrt(3)*2*(1+sqrt(2)));
            double Aflat=36*r*r;
            double delta=sqrt(1/(M_PI*frequency_*4e-7*M_PI*mur_*loss_));
            double factor;
            if (Rz_ == 0) factor=1;
            else factor=1+84*(M_PI*r*r/Aflat)/(1+delta/r+delta*delta/(2*r*r));

            Rs*=factor;
         }
      }
   } else {

      // Debye model
      PetscPrintf(PETSC_COMM_WORLD,"%s%s%sERROR552: Attempt to use a Debye model for an Rs calculation.\n",indent.c_str(),indent.c_str(),indent.c_str());
   }

   return Rs;
}

Temperature::~Temperature()
{
   long unsigned int i=0;
   while (i < frequencyList.size()) {
      delete frequencyList[i];
      i++;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////////////////

Source::Source (int startLine_, int endLine_)
{
   startLine=startLine_;
   endLine=endLine_;
}

void Source::print(string indent)
{
   PetscPrintf(PETSC_COMM_WORLD,"%d: %sSource\n",startLine,indent.c_str());

   long unsigned int i=0;
   while (i < lineNumberList.size()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s%s\n",lineNumberList[i],indent.c_str(),indent.c_str(),lineList[i].c_str());
      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"%d: %sEndSource\n",endLine,indent.c_str());
}

bool Source::load(inputFile *inputs)
{
   int start_lineNumber=inputs->get_next_lineNumber(startLine);
   int stop_lineNumber=inputs->get_previous_lineNumber(endLine);

   while (start_lineNumber <= stop_lineNumber) {
      lineNumberList.push_back(start_lineNumber);
      lineList.push_back(inputs->get_line(start_lineNumber));
      start_lineNumber=inputs->get_next_lineNumber(start_lineNumber);
   }
   return false;
}

bool Source::inSourceBlock (int lineNumber)
{
   if (lineNumber >= startLine && lineNumber <= endLine) return true;
   return false;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Material
///////////////////////////////////////////////////////////////////////////////////////////


Material::Material (int startLine_, int endLine_)
{
   startLine=startLine_;
   endLine=endLine_;

   name.push_alias("name");
   name.set_loaded(false);
   name.set_positive_required(false);
   name.set_non_negative_required(false);
   name.set_lowerLimit(0);
   name.set_upperLimit(0);
   name.set_checkLimits(false);
}

void Material::print(string indent)
{
   PetscPrintf(PETSC_COMM_WORLD,"%d: Material",startLine);
   if (merged) PetscPrintf(PETSC_COMM_WORLD," (merged)");
   PetscPrintf(PETSC_COMM_WORLD,"\n");

   if (name.is_loaded()) {
      PetscPrintf(PETSC_COMM_WORLD,"%d: %s%s\n",name.get_lineNumber(),indent.c_str(),name.get_value().c_str());
   }

   long unsigned int i=0;
   while (i < temperatureList.size()) {
      temperatureList[i]->print(indent);
      i++;
   }

   i=0;
   while (i < sourceList.size()) {
      sourceList[i]->print(indent);
      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"%d: EndMaterial\n",endLine);
}

bool Material::findTemperatureBlocks(inputFile *inputs, bool checkLimits)
{
   bool fail=false;
   int start_lineNumber=startLine;
   int stop_lineNumber=endLine;
   int block_start,block_stop;

   while (start_lineNumber < stop_lineNumber) {
      if (inputs->findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                 "Temperature", "EndTemperature", false)) {
         fail=true;
      } else {
         if (block_start >=0 && block_stop >= 0) {
            Temperature *newTemperature=new Temperature(block_start,block_stop,checkLimits);
            temperatureList.push_back(newTemperature);
         }
      }
      start_lineNumber=inputs->get_next_lineNumber(block_stop);
   }
   return fail;
}

bool Material::findSourceBlocks(inputFile *inputs)
{
   bool fail=false;
   int start_lineNumber=startLine;
   int stop_lineNumber=endLine;
   int block_start,block_stop;

   while (start_lineNumber < stop_lineNumber) {
      if (inputs->findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                            "Source", "EndSource", false)) {
         fail=true;
      } else {
         if (block_start >=0 && block_stop >= 0) {
            Source *newSource=new Source(block_start,block_stop);
            sourceList.push_back(newSource);
         }
      }
      start_lineNumber=inputs->get_next_lineNumber(block_stop);
   }
   return fail;
}

bool Material::inTemperatureBlocks(int lineNumber)
{
   long unsigned int i=0;
   while (i < temperatureList.size()) {
      if (temperatureList[i]->inTemperatureBlock(lineNumber)) return true;
      i++;
   }
   return false;
}

bool Material::inSourceBlocks(int lineNumber)
{
   long unsigned int i=0;
   while (i < sourceList.size()) {
      if (sourceList[i]->inSourceBlock(lineNumber)) return true;
      i++;
   }
   return false;
}

bool Material::load(string *indent, inputFile *inputs, bool checkInputs)
{
   bool fail=false;

   // Temperature blocks
   if (findTemperatureBlocks(inputs, checkInputs)) fail=true;

   long unsigned int i=0;
   while (i < temperatureList.size()) {
      if (temperatureList[i]->load(indent,inputs,checkInputs)) fail=true;
      i++;
   }

   // Souce blocks
   if (findSourceBlocks(inputs)) fail=true;

   i=0;
   while (i < sourceList.size()) {
      if (sourceList[i]->load(inputs)) fail=true;
      i++;
   }

   //  now the keyword pairs

   int lineNumber=inputs->get_next_lineNumber(startLine);
   int stopLineNumber=inputs->get_previous_lineNumber(endLine);
   while (lineNumber <= stopLineNumber) {

      if (!inTemperatureBlocks(lineNumber) && !inSourceBlocks(lineNumber)) {
         string token,value,line;
         line=inputs->get_line(lineNumber);
         get_token_pair(&line,&token,&value,&lineNumber,*indent);

         int recognized=0;
         if (name.match_alias(&token)) {
            if (name.is_loaded()) {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR553: Duplicate entry at line %d for previous entry at line %d.\n",
                                            indent->c_str(),lineNumber,name.get_lineNumber());
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
            PetscPrintf(PETSC_COMM_WORLD,"%sERROR554: Unrecognized keyword at line %d.\n",indent->c_str(),lineNumber);
            fail=true;
         }
      }
      lineNumber=inputs->get_next_lineNumber(lineNumber);
   }

   return fail;
}

bool Material::check(string indent)
{
   bool fail=false;

   if (name.is_loaded()) {
      if (name.get_value().compare("") == 0) {
         PetscPrintf(PETSC_COMM_WORLD,"%sERROR555: name is blank at line %d.\n",indent.c_str(),name.get_lineNumber());
         fail=true;
      }
   } else {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR556: Material block at line %d must specify a name.\n",indent.c_str(),startLine);
      fail=true;
   }

   // no Temperature blocks
   if (temperatureList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR557: Material block at line %d must specify at least one Temperature block.\n",indent.c_str(),startLine);
      fail=true;
   }

   // no Source blocks
   if (sourceList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR558: Material block at line %d must specify at least one Source block.\n",indent.c_str(),startLine);
      fail=true;
   }

   // temperature block checks
   long unsigned int i=0;
   while (i < temperatureList.size()) {

      // single block checks
      if (temperatureList[i]->check(indent)) fail=true;

      // cross-block checks
      long unsigned int j=i+1;
      while (i < temperatureList.size()-1 && j < temperatureList.size()) {
         if (temperatureList[i]->get_temperature()->is_any()) {
            if (temperatureList[j]->get_temperature()->is_any()) {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR559: Material block at line %d incorrectly specifies another temperature block with temperature=any at line %d\n",
                                            indent.c_str(),startLine,temperatureList[j]->get_startLine());
               fail=true;
            } else {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR560: Material block at line %d incorrectly specifies a temperature block at line %d after specifying \"any\" at line %d.\n",
                                            indent.c_str(),startLine,temperatureList[j]->get_startLine(),temperatureList[i]->get_temperature()->get_lineNumber());
               fail=true;
            }
         } else {
            if (temperatureList[j]->get_temperature()->is_any()) {
               PetscPrintf(PETSC_COMM_WORLD,"%sERROR561: Material block at line %d incorrectly specifies a temperature block at line %d after specifying \"any\" at line %d.\n",
                                            indent.c_str(),startLine,temperatureList[i]->get_startLine(),temperatureList[j]->get_temperature()->get_lineNumber());
               fail=true;
            } else {
               if (temperatureList[i]->get_temperature()->dbl_compare(temperatureList[j]->get_temperature())) {
                  PetscPrintf(PETSC_COMM_WORLD,"%sERROR562: Material block at line %d has temperature blocks with duplicated temperatures at lines %d and %d.\n",
                                               indent.c_str(),startLine,temperatureList[i]->get_temperature()->get_lineNumber(),temperatureList[j]->get_temperature()->get_lineNumber());
                  fail=true;
               }
            }
         }
         j++;
      }
      i++;
   }

   return fail;
}

// Find the temperature block matching the given temperature.
// Do not interpolate or extrapolate.
Temperature* Material::get_temperature(double temperature_, double tolerance, string indent)
{
   long unsigned int i=0;
   while (i < temperatureList.size()) {

      // any
      if (temperatureList[i]->get_temperature()->get_value().compare("any") == 0) {
         return temperatureList[i];
      }

      // exact match
      if (double_compare (temperature_, temperatureList[i]->get_temperature()->get_dbl_value(), tolerance)) {
         return temperatureList[i];
      }

      i++;
   }

   PetscPrintf(PETSC_COMM_WORLD,"%s%s%sERROR563: Failed to find a temperature entry for material \"%s\" at a temperature of %g.\n",
                           indent.c_str(),indent.c_str(),indent.c_str(),name.get_value().c_str(),temperature_);

   return nullptr;
}

complex<double> Material::get_eps(double temperature_, double frequency, double tolerance, string indent)
{
   Temperature *temperature;
   complex<double> complex_eps=complex<double>(-DBL_MAX,0);

   temperature=get_temperature(temperature_, tolerance, indent);
   if (temperature) complex_eps=temperature->get_eps(frequency, tolerance, indent);

   return complex_eps;
}

double Material::get_mu(double temperature_, double frequency, double tolerance, string indent)
{
   Temperature *temperature;
   double mu=-DBL_MAX;

   temperature=get_temperature(temperature_, tolerance, indent);
   if (temperature) mu=temperature->get_mu(frequency, tolerance, indent);

   return mu;
}

double Material::get_Rs(double temperature_, double frequency, double tolerance, string indent)
{
   Temperature *temperature;
   double Rs=-DBL_MAX;

   temperature=get_temperature(temperature_, tolerance, indent);
   if (temperature) Rs=temperature->get_Rs(frequency, tolerance, indent);

   return Rs;
}

Material::~Material ()
{
   long unsigned int i=0;
   while (i < temperatureList.size()) {
      delete temperatureList[i];
      i++;
   }

   i=0;
   while (i < sourceList.size()) {
      delete sourceList[i];
      i++;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////
// MaterialDatabase
///////////////////////////////////////////////////////////////////////////////////////////

bool MaterialDatabase::findMaterialBlocks()
{
   bool fail=false;
   int start_lineNumber=inputs.get_first_lineNumber();
   int stop_lineNumber=inputs.get_last_lineNumber();
   int block_start,block_stop;

   // skip the version line, which must be the first line
   start_lineNumber=inputs.get_next_lineNumber(start_lineNumber);

   while (start_lineNumber < stop_lineNumber) {
      if (inputs.findBlock(start_lineNumber,stop_lineNumber, &block_start, &block_stop,
                                "Material", "EndMaterial", true)) {
         fail=true;
      } else {
         if (block_start >= 0 && block_stop >= 0) {
            Material *newMaterial=new Material(block_start,block_stop);
            materialList.push_back(newMaterial);
         }
      }
      start_lineNumber=inputs.get_next_lineNumber(block_stop);
   }
   return fail;
}

// return true on fail
bool MaterialDatabase::load(const char *path, const char *filename, bool checkInputs)
{
   // assemble the full path name
   char *fullPathName=(char *)malloc((strlen(path)+strlen(filename)+1)*sizeof(char));
   if (!fullPathName) return 1;
   sprintf (fullPathName,"%s%s",path,filename);
   PetscPrintf(PETSC_COMM_WORLD,"%sloading materials file \"%s\"\n",indent.c_str(),fullPathName);
 
   bool fail=false;
   if (inputs.load(fullPathName)) {if (fullPathName) free(fullPathName); fullPathName=nullptr; return true;}
   if (fullPathName) {free(fullPathName); fullPathName=nullptr;}
   inputs.createCrossReference();

   if (inputs.checkVersion(version_name, version_value)) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR564: Version mismatch.  Expecting the first line to be: %s %s\n",
                                   indent.c_str(),indent.c_str(),version_name.c_str(),version_value.c_str());
      return true;
   }

   // load the materials file

   if (findMaterialBlocks()) fail=true;

   long unsigned int i=0;
   while (i < materialList.size()) {
      if (materialList[i]->load(&indent, &inputs, checkInputs)) fail=true;
      i++;
   }

   if (check()) fail=true;

   if (fail) {
      PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR565: Failed to load materials.\n",indent.c_str(),indent.c_str());
      return fail;
   }

   return fail;
};

Material* MaterialDatabase::get(string name)
{
   long unsigned int i=0;
   while (i < materialList.size()) {
      if (materialList[i]->get_name()->get_value().compare(name) == 0) return materialList[i];
      i++;
   }
   return nullptr;
}

void MaterialDatabase::print(string indent)
{
   long unsigned int i=0;
   while (i < materialList.size()) {
      materialList[i]->print(indent);
      i++;
   }
}

bool MaterialDatabase::check()
{
   bool fail=false;

   if (materialList.size() == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR566: No materials loaded.\n",indent.c_str());
      fail=true;
   }

   long unsigned int i=0;
   while (i < materialList.size()) {

      // single block checks
      if (materialList[i]->check(indent)) fail=true;

      // cross-block checks
      long unsigned int j=i+1;
      while (i < materialList.size()-1 && j < materialList.size()) {
         if (materialList[i]->get_name()->get_value().compare(materialList[j]->get_name()->get_value()) == 0) {
            PetscPrintf(PETSC_COMM_WORLD,"%sERROR567: name at line %d duplicates the name at line %d.\n",
                                         indent.c_str(),materialList[j]->get_name()->get_lineNumber(),materialList[i]->get_name()->get_lineNumber());
            fail=true;
         }
         j++;
      }
      i++;
   }

   return fail;
}

// MaterialDatabase takes ownership of the contents of db
bool MaterialDatabase::merge(MaterialDatabase *db, string indent)
{
   // checks for nothing to do
   if (!db) return false;
   if (db->materialList.size() == 0) return false;

   // check for version alignment
   if (version_name.compare(db->version_name) != 0 ||
       version_value.compare(db->version_value) != 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%sERROR568: Material database merge attempted with different versions.\n",indent.c_str());
      return true;
   }

   // merge in db

   if (! double_compare(tol,db->tol,db->tol)) {
      PetscPrintf(PETSC_COMM_WORLD,"%sReplacing tol with the merged database value.\n",indent.c_str());
   }

   // add the new materials while checking for duplicates
   long unsigned int i=0;
   while (i < db->materialList.size()) {
      bool found=false;
      long unsigned int j=0;
      while (j < materialList.size()) {
         if (! materialList[j]->get_merged()) {
            if (db->materialList[i]->get_name()->get_value().compare(materialList[j]->get_name()->get_value()) == 0) {
               PetscPrintf(PETSC_COMM_WORLD,"%sReplacing duplicate material \"%s\" with the material from the local material database.\n",
                                            indent.c_str(),db->materialList[i]->get_name()->get_value().c_str());
               delete materialList[j];
               materialList[j]=db->materialList[i];
               materialList[j]->set_merged(true);
               found=true;
               break;
            }
         }
         j++;
      }

      if (!found) {
         push(db->materialList[i]);
         materialList[materialList.size()-1]->set_merged(true);   // data ownership stays with db
      }
      i++;
   }

   db->isTransferred=true;
   return 0;
}

MaterialDatabase::~MaterialDatabase()
{
   if (isTransferred) return;

   long unsigned int i=0;
   while (i < materialList.size()) {
      delete materialList[i];
      i++;
   }
}

// last used is 568


