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

#ifndef FREQUENCYPLAN_H
#define FREQUENCYPLAN_H

#include "mfem.hpp"
#include "petscsys.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "project.h"

using namespace std;
using namespace mfem;


// meshRefinementType:
// 0 - refine at each frequency individually
// 1 - refine at a given refinementFrequency [must be > 0]
// 2 - do not refine at any frequencies

class FrequencyPlanPoint {
   private:
      double frequency;
      int refinementPriority;                // 0 means no refinement
      bool restart;                          // true: restart refinement with initial mesh
      bool simulated;                        // true: used for simulation
      bool active;
   public:
      void set_frequency (double frequency_) {frequency=frequency_;}
      void set_refinementPriority (int refinementPriority_) {refinementPriority=refinementPriority_;}
      void set_restart (bool restart_) {restart=restart_;}
      void set_simulated (bool simulated_) {simulated=simulated_;}
      void set_active (bool active_) {active=active_;}
      double get_frequency () {return frequency;}
      int get_refinementPriority () {return refinementPriority;}
      bool get_restart () {return restart;}
      bool get_simulated () {return simulated;}
      bool get_active () {return active;}
      void print();
};

class FrequencyPlan {
   private:
      vector<FrequencyPlanPoint *> planList;
      int refinedCount;
      bool hasRefined;
      long unsigned int lastRefinedIndex;
   public:
      ~FrequencyPlan();
      bool assemble(struct projectData *);
      bool get_frequency(char *, double *, bool *, bool *);
      bool is_refining();
      void sort();
      void eliminateDuplicates();
      void setAllRefineRestart();
      void setLowRefinementPriority (int);
      void setHighRefinementPriority (int);
      void print();
};


#endif
