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

#include "frequencyPlan.hpp"

void FrequencyPlanPoint::print()
{
   if (active) {
      PetscPrintf(PETSC_COMM_WORLD,"   %10g",frequency);
      if (refinementPriority > 0) {
         PetscPrintf(PETSC_COMM_WORLD,"   refine    %5d",refinementPriority);
         if (restart) PetscPrintf(PETSC_COMM_WORLD,"     restart");
      }
      PetscPrintf(PETSC_COMM_WORLD,"\n");
   }
}

// sort from lowest to highest frequency
void FrequencyPlan::sort()
{
   bool found=true;
   while (found) {
      found=false;
      unsigned long int i=0;
      while (i < planList.size()-1) {
         unsigned long int j=i+1;
         while (planList[i]->get_active() && j < planList.size()) {
            if (planList[j]->get_active() && planList[i]->get_frequency() > planList[j]->get_frequency()) {
               found=true;
               FrequencyPlanPoint *temp;
               temp=planList[j];
               planList[j]=planList[i];
               planList[i]=temp;
            }
            j++;
         }
         i++;
      }
   }
}

// eliminate duplicates
void FrequencyPlan::eliminateDuplicates()
{
   long unsigned int i=0;
   while (i < planList.size()-1) {
      unsigned long int j=i+1;
      while (planList[i]->get_active() && j < planList.size()) {
         if (planList[j]->get_active() && abs(planList[i]->get_frequency()-planList[j]->get_frequency())/planList[j]->get_frequency() < 1e-14) {
            if (planList[j]->get_refinementPriority() > 0) {
              if (planList[i]->get_refinementPriority() == 0) {
                 planList[i]->set_active(false);
                 break;
              } else {
                 if (planList[j]->get_refinementPriority() < planList[i]->get_refinementPriority()) {
                    planList[i]->set_active(false);
                    break;
                 }
              }
            }
            planList[j]->set_active(false);
         }
         j++;
      }
      i++;
   }
}

void FrequencyPlan::setAllRefineRestart()
{
   long unsigned int i=0;
   while (i < planList.size()) {
      if (planList[i]->get_active()) {
         planList[i]->set_refinementPriority((int)i+1);
         planList[i]->set_restart(true);
      }
      i++;
   }
}

void FrequencyPlan::setLowRefinementPriority (int priority)
{
   long unsigned int i=0;
   while (i < planList.size()) {
      if (planList[i]->get_active()) {
         planList[i]->set_refinementPriority(priority);
         if (priority == 1) planList[i]->set_restart(true);
         else planList[i]->set_restart(false);
         break;
      }
      i++;
   }
}

void FrequencyPlan::setHighRefinementPriority (int priority)
{
   long unsigned int i=planList.size()-1;
   while (i >= 0) {
      if (planList[i]->get_active()) {
         planList[i]->set_refinementPriority(priority);
         if (priority == 1) planList[i]->set_restart(true);
         else planList[i]->set_restart(false);
         break;
      }
      if (i == 0) break;
      i--;
   }
}

bool FrequencyPlan::assemble(struct projectData *projData)
{
   unsigned long int LIMIT=10000000;
   refinedCount=0;
   hasRefined=false;
   lastRefinedIndex=INT_MAX;
   int refinementPriority=1;

   // add in the frequencies from the frequency plans in projData
   unsigned long int i=0;
   while (i < projData->inputFrequencyPlansCount) {
      if (projData->inputFrequencyPlans[i].type == 0) {
         double frequency=projData->inputFrequencyPlans[i].start;
         while (frequency <= projData->inputFrequencyPlans[i].stop*(1+1e-12)) {

            FrequencyPlanPoint *planPoint=new FrequencyPlanPoint;
            planList.push_back(planPoint);

            if (planList.size() > LIMIT) {
               PetscPrintf(PETSC_COMM_WORLD,"ERROR318: Excessive frequency count > %ld.\n",LIMIT);
               return true;
            }

            planPoint->set_frequency(frequency);

            if (strcmp(projData->refinement_frequency,"plan") == 0 && 
                projData->inputFrequencyPlans[i].refine == 1) {
               planPoint->set_refinementPriority(refinementPriority);
               if (refinementPriority == 1) planPoint->set_restart(true);
               else planPoint->set_restart(false);
               refinementPriority++;
            } else {
               planPoint->set_refinementPriority(0);
               planPoint->set_restart(true);
            }

            planPoint->set_simulated(false);
            planPoint->set_active(true);

            frequency+=projData->inputFrequencyPlans[i].step;
         }
      } else if (projData->inputFrequencyPlans[i].type == 1) {
         double factor=pow(10,1.0/(double)projData->inputFrequencyPlans[i].pointsPerDecade);
         double frequency=projData->inputFrequencyPlans[i].start;
         while (frequency <= projData->inputFrequencyPlans[i].stop*(1+1e-12)) {

            FrequencyPlanPoint *planPoint=new FrequencyPlanPoint;
            planList.push_back(planPoint);

            if (planList.size() > LIMIT) {
               PetscPrintf(PETSC_COMM_WORLD,"ERROR319: Excessive frequency count > %ld.\n",LIMIT);
               return true;
            }

            planPoint->set_frequency(frequency);

            if (strcmp(projData->refinement_frequency,"plan") == 0 &&
                projData->inputFrequencyPlans[i].refine == 1) {
               planPoint->set_refinementPriority(refinementPriority);
               if (refinementPriority == 1) planPoint->set_restart(true);
               else planPoint->set_restart(false);
               refinementPriority++;
            } else {
               planPoint->set_refinementPriority(0);
               planPoint->set_restart(true);
            }

            planPoint->set_simulated(false);
            planPoint->set_active(true);

            frequency*=factor;
         }
      } else if (projData->inputFrequencyPlans[i].type == 2) {

            FrequencyPlanPoint *planPoint=new FrequencyPlanPoint;
            planList.push_back(planPoint);

            if (planList.size() > LIMIT) {
               PetscPrintf(PETSC_COMM_WORLD,"ERROR320: Excessive frequency count > %ld.\n",LIMIT);
               return true;
            }

            planPoint->set_frequency(projData->inputFrequencyPlans[i].frequency);

            if (strcmp(projData->refinement_frequency,"plan") == 0 &&
                projData->inputFrequencyPlans[i].refine == 1) {
               planPoint->set_refinementPriority(refinementPriority);
               if (refinementPriority == 1) planPoint->set_restart(true);
               else planPoint->set_restart(false);
               refinementPriority++;
            } else {
               planPoint->set_refinementPriority(0);
               planPoint->set_restart(true);
            }

            planPoint->set_simulated(false);
            planPoint->set_active(true);
      }

      i++;
   }

   eliminateDuplicates();
   sort();

   // set remaining cases ("plan" and "none" are taken care of above)

   if (strcmp(projData->refinement_frequency,"all") == 0) setAllRefineRestart();

   if (strcmp(projData->refinement_frequency,"low") == 0) 
      setLowRefinementPriority(refinementPriority++);

   if (strcmp(projData->refinement_frequency,"lowhigh") == 0) {
      setLowRefinementPriority(refinementPriority++);
      setHighRefinementPriority(refinementPriority++);
   }

   if (strcmp(projData->refinement_frequency,"high") == 0) 
      setHighRefinementPriority(refinementPriority++);

   if (strcmp(projData->refinement_frequency,"highlow") == 0) {
      setHighRefinementPriority(refinementPriority++);
      setLowRefinementPriority(refinementPriority++);
   }

   return false;
}

bool FrequencyPlan::is_refining ()
{
   long unsigned int i=0;
   while (i < planList.size()) {
      if (planList[i]->get_active() && planList[i]->get_refinementPriority() > 0) return true;
      i++;
   }
   return false;
}

bool FrequencyPlan::get_frequency (char *refinement_frequency, double *frequency, bool *refine, bool *restart) {

   // get the next refinement case, if any

   if (! hasRefined) {

      // get the next frequency that is refining with the lowest priority

      long unsigned int refinementIndex=0;
      int refinementPriority=INT_MAX;
      bool foundRefinement=false;
      long unsigned int i=0;
      while (i < planList.size()) {
         if (planList[i]->get_active() && !planList[i]->get_simulated() && planList[i]->get_refinementPriority() > 0) {
            foundRefinement=true;
            if (refinementPriority == INT_MAX) {
               refinementPriority=planList[i]->get_refinementPriority();
               refinementIndex=i;
            } else {
               if (planList[i]->get_refinementPriority() < refinementPriority) {
                  refinementPriority=planList[i]->get_refinementPriority();
                  refinementIndex=i;
               }
            }
         }
         i++;
      }

      if (foundRefinement) {
         planList[refinementIndex]->set_simulated(true);
         *frequency=planList[refinementIndex]->get_frequency();
         *refine=true;
         *restart=planList[refinementIndex]->get_restart();
         refinedCount++;
         lastRefinedIndex=refinementIndex;
         return true;
      }
   }

   // finished for case all
   if (strcmp(refinement_frequency,"all") == 0) return false;

   // refinement is over

   // re-run frequencies if there is more than one mesh refinement
   if (!hasRefined && refinedCount > 1) {
      long unsigned int i=0;
      while (i < planList.size()) {
         if (planList[i]->get_active()) planList[i]->set_simulated(false);
         i++;
      }
      planList[lastRefinedIndex]->set_simulated(true);   // do not need to re-run this frequency
   }

   // always true for case "none"
   hasRefined=true;

   // get the next available frequency
   long unsigned int i=0;
   while (i < planList.size()) {
      if (planList[i]->get_active() && !planList[i]->get_simulated()) {
         *frequency=planList[i]->get_frequency();
         *refine=false;
         *restart=false;
         planList[i]->set_simulated(true);
         return true;
      }
      i++;
   }

   return false;
}

void FrequencyPlan::print(){
   PetscPrintf(PETSC_COMM_WORLD,"Frequency Plan:\n");
   PetscPrintf(PETSC_COMM_WORLD,"   ----------------------------------------\n");
   PetscPrintf(PETSC_COMM_WORLD,"    frequency   refinement priority restart\n");
   PetscPrintf(PETSC_COMM_WORLD,"   ----------------------------------------\n");
   long unsigned int i=0;
   while (i < planList.size()) {
      planList[i]->print();
      i++;
   }
   PetscPrintf(PETSC_COMM_WORLD,"   ----------------------------------------\n");
}

FrequencyPlan::~FrequencyPlan()
{
   long unsigned int i=0;
   while (i < planList.size()) {
      delete planList[i];
      i++;
   }
}

// last ERROR 320
// last ERROR - shared with fieldPoints.cpp and frequencyPlan.cpp

