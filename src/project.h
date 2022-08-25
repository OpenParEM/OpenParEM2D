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

// Structure for maintaining and using project setup information in OpenParEM tools.
// This is a structure rather than a class so that it can be passed to eigensolve.c.

#ifndef PROJECT_H
#define PROJECT_H

#include "petscsys.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// for the refinement variable
#define NONE 1
#define ALPHA 2
#define BETA 3
#define ABSGAMMA 4
#define ABSZO 5
#define REZO 6
#define IMZO 7

struct inputFrequencyPlan {
   int type;                 // 0 => linear, 1 => log, 2 => point
   double frequency;         // for point
   double start;             // for linear and log
   double stop;              // for linear and log
   double step;              // for linear
   int pointsPerDecade;      // for log
   int refine;               // 0 => do not refine the mesh at the frequency point(s), 2 => refine the mesh
   int lineNumber;
};

struct projectData {
   char* version_name;       // set in init_project 
   char* version_value;      // set in init_project

   // user settable

   int project_save_fields;

   char *mesh_file;
   int mesh_order;
   int mesh_uniform_refinement_count;
   double mesh_refinement_cutoff;
   double mesh_refinement_fraction;

   char *mode_definition_file;

   char *materials_global_path;
   char *materials_global_name;
   char *materials_local_path;
   char *materials_local_name;
   int materials_check_limits;

   char *refinement_frequency;
   char *refinement_variable;
   int refinement_iteration_min;
   int refinement_iteration_max;
   int refinement_required_passes;
   char* refinement_tolerance;

   struct inputFrequencyPlan *inputFrequencyPlans;
   unsigned long int inputFrequencyPlansAllocated;
   unsigned long int inputFrequencyPlansCount;

   int solution_modes;
   double solution_temperature;
   double solution_tolerance;
   int solution_iteration_limit;
   int solution_modes_buffer;
   char *solution_impedance_definition;
   char *solution_impedance_calculation;
   int solution_check_closed_loop;
   int solution_accurate_residual;
   int solution_shift_invert;
   int solution_use_initial_guess;
   double solution_shift_factor;

   int output_show_refining_mesh;
   int output_show_postprocessing;
   int output_show_iterations;
   int output_show_license;

   int test_create_cases;
   int test_show_audit;
   int test_show_detailed_cases;

   int debug_show_memory;
   int debug_show_project;
   int debug_show_frequency_plan;
   int debug_show_materials;
   int debug_show_mode_definitions;
   int debug_show_impedance_details;
   int debug_skip_solve;
   int debug_tempfiles_keep;

   int field_points_count;
   int field_points_allocated;
   double *field_points_x;
   double *field_points_y;

   // not user settable

   char *project_name;
   int refinement_refine_converged_modes; // ToDo - needs work.  Produces knots of highly refine mesh in WR90.
   int solution_active_mode_count;
};

char* allocCopyString (char *);
char* removeNewLineChar (char *);
char *removeComment (char *);
int removeQuote (char *);
void removeQuotes (char *);

#endif
