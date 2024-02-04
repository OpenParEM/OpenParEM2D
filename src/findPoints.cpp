/******************************************************************************************

BSD 3-Clause License

Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

******************************************************************************************/

#include "mfem.hpp"

using namespace mfem;

// This is a simplified version of the MFEM function Mesh::FindPoints to get access to options.
// Currently set up for reliability at the expense of run time.  For 2D meshes (dim=2), run time should never be an issue.

void findPoints(ParMesh *pmesh, DenseMatrix &point_mat, Array<int>& elem_ids, Array<IntegrationPoint>& ips, int dim)
{
   int num_mesh_elements=pmesh->GetNE();

   const int npts = point_mat.Width();
   double *data=point_mat.GetData();

   elem_ids.SetSize(npts);
   ips.SetSize(npts);
   elem_ids = -1;

   Vector pt(dim);
   pt.NewDataAndSize(nullptr, dim);

   InverseElementTransformation *inv_tr;
   inv_tr=new InverseElementTransformation;

   // tolerances, defaults shown to the right
   inv_tr->SetMaxIter(64);            // 16     // Set the maximum number of iterations when solving for a reference point.
   inv_tr->SetReferenceTol(1e-15);    // 1e-15  // Set the reference-space convergence tolerance.
   inv_tr->SetPhysicalRelTol(1e-15);  // 1e-15  // Set the relative physical-space convergence tolerance.
   inv_tr->SetElementTol(1e-8);       // 1e-8   // Set the tolerance used to determine if a point lies inside, Newton solver only

   inv_tr->SetPrintLevel(-1);         // -1 - never print
                                      // 0 - print only errors;
                                      // 1 - print the first and last iterations
                                      // 2 - print every iteration
                                      // 3 - print every iteration including point coordinates

   // Highest reliability (slowest) with ClosestRefNode and NewtonElementProject
   // All other options shown

   inv_tr->SetInitialGuessType(InverseElementTransformation::Center);
   //inv_tr->SetInitialGuessType(InverseElementTransformation::ClosestPhysNode);
//   inv_tr->SetInitialGuessType(InverseElementTransformation::ClosestRefNode);   // OpenParEM2D
   //inv_tr->SetInitialGuessType(InverseElementTransformation::GivenPoint); // requires SetInitialGuess, below

   inv_tr->SetSolverType(InverseElementTransformation::Newton);
   //inv_tr->SetSolverType(InverseElementTransformation::NewtonSegmentProject);
//   inv_tr->SetSolverType(InverseElementTransformation::NewtonElementProject);   // OpenParEM2D

   ElementTransformation *eltransf;

   int pts_found=0;
   for (int k = 0; k < npts; k++)
   {
      if (elem_ids[k] < 0) {

         pt.SetData(data+k*dim);

         int element=0;
         while (element < num_mesh_elements) {
            eltransf=pmesh->GetElementTransformation(element);
            inv_tr->SetTransformation(*eltransf);
            //inv_tr->SetInitialGuess(ips[k]);  // required for GivenPoint
            int res = inv_tr->Transform(pt, ips[k]);
            if (res == InverseElementTransformation::Inside)
            {
               elem_ids[k] = element;
               pts_found++;
               break;  // keep the first one found; might be found in multiple elements (i.e. on border)
            }
            element++;
         }
      }
   }

   delete inv_tr;

   return;
}

