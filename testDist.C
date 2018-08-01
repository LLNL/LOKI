/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 ************************************************************************/
#include <stdlib.h>
#include <iostream>
#include "Overture.h"
#include "ParmParse.H"
#include "PerlPreProcessor.H"

using namespace std;
using namespace Loki;

/*
 * Serial utility to find the likely size of the smallest dimension of the
 * domain after Overture decomposes the problem.  The utility issues a warning
 * if this size is less than the minimum necessary for the LOKI algorithms.
 * The decomposition scheme used by Overture could not be found so this is an
 * estimate although it seems to be accurate.  We can't just use Overture
 * because it requires you to be running on the number of processors you're
 * interested in which defeats the purpose of a quick utility to see if that
 * decomposition will work.
 */

// Allocated quantities.  These are global variables to facilitate freeing them
// upon exit.

// The size of the 4 phase space dimensions for each species.
int** dimVec;

// The number of processors decomposing the 4 phase space directions for each
// species.
int** dimProc;

// The total number of decomposing processors for each species.  This is the
// product of the 4 dimProc values for each species.
int* nProc;

// A measure of the work load for each species.
float* work;

// Free allocated quantities.  Made this a function in lieu of a goto statement.
void
cleanup(int number_of_species)
{
   for (int species = 0; species < number_of_species; ++species) {
      delete [] dimVec[species];
      delete [] dimProc[species];
   }
   delete [] dimVec;
   delete [] dimProc;
   delete [] nProc;
   delete [] work;
}

// Returns the smallest prime factor of the supplied integer.
int
leastPrimeFactor(int n)
{
   if (n < 0) {
      cout << "leastPrimeFactor:ERROR: can't factor a negative number" << endl;
      return n;
   }

   if (n < 2) {
      return n;
   }
   if (n % 2 == 0) {
      return 2;
   }

   int sqrtn = int(sqrt(double(n+.5)));
   for (int i = 3; i <= sqrtn; i += 2) {
      if (n % i == 0) {
         return i;
      }
   }
   return n;
}

// Returns the largest prime factor of the supplied integer.
int
greatestPrimeFactor(int n)
{
   if (n < 0) {
      cout << "greatestPrimeFactor:ERROR: can't factor a negative number" << endl;
      return n;
   }

   // After dividing out all the smaller prime factors, what's left is the
   // largest one.
   for (int lpf = leastPrimeFactor(n); lpf != n; lpf = leastPrimeFactor(n)) {
      n = n/lpf;
   }

   return n;
}

// Finds the likely size of the smallest dimension of the domain after Overture
// decomposes the problem.  This utility issues a warning if this size is less
// than the minimum necessary for the LOKI algorithms.  It takes 2 input
// arguments, the total number of Vlasov processors and the input deck for the
// problem of interest.
int
main(int argc, char* argv[])
{
   // Number of phase space dimensions.
   const int nDims = 4;

   // Check that we got the right number of input args.
   if (argc < 3) {
      cout << "Usage: total_number_of_vlasov_processors input_file" << endl;
      cout << "Missing input." << endl;
      return -1;
   }
   if (argc > 3) {
      cout << "Usage: total_number_of_vlasov_processors input_file" << endl;
      cout << "Extra input will be ignored." << endl;
   }

   // Get the total number of Vlasov processors from the input args.
   int nProcTotal = atoi(argv[1]);

   // Read the input file name from the input args and construct the input db.
   PerlPreProcessor ppp(argc, argv);
   string input_file_name(argv[2]);
   printF("\n\n*** Reading input file %s\n\n", input_file_name.c_str());
   string output_file_name(ppp.process(input_file_name));
   char* cstr;
   ParmParse pp(0, &cstr, "", output_file_name.c_str());

   // Read the spatial solution order from the input deck and set nGhosts and
   // stencilWidth based on that value.
   int spatialOrder = 4;
   pp.query("spatial_solution_order", spatialOrder);
   int nGhosts;
   if (spatialOrder == 4) {
      nGhosts = 2;
   }
   else if (spatialOrder == 6) {
      nGhosts = 3;
   }
   else {
      cout << "Spatial solution order must be 4 or 6" << endl;
      return -1;
   }
   int stencilWidth = spatialOrder + 1;

   // Read the number of kinetic species from the input deck.
   int number_of_species = 1;
   pp.query("number_of_species", number_of_species);
   if (number_of_species <= 0) {
      cout << "There must be at least one kinetic species in the problem."
           << endl;
      return -1;
   }

   // For each species allocate an array for the phase space dimension, an array
   // for the partitioning of each phase space dimension, an array for the
   // number of processors allocated to it, and an array for the work load.
   dimVec = new int* [number_of_species];
   dimProc = new int* [number_of_species];
   nProc = new int [number_of_species];
   work = new float[number_of_species];
   for (int species = 0; species < number_of_species; ++species) {
      dimVec[species] = new int [nDims];
      dimProc[species] = new int [nDims];
   }

   // Read the number of spatial points from the input deck.  This is the
   // size of the first 2 dimensions of each species.  Also, start computing
   // species work.
   if (!pp.contains("N")) {
      cout << "The input deck does not contain the number of spatial points."
           << endl;
      cleanup(number_of_species);
      return -1;
   }
   Array<int> tmp(2);
   pp.queryarr("N", tmp, 0, 2);
   for (int species = 0; species < number_of_species; ++species) {
      dimVec[species][0] = tmp[0];
      dimVec[species][1] = tmp[1];
      work[species] = tmp[0] * tmp[1];
   }

   // Read the number of velocity points from the input deck for each species.
   // These are the 3rd and 4th dimensions of each species.  Also, finish
   // computing species work.
   // The net_work is the combined work load of all the species.
   float net_work = 0.0;
   for (int species = 0; species < number_of_species; ++species) {
      char buffer[80];
      sprintf(buffer, "kinetic_species.%i.Nv", species+1);
      if (!pp.contains(buffer)) {
         cout << "The input deck does not contain Nv for species " << species+1
              << "." << endl;
         cleanup(number_of_species);
         return -1;
      }
      pp.queryarr(buffer, tmp, 0, 2);
      dimVec[species][2] = tmp[0];
      dimVec[species][3] = tmp[1];
      work[species] *= tmp[0] * tmp[1];
      net_work += work[species];
   }

   // Figure out how many Vlasov processors are assigned to each species.  This
   // code is essentially a copy of what LoadBalancer does.
   if (number_of_species > nProcTotal) {
      cout << "There must be at least 1 Vlasov processor for each species."
           << endl;
      cleanup(number_of_species);
      return -1;
   }
   int next_proc = 0;
   for (int species = 0; species < number_of_species; ++species) {
      const float relative_work =
         (work[species]*static_cast<float>(nProcTotal)) / net_work;
      const int n_proportional = static_cast<int>(floor(relative_work + 0.5));
      const int n_left = nProcTotal - next_proc;
      const int n = std::min(n_proportional, n_left);
      if (n > 0) {
         nProc[species] = n;
         next_proc += n;
      }
      else {
         cout << "There are not enough Vlasov processors." << endl;
         cleanup(number_of_species);
         return -1;
      }
   }

   // Initialize dimProc to 1 for each dimension of each species.
   for (int species = 0; species < number_of_species; ++species) {
      for (int i = 0; i < nDims; ++i) {
         dimProc[species][i] = 1;
      }
   }

   // Figure out how each dimension of each species is partitioned.  This is
   // essentially Overture's GridDistribution::computeParallelArrayDistribution.
   // This SEEMS to be what Overture uses to figure out the decomposition but I
   // honestly don't see how or even where Overture uses it to decompose one of
   // our 4D arrays.  Someone said that this functionality is in BoxLib but I
   // find no trace of it there.  This is where the guessing and assuming begin.
   // Heuristically, this matches what Overture seems to do.
   for (int species = 0; species < number_of_species; ++species) {
      for (int used = 1; used < nProc[species];) {
         int bigRatioDim = 0;
         float ratio = ((float)dimVec[species][0])/((float)dimProc[species][0]);
         for (int i = 1; i < nDims; ++i)   {
            float nRatio =
               ((float)dimVec[species][i])/((float) dimProc[species][i]);
            if (nRatio > ratio ||
                (nRatio == ratio && nRatio > 1 &&
                 dimProc[species][i] < dimProc[species][bigRatioDim])) {
               ratio = nRatio;
               bigRatioDim = i;
            }
         }
         int factor = greatestPrimeFactor(nProc[species]/used);
         dimProc[species][bigRatioDim] *= factor;
         used *= factor;
      }
   }

   // For each species, find what we think will be the smallest number of zones
   // in each dimension and see if it is large enough to run.  More assumptions
   // here.  Since I couldn't find where the processor decomposition is done
   // I definitely couldn't find how that decomposition is used to actually
   // divide the work.  Again, this seems to match what Overture does.
   // Essentially, there are 2 conditions that Overture seems to try to meet:
   // 1) Equal numbers of zones (real+ghost) among all the partitions of a
   //    given direction.
   // 2) The total number of real zones has to add up.
   //
   // So we need the solution to 2 linear equations with 2 variables.
   // Partitions on the domain boundary end up with Overture's extra layer of
   // ghosts at that boundary.  So those partitions will have fewer real zones
   // and be the determining factor.
   //
   // Let numEdgeZones be the number of zones in a dimension for a domain at a
   // physical boundary and numInteriorZones be the number of zones in a
   // dimension for a domain in the interior of the mesh.  Then the equations
   // for each dimension of each species are:
   // numEdgeZones + nGhosts = numInteriorZones (condition 1)
   // 2*numEdgeZones + (dimProc - 2)*numInteriorZones = dimVec (condition 2)
   // So for a given dimension of a given species:
   // numEdgeZones = (dimVec - (dimProc - 2)*nGhosts)/dimProc
   // Overture seems to truncate numEdgeZones to an integer and then fiddle
   // with the remaining interior zones to get the total zone count right.
   bool likelyWorks = true;
   string outputString[4] = {"x", "y", "vx", "vy"};
   cout << endl << endl << endl;
   for (int species = 0; species < number_of_species; ++species) {
      cout << "Likely decomposition for species " << species + 1
           << " using " << nProc[species] << " processors:" << endl;
      for (int i = 0; i < nDims; ++i) {
         int numEdgeZones;
         // Need to handle the corner case that there is only 1 partition in
         // this dimension.  In that case we already know the number of zones.
         if (dimProc[species][i] == 1) {
            numEdgeZones = dimVec[species][i];
         }
         else {
            numEdgeZones =
               (dimVec[species][i]-(dimProc[species][i]-2)*nGhosts)/
               dimProc[species][i];
         }
         cout << "Smallest number of zones in " << outputString[i]
              << " direction is likely " << numEdgeZones << ".";
         // If the number of real zones on the domain boundary, which we think
         // is the smallest partition, is less than the stencil width then we
         // don't have a valid partition.
         if (numEdgeZones < stencilWidth) {
            cout << "  Too few interior points!" << endl;
            likelyWorks = false;
         }
         else {
            cout << endl;
         }
      }
      cout << endl;
   }
   // Summarize if this partition will or will not work.
   if (likelyWorks) {
      cout << "This decomposition is likely to work." << endl;
   }
   else {
      cout << "This decomposition is likely to fail." << endl;
   }
   cleanup(number_of_species);
   return 0;
}
