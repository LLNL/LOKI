/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include "Loki_Utilities.H"
#include "LokiInputParser.H"
#include "PerlPreProcessor.H"

using namespace std;
using namespace Loki;

/*
 * Serial utility to find the smallest and largest number of interior zones
 * owned by a processor in each dimension after ParallelArray decomposes the
 * problem.  The utility issues a warning if there are fewer than the minimum
 * number of interior zones necessary for the Loki algorithms.  Note that the
 * information reported by this utility is different that the max and min
 * dimensions reported by the main code.  That includes the ghost and boundary
 * zones.  This is pretty much what is done in ParallelArray::partition.
 */

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
      cout << "greatestPrimeFactor:ERROR: can't factor a negative number"
           << endl;
      return n;
   }

   // After dividing out all the smaller prime factors, what's left is the
   // largest one.
   for (int lpf = leastPrimeFactor(n); lpf != n; lpf = leastPrimeFactor(n)) {
      n = n/lpf;
   }

   return n;
}

// Finds the smallest and largest number of zones owned by a processor in each
// dimension after ParallelArray decomposes the problem.  This utility issues a
// warning if there are fewer than the minimum number of zones necessary for the
// Loki algorithms.  It takes 2 input arguments, the total number of Vlasov
// processors and the input deck for the problem of interest.
int
main(int argc, char* argv[])
{
   MPI_Init(&argc, &argv);
   Loki_Utilities::initialize();

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
   Loki_Utilities::printF("\n\n*** Reading input file %s\n\n",
      input_file_name.c_str());
   string output_file_name(ppp.process(input_file_name));
   char* cstr;
   LokiInputParser pp(0, &cstr, "", output_file_name.c_str());

   // Read the spatial solution order from the input deck and set stencilWidth
   // based on that value.
   int spatialOrder = 4;
   pp.query("spatial_solution_order", spatialOrder);
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
   vector<vector<int> > dimVec(number_of_species);
   vector<vector<int> > dimProc(number_of_species);
   vector<int> nProc(number_of_species);
   vector<float> work(number_of_species);
   for (int species = 0; species < number_of_species; ++species) {
      dimVec[species].resize(nDims);
      dimProc[species].resize(nDims);
   }

   // Read the number of spatial points from the input deck.  This is the
   // size of the first 2 dimensions of each species.  Also, start computing
   // species work.
   if (!pp.contains("N")) {
      cout << "The input deck does not contain the number of spatial points."
           << endl;
      return -1;
   }
   vector<int> tmp(2);
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
      ostringstream input_string;
      input_string << "kinetic_species." << species+1 << ".Nv";
      if (!pp.contains(input_string.str().c_str())) {
         cout << "The input deck does not contain Nv for species " << species+1
              << "." << endl;
         return -1;
      }
      pp.queryarr(input_string.str().c_str(), tmp, 0, 2);
      input_string.str("");
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
      return -1;
   }
   int next_proc = 0;
   for (int species = 0; species < number_of_species; ++species) {
      const float relative_work =
         (work[species]*static_cast<float>(nProcTotal)) / net_work;
      const int n_proportional = static_cast<int>(floor(relative_work + 0.5));
      const int n_left = nProcTotal - next_proc;
      const int n = min(n_proportional, n_left);
      if (n > 0) {
         nProc[species] = n;
         next_proc += n;
      }
      else {
         cout << "There are not enough Vlasov processors." << endl;
         return -1;
      }
   }

   // Figure out how many partitions there are in each dimension of each
   // species.  This is essentially the assignProcs method in ParallelArray.

   // Initialize dimProc to 1 for each dimension of each species.
   for (int species = 0; species < number_of_species; ++species) {
      for (int i = 0; i < nDims; ++i) {
         dimProc[species][i] = 1;
      }
   }

   // For each species assign the prime factors of its nProc to the 4
   // dimensions.  This tries to create as "square" a partition as possible
   // given the problem dimensions and the prime factorization of the number
   // of processors.
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

   // For each species, find the smallest and largest number of zones in each
   // dimension and see if the smallest is large enough to run the problem.
   // For each species all we need to do is to divide dimVec by dimProc for each
   // dimension.  This is the minimum number of zones in each dimension.  If
   // the division in a dimension is exact, the maximum number of zone in that
   // dimension is the minimum.  Otherwise the maximum number of zones in that
   // dimension is one greater than the minimum.
   bool works = true;
   string outputString[4] = {"x", "y", "vx", "vy"};
   cout << endl << endl << endl;
   for (int species = 0; species < number_of_species; ++species) {
      cout << "Decomposition for species " << species + 1 << " using "
           << nProc[species] << " processors:" << endl;
      for (int i = 0; i < nDims; ++i) {
         cout << dimProc[species][i] << " partitions in " << outputString[i]
              << " direction." << endl;
      }
      for (int i = 0; i < nDims; ++i) {
         int minZones = dimVec[species][i]/dimProc[species][i];
         int maxZones = minZones;
         if (dimVec[species][i]%dimProc[species][i] != 0) {
            ++maxZones;
         }
         cout << "Smallest number of interior zones in " << outputString[i]
              << " direction is " << minZones << ".";
         // If the smallest partition, is less than the stencil width then we
         // don't have a valid partition.
         if (minZones < stencilWidth) {
            cout << "  Too few interior points!" << endl;
            works = false;
         }
         else {
            cout << endl;
         }
         cout << "Largest number of interior zones in " << outputString[i]
              << " direction is " << maxZones << "." << endl;
      }
      cout << endl;
   }
   // Summarize if this partition will or will not work.
   if (works) {
      cout << "This decomposition will work." << endl << endl;
   }
   else {
      cout << "This decomposition will fail." << endl << endl;
   }

   // Suggest a number of EM processors.
   int suggestion = 0;
   for (int species = 0; species < number_of_species; ++species) {
      suggestion = max(suggestion, dimProc[species][0]*dimProc[species][1]);
   }
   suggestion = ceil(suggestion/50.0);
   cout << "To reduce communication bottlenecks between the EM and" << endl;
   cout << "Vlasov processors you should use enough EM processors" << endl;
   cout << "such that (number of x partitions * number of y partitions)" << endl;
   cout << "is about 50x the number of EM processors.  For this problem" << endl;
   cout << "that number is " << suggestion << "." << endl;
   cout << "You should adjust this value in order to get something with" << endl;
   cout << "prime factors suitable to the number of x and y points." << endl;
   MPI_Finalize();
   return 0;
}
