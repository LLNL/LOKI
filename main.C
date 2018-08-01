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
#define BOUNDS_CHECK

#include "Overture.h"
#include "Oges.h"
#include "ParmParse.H"

#include "ParallelUtility.h"
#include "PerlPreProcessor.H"
#include "RestartManager.H"
#include "Simulation.H"
#include "TimerManager.H"
#include "tbox/Box.H"

using namespace Loki;

/**
 * Preprocess input deck so that all perl expressions have been evaluated and
 * create resulting Loki readable file.
 */
ParmParse
createInputDatabase(
   int a_argc,
   char** a_argv);

/**
 * Free all allocated static data.
 */
void
shutdown();

/**
 * Main driver for vlasovPoisson4D.
 */
int
main(
   int a_argc,
   char* a_argv[])
{
   time_t wall_time = time(NULL);
   Overture::start(a_argc, a_argv);
   printF("\n************************************************\n");
   printF("Wall clock at start up: %s", asctime(localtime(&wall_time)));
   printF("************************************************\n\n");

   char* job_id = getenv("SLURM_JOB_ID");
   if (job_id) {
      printF("\n***************\n");
      printF("Job ID: %s\n", job_id);
      printF("***************\n");
   }

   {
      INIT_PETSC_SOLVER();
      GenericDataBase::setMaximumNumberOfFilesForWriting(16);
      ParmParse pp(createInputDatabase(a_argc, a_argv));
      Simulation simulation(pp);
      wall_time = time(NULL);
      printF("\n*****************************************************************\n");
      printF("Wall clock after problem initialization: %s",
             asctime(localtime(&wall_time)));
      printF("*****************************************************************\n\n");
      while (simulation.notDone()) {
         simulation.advance();
      }
      simulation.finalize();
   }

   Communication_Manager::Sync();
   Overture::finish();

   // This is lame and should be handled by a Manager class.  Classes needing
   // initialization/finalization should register with the Manager and the
   // Manager should be initialized/finalized.
   shutdown();
   return 0;
}


ParmParse
createInputDatabase(
   int a_argc,
   char** a_argv)
{
   // Make the pre-processor and have it process the file.  The preprocessor
   // returns the name of the processed file that Loki will actually parse.
   // Build the ParmParse object that contains the database Loki parses.
   PerlPreProcessor ppp(a_argc, a_argv);
   std::string input_file_name;
   if (a_argc == 2) {
      input_file_name = std::string(a_argv[1]);
      printF("\n\n*** Reading input file %s\n\n", input_file_name.c_str());
   }
   else {
      Overture::abort("Usage: 'srun -np N vlasovPoisson4D <inputfile>'");
   }
   std::string output_file_name(ppp.process(input_file_name));
   Communication_Manager::Sync();
   char* cstr;
   ParmParse pp(0, &cstr, "", output_file_name.c_str());
   return pp;
}

void
shutdown()
{
   // Delete all singletons, static objects, and remaining reference counted
   // objects.
   RestartManager::shutdown();
   TimerManager::shutdown();
   tbox::Box::shutdown();
   tbox::IntVector::shutdown();
   tbox::ReferenceCounter::finalizeCallback();
}
