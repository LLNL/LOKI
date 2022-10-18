/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "PerlPreProcessor.H"
#include "Loki_Utilities.H"
#include "RestartManager.H"
#include "Simulation.H"
#include "TimerManager.H"
#include "tbox/Box.H"

using namespace Loki;

/**
 * Preprocess input deck so that all perl expressions have been evaluated and
 * create resulting Loki readable file.
 */
LokiInputParser
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
   MPI_Init(&a_argc, &a_argv);
   Loki_Utilities::initialize();
   Loki_Utilities::printF("\n************************************************\n");
   Loki_Utilities::printF("Wall clock at start up: %s",
      asctime(localtime(&wall_time)));
   Loki_Utilities::printF("************************************************\n\n");

   char* job_id = getenv("SLURM_JOB_ID");
   if (job_id) {
      Loki_Utilities::printF("\n***************\n");
      Loki_Utilities::printF("Job ID: %s\n", job_id);
      Loki_Utilities::printF("***************\n");
   }

   {
      LokiInputParser pp(createInputDatabase(a_argc, a_argv));
      Simulation simulation(pp);
      wall_time = time(NULL);
      Loki_Utilities::printF("\n*****************************************************************\n");
      Loki_Utilities::printF("Wall clock after problem initialization: %s",
         asctime(localtime(&wall_time)));
      Loki_Utilities::printF("*****************************************************************\n\n");
      while (simulation.notDone()) {
         simulation.advance();
      }
      simulation.finalize();
   }

   MPI_Finalize();

   // This is lame and should be handled by a Manager class.  Classes needing
   // initialization/finalization should register with the Manager and the
   // Manager should be initialized/finalized.
   shutdown();

   wall_time = time(NULL);
   Loki_Utilities::printF("\n************************************************\n");
   Loki_Utilities::printF("Wall clock at shut down: %s",
      asctime(localtime(&wall_time)));
   Loki_Utilities::printF("************************************************\n\n");
   return 0;
}


LokiInputParser
createInputDatabase(
   int a_argc,
   char** a_argv)
{
   // Make the pre-processor and have it process the file.  The preprocessor
   // returns the name of the processed file that Loki will actually parse.
   // Build the LokiInputParser object that contains the database Loki parses.
   PerlPreProcessor ppp(a_argc, a_argv);
   string input_file_name;
   if (a_argc == 2) {
      input_file_name = string(a_argv[1]);
      Loki_Utilities::printF("\n\n*** Reading input file %s\n\n",
         input_file_name.c_str());
   }
   else {
      LOKI_ABORT("Usage: 'srun -np N vlasovPoisson4D <inputfile>'");
   }
   string output_file_name(ppp.process(input_file_name));
   fflush(0);
   MPI_Barrier(MPI_COMM_WORLD);
   char* cstr;
   LokiInputParser pp(0, &cstr, "", output_file_name.c_str());
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
