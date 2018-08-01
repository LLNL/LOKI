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
#include "Overture.h"
#include "HDF_DataBase.h"
#include "ShowFileReader.h"
#include "ParallelUtility.h"
#include "KineticSpeciesPtrVect.H"
#include <sstream>

using namespace Loki;

void
parseCommandLine(
   aString& a_prefix,
   bool& a_skip_dists,
   int a_argc,
   char* a_argv[])
{
   int len(0);
   aString line;
   bool found_prefix = false;
   for (int i = 1; i < a_argc; ++i) {
      line = a_argv[i];

      // If this is the -prefix arg, then get the prefix and set a_prefix.
      len = line.matches("-prefix=");
      if (len) {
        a_prefix = line(len, static_cast<int>(line.length())-1);
        printF("\n\n*** using the prefix %s\n\n", (const char*)a_prefix);
        found_prefix = true;
        continue;
      }

      // If this is the -skip_dists arg, then set a_skip_dists.
      len = line.matches("-skip_dists");
      if (len) {
         a_skip_dists = true;
         continue;
      }

      // If it's neither of the above print error and exit.
      printF("Error: unknonwn command line option ... quitting1\n");
      printF("Usage: 'mpirun -np N vp4DPostProcess -prefix=<prefix> [-skip_dists]'\n");
      Overture::finish();
      exit(-1);
   }

   // The prefix is required.  If none supplied print error and exit.
   if (!found_prefix) {
      printF("Error: missing required -prefix command ... quitting1\n");
      printF("Usage: 'mpirun -np N vp4DPostProcess -prefix=<prefix> [-skip_dists]'\n");
      Overture::finish();
      exit(-1);
   }
   return;
}


void
writeTimeHistories(
   const aString& a_prefix,
   bool a_isMaxwell,
   int a_species_list_size)
{
   char buffer[100];
   aString showFileName;

   // Try to open the root show file and get the first frame which should exist.
   // If unsuccessful, exit.
   sprintf(buffer, "%s.show", (const char*)a_prefix);
   showFileName = buffer;
   printF("\nreading from %s\n", (const char*)showFileName);
   ShowFileReader show;

   show.open(showFileName);

   int numProbes;
   HDF_DataBase* dbp;
   dbp = 0;
   dbp = show.getFrame(1);
   if (dbp == 0) {
      OV_ABORT("Can not get frame from show file");
   }

   // Root show file looks good.  Read the number of probes and their locations.
   dbp->get(numProbes, "numProbes");
   printF("nProbes=%i\n", numProbes);
   RealArray fieldProbes(2, numProbes);
   IntegerArray probes_x_index(numProbes);
   IntegerArray probes_y_index(numProbes);
   char buff[12];
   for (int ip(0); ip < numProbes; ++ip) {
      sprintf(buff, "ix_pos%i", ip);
      dbp->get(probes_x_index(ip), buff);

      sprintf(buff, "iy_pos%i", ip);
      dbp->get(probes_y_index(ip), buff);
   }

   // Read the number of tracking particles.
   int numTrackingParticles;
   dbp->get(numTrackingParticles, "numTrackingParticles");
   printF("nTrackingParticles=%i\n", numTrackingParticles);

   // We're going to get the time histories from frame series 0.  Figure out
   // how many time histories there are.
   show.setCurrentFrameSeries(0);
   aString name;
   RealArray seriesTime, seriesValue;
   int maxComponentName1;
   int maxComponentName2 = 1;
   if (a_isMaxwell) {
      maxComponentName1 = 12 + (6 + a_species_list_size) * numProbes +
         numTrackingParticles * 4 + a_species_list_size * 5;
   }
   else {
      maxComponentName1 =  5 + 2 * numProbes +
         numTrackingParticles * 4 + a_species_list_size * 5;
   }
   aString componentName1[maxComponentName1], componentName2[maxComponentName2];

   // Get the time history.
   show.getSequence(0,
      name,
      seriesTime,
      seriesValue,
      componentName1,
      maxComponentName1,
      componentName2,
      maxComponentName2);

   // Create the post-processed time history file and write the probe locations,
   // time histories, and history times to it.
   sprintf(buffer, "%s_timeSeries.hdf", (const char*)a_prefix);
   HDF_DataBase db;

   db.mount(buffer, "I");                      // mount, I = initialize
   db.put(probes_x_index, "x_index");
   db.put(probes_y_index, "y_index");
   db.put(seriesValue, "series_data");
   db.put(seriesTime, "series_time");

   db.unmount();

   show.close();

   return;
}


void
writeFields(
   const aString& a_prefix,
   int nGhost,
   bool isMaxwell,
   int num_kinetic_species,
   const std::vector<string> species_names)
{
   char buffer[100];
   aString showFileName;

   // Open the root show file.
   printF("\n\n\n\n STARTING*******\n\n\n\n");
   sprintf(buffer, "%s.show", (const char*)a_prefix);
   showFileName = buffer;
   printF("\nreading from %s\n", (const char*)showFileName);
   ShowFileReader show;

   show.open(showFileName);

   CompositeGrid cg;
   realCompositeGridFunction u;
   int solnNumber = 1;
   int NX, NY;

   // Get the number of frames.  This is the number of time steps for each plot.
   // Then get the first solution (it reads the complete solution which is
   // splayed across multiple show files) and the extents of that solution which
   // must include the number of ghosts.
   int numberOfFrames = show.getNumberOfFrames();

   int nd1a, nd1b, nd2a, nd2b;
   int n1a, n1b, n2a, n2b;
   real x1a, x1b, x2a, x2b;
   real t, dt;
   show.getASolution(solnNumber, cg, u);
   cg.update();

   nd1a = u[0].getBase (axis1);
   nd1b = u[0].getBound(axis1);
   nd2a = u[0].getBase (axis2);
   nd2b = u[0].getBound(axis2);

   n1a = nd1a+nGhost;
   n1b = nd1b-nGhost;
   n2a = nd2a+nGhost;
   n2b = nd2b-nGhost;

   NX = n1b-n1a+1;
   NY = n2b-n2a+1;

   // Get the grid vertices which is common to all plots and the upper and lower
   // x and y coordinates.
   int grid = 0; // there is only 1 grid here since this is not a compisite grid example
   const realArray& vertex = cg[grid].vertex();
   x1a = vertex(n1a, n2a, 0, axis1);
   x1b = vertex(n1b, n2a, 0, axis1);
   x2a = vertex(n1a, n2a, 0, axis2);
   x2b = vertex(n1a, n2b, 0, axis2);

   // Create x and y ranges based on the number of x and y grid points.  Also
   // create a range for the number of plot times.  Then create arrays to hold
   // the x and y coordinates and the plot times.
   Range RX(0, NX-1), RY(0, NY-1), RT(0, numberOfFrames-1);
   RealArray x(RX), y(RY), time(RT);

   // Generate the coordinates.
   real dx = (x1b-x1a)/(NX-1);
   real dy = (x2b-x2a)/(NY-1);
   for (int i1 = 0; i1 < NX; ++i1) {
      x(i1) = x1a+(i1)*dx;
   }
   for (int i2 = 0; i2 < NY; ++i2) {
      y(i2) = x2a+(i2)*dy;
   }

   // Open the post-processed fields file.
   sprintf(buffer, "%s_fields.hdf", (const char*)a_prefix);
   HDF_DataBase db;

   db.mount(buffer, "I");                      // mount, I = initialize

   HDF_DataBase* dbp;

   // Generate plot names for the case of a Maxwell or a Poisson system.  This
   // is not so nice as we need to keep this sync'd up with their respective
   // plot methods.
   int numPlots;
   char** plotNames;
   if (isMaxwell) {
      numPlots = 6+5*num_kinetic_species;
      plotNames = new char* [numPlots];
      for (int plot = 0; plot < numPlots; ++plot) {
         plotNames[plot] = new char [80];
      }
      sprintf(plotNames[0], "EX");
      sprintf(plotNames[1], "EY");
      sprintf(plotNames[2], "EZ");
      sprintf(plotNames[3], "BX");
      sprintf(plotNames[4], "BY");
      sprintf(plotNames[5], "BZ");
      for (int ks = 0; ks < num_kinetic_species; ++ks) {
         sprintf(plotNames[6+5*ks], "%s VZ",
                 species_names[ks].c_str());
         sprintf(plotNames[7+5*ks], "%s ke flux vx lo",
                 species_names[ks].c_str());
         sprintf(plotNames[8+5*ks], "%s ke flux vx hi",
                 species_names[ks].c_str());
         sprintf(plotNames[9+5*ks], "%s ke flux vy lo",
                 species_names[ks].c_str());
         sprintf(plotNames[10+5*ks], "%s ke flux vy hi",
                 species_names[ks].c_str());
      }
   }
   else {
      numPlots = 2+4*num_kinetic_species;
      plotNames = new char* [numPlots];
      for (int plot = 0; plot < numPlots; ++plot) {
         plotNames[plot] = new char [80];
      }
      sprintf(plotNames[0], "EX");
      sprintf(plotNames[1], "EY");
      for (int ks = 0; ks < num_kinetic_species; ++ks) {
         sprintf(plotNames[2+4*ks], "%s ke flux vx lo",
                 species_names[ks].c_str());
         sprintf(plotNames[3+4*ks], "%s ke flux vx hi",
                 species_names[ks].c_str());
         sprintf(plotNames[4+4*ks], "%s ke flux vy lo",
                 species_names[ks].c_str());
         sprintf(plotNames[5+4*ks], "%s ke flux vy hi",
                 species_names[ks].c_str());
      }
   }

   // Now we know how may plots to read.  Loop through them and read each time
   // step for each plot.
   RealArray plotData(RT, RY, RX);
   for (int plot = 0; plot < numPlots; ++plot) {
      show.setCurrentFrameSeries(plot);
      for (int soln = 1; soln <= numberOfFrames; ++soln) {
         solnNumber = soln;
         show.getASolution(solnNumber, cg, u);

         // We need to get the plot times and dt once for the first plot.
         // Seems like this could be done outside these loops.  Also seems like
         // we never use dt.
         if (plot == 0) {
            dbp = 0;
            dbp = show.getFrame(soln);
            if (dbp == 0) {
               OV_ABORT("Can not get frame from show file");
            }
            dbp->get(t, "time");
            dbp->get(dt, "dt");
            time(soln-1) = t;
         }

         // Move the grid function plot data into our plotData array.
         for (int i1 = 0; i1 < NX; ++i1) {
            for (int i2 = 0; i2 < NY; ++i2) {
               plotData(soln-1, i2, i1) = u[0](i1, i2);
            }
         }
      }
      // Write this plot to the post-processed fields file.  Write the time and
      // coordinates once as they are common to all plots.
      if (plot == 0) {
         db.put(time,   "time");
         db.put(x,      "x");
         db.put(y,      "y");
      }
      db.put(plotData, plotNames[plot]);
   }

   // Close files and clean up.
   db.unmount();

   show.close();

   for (int plot = 0; plot < numPlots; ++plot) {
      delete [] plotNames[plot];
   }
   delete [] plotNames;

   return;
}

/**
 * Main driver for vp4DPostProcess.  There are only 2 command line arg:
 * -prefix, this is the restart.write_directory entry in the input file
 * -skip_dists, optionally skip processing the distribution functions
 */
int
main(
   int a_argc,
   char* a_argv[])
{
   // initialize Overture
   Overture::start(a_argc, a_argv);
   GenericDataBase::setMaximumNumberOfFilesForWriting(1);

   bool isMaxwell = true;
   int species_list_size;
   std::vector<string> species_names;
   // These will be read in so their initializations probably are not needed.
   int nGhost = 3;
   int spatial_solution_order = 6;
   int temporal_solution_order = 6;

   // Figure out what we're reading and how to process it.
   aString prefix;
   bool skip_dists = false;
   parseCommandLine(prefix, skip_dists, a_argc, a_argv);

   // Maximum number of check points of the distribution function we will read.
   // This is probably here just to put a cap on the amount of data this thing
   // will write.
   const int iMax(1000);

   if (!skip_dists) {
      // Loop over the maximum allowed number of checkpoint files.
      printF("starting with distribution functions\n");
      for (int i(0); i < iMax; ++i) {

         // Try to open this checkpoint file.  If we can't then we're done
         // with the checkpoint files that have been written.
         char buffer[100];
         sprintf(buffer, "%s/dist_%i.hdf", (const char*)prefix, i);

         HDF_DataBase db;
         if (db.mount(buffer, "R") != 0) { // mount, R = read-only
            break;
         }
         else {
            printF("opening dist_%i\n", i);

            // Get the time stamp of this checkpoint and the number of ghosts.
            // Then figure out the solution orders.
            real t(0.0);
            db.get(t, "time");

            db.get(nGhost, "nGhost");

            if (nGhost == 2) {
               spatial_solution_order = 4;
               temporal_solution_order = 4;
            }
            else {
               spatial_solution_order = 6;
               temporal_solution_order = 6;
            }

            // See if this is a Maxwell run.
            HDF_DataBase maxwell_db;
            if (db.locate(maxwell_db, "Maxwell") != 0) {
               isMaxwell = false;
            }

            // Get the list of species names.
            HDF_DataBase sub_db;
            if (db.locate(sub_db, "species_list") != 0) {
               std::cout << "ERROR: vp4DPostProcess.C: Unable to locate 'species_list' GROUP"
                         << std::endl;
               exit(-1);
            }
            std::cout << "Located species_list sub_db..." << std::endl;
            if (db.get(species_list_size, "species_list_size") != 0) {
               std::cout << "ERROR: vp4DPostProcess.C: Unable to locate 'species_list_size' DATASET"
                         << std::endl;
               exit(-1);
            }
            std::cout << "Obtained species_list_size of "
                      << species_list_size << std::endl;

            for (int s(0); s < species_list_size; ++s) {
               aString tmp;
               std::stringstream tag;
               tag << "species." << s + 1;
               sub_db.get(tmp, (aString)tag.str());
               species_names.push_back(tmp);
            }

            // Construct each species.
            KineticSpeciesPtrVect kinetic_species;
            kinetic_species.resize(species_list_size);
            for (int s(0); s < species_list_size; ++s) {
               kinetic_species[s] =
                  new KineticSpecies(db,
                                     s+1,
                                     species_list_size,
                                     spatial_solution_order,
                                     temporal_solution_order,
                                     false,
                                     species_names[s]);
            }
            std::cout << "**Created vector of new KineticSpecies" << std::endl;

            // Done reading from this checkpoint file.
            db.unmount();

            // Open the post-processed distribution function file for this
            // checkpoint and write the distribution functions to it.
            sprintf(buffer, "%s_dist_%i.hdf", (const char*)prefix, i);
            HDF_DataBase db2;
            db2.mount(buffer, "I");                   // mount, I = initialize

            // Save the time stamp and the number of ghosts.
            db2.put(t, "time");
            db2.put(nGhost, "nGhost");

            // Write the number of species.
            HDF_DataBase sub_db2;
            db2.create(sub_db2, "species_list", "directory");
            db2.put(species_list_size, "species_list_size");

            // Have each species write it distribution function and other
            // basic parameters.  We don't care about Krook layers.
            for (int s(0); s < species_list_size; ++s) {
               kinetic_species[s]->putToRestart_SkipKrook(db2, 0.0);
               // TODO: Replace below with call to VPSystem::putToRestart()...?
               std::stringstream tag;
               tag << "species." << s + 1;
               sub_db2.put(species_names[s], tag.str().c_str());
            }
            std::cout << "**Put KineticSpecies vector to Restart" << std::endl;
            db2.unmount();
         }

         std::cout << "Completed iteration " << i << " of loop" << std::endl;
      }
      printF("\nfinished with distribution functions\n");
   }
   else {
      // If we're not processing the distribution functions all we need to do
      // is to get some basic info from the first checkpoint file.

      // Open the first checkpoint file and get the number of ghosts.
      char buffer[100];
      sprintf(buffer, "%s/dist_0.hdf", (const char*)prefix);

      HDF_DataBase db;
      db.mount(buffer, "R");
      db.get(nGhost, "nGhost");

      // See if this is a Maxwell run.
      HDF_DataBase maxwell_db;
      if (db.locate(maxwell_db, "Maxwell") != 0) {
         isMaxwell = false;
      }

      // Get list of species names.
      HDF_DataBase sub_db;
      if (db.locate(sub_db, "species_list") != 0) {
         std::cout << "ERROR: vp4DPostProcess.C: Unable to locate 'species_list' GROUP"
                   << std::endl;
         exit(-1);
      }
      std::cout << "Located species_list sub_db..." << std::endl;
      if (db.get(species_list_size, "species_list_size") != 0) {
         std::cout << "ERROR: vp4DPostProcess.C: Unable to locate 'species_list_size' DATASET"
                   << std::endl;
         exit(-1);
      }
      std::cout << "Obtained species_list_size of "
                << species_list_size << std::endl;

      for (int s(0); s < species_list_size; ++s) {
         aString tmp;
         std::stringstream tag;
         tag << "species." << s + 1;
         sub_db.get(tmp, (aString)tag.str());
         species_names.push_back(tmp);
      }
   }

   // Finished with distribution functions so process the time histories and
   // fields.
   writeTimeHistories(prefix, isMaxwell, species_list_size);

   writeFields(prefix, nGhost, isMaxwell, species_list_size, species_names);

   Overture::finish();
   return 0;
}
