/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "RestartReader.H"
#include "RestartWriter.H"
#include "FieldReader.H"
#include "FieldWriter.H"
#include "TimeHistReader.H"
#include "TimeHistWriter.H"
#include "KineticSpeciesPtrVect.H"
#include "Maxwell.H"
#include "Poisson.H"
#include <deque>
#include <sstream>
#include <sys/stat.h>

using namespace Loki;

void
parseCommandLine(
   string& a_prefix,
   bool& a_skip_dists,
   deque<int>& a_process_dists,
   bool& a_coll,
   int a_argc,
   char* a_argv[])
{
   int len(0);
   string line;
   bool found_prefix = false;
   for (int i = 1; i < a_argc; ++i) {
      line = a_argv[i];

      // If this is the -prefix arg, then get the prefix and set a_prefix.
      if (line.substr(0, 8) == "-prefix=") {
        a_prefix = line.substr(8);
        Loki_Utilities::printF("\n\n*** using the prefix %s\n\n",
          a_prefix.c_str());
        found_prefix = true;
        continue;
      }

      // If this is the -skip_dists arg, then set a_skip_dists.
      if (line.substr(0, 11) == "-skip_dists") {
         a_skip_dists = true;
         continue;
      }

      // If this is the -process_dists arg, then set the code up to process
      // specific distribution dumps.
      if (line.substr(0, 15) == "-process_dists=") {
         if (a_skip_dists) {
            Loki_Utilities::printF("Error: Attempting to both skip and process specific dists\n");
            exit(-1);
         }
         string dist_list = line.substr(15);
         while (!dist_list.empty()) {
            size_t pos = dist_list.find(",");
            string dist_num = dist_list.substr(0, pos);
            a_process_dists.push_back(atoi(dist_num.c_str()));
            if (pos == string::npos) {
               dist_list.erase(0, pos);
            }
            else {
               dist_list.erase(0, pos+1);
            }
         }
         continue;
      }

      // If this is the -coll arg, then set a_coll.
      if (line.substr(0, 5) == "-coll") {
         a_coll = true;
         continue;
      }

      // If it's neither of the above print error and exit.
      Loki_Utilities::printF("Error: unknonwn command line option ... quitting\n");
      Loki_Utilities::printF("Usage: 'vp4DPostProcess -prefix=<prefix> [-skip_dists]'\n");
      exit(-1);
   }

   // The prefix is required.  If none supplied print error and exit.
   if (!found_prefix) {
      Loki_Utilities::printF("Error: missing required -prefix command ... quitting\n");
      Loki_Utilities::printF("Usage: 'vp4DPostProcess -prefix=<prefix> [-skip_dists]'\n");
      exit(-1);
   }
   return;
}


void
writeTimeHistories(
   const string& a_prefix,
   int a_isMaxwell,
   bool a_isCollision,
   vector<string>& a_species_names)
{
   // Find the last time history file written.  This contains all the
   // accumulated time histories.  Read them and write to the post processed
   // time history file.
   vector<vector<double> > time_hists;
   vector<double> time_hist_times;
   vector<string> time_hist_names;
   int time_hist_idx = 0;
   int tmp_idx = 0;
   ostringstream this_time_hist_file_name;
   Loki_Utilities::printF("\n\n\n\n STARTING TIME HISTORIES*******\n\n\n\n");

   // Find the last time history file.
   while (true) {
      // Construct a time history file name and try to open it.  If it can't
      // be opened then we've exhausted all the time history files.
      struct stat stat_buf;
      if (a_isCollision) {
         this_time_hist_file_name << a_prefix << "_coll.time_hists_" << tmp_idx
                                  << ".hdf";
      }
      else {
         this_time_hist_file_name << a_prefix << ".time_hists_" << tmp_idx
                                  << ".hdf";
      }
      if (stat(this_time_hist_file_name.str().c_str(), &stat_buf) == 0) {
         time_hist_idx = tmp_idx;
         ++tmp_idx;
         this_time_hist_file_name.str("");
      }
      else {
         this_time_hist_file_name.str("");
         break;
      }
   }

   // Now open the last time history file.
   if (a_isCollision) {
      this_time_hist_file_name << a_prefix << "_coll.time_hists_"
                               << time_hist_idx << ".hdf";
   }
   else {
      this_time_hist_file_name << a_prefix << ".time_hists_"
                               << time_hist_idx << ".hdf";
   }
   TimeHistReader th_reader(this_time_hist_file_name.str());

   // Read the number of probes and the number of tracking particles from the
   // file and figure out the number of time histories.  Allocate the vectors
   // of post-processed time histories and time history names.  Fill in the
   // time history names.
   int numProbes = 0, numTrackingParticles = 0;
   if (!a_isCollision) {
      th_reader.readNumProbes(numProbes);
      th_reader.readNumTrackingParticles(numTrackingParticles);
   }
   if (a_isMaxwell == 1) {
      if (a_isCollision) {
         Maxwell::buildCollTimeHistoryNames(a_species_names, time_hist_names);
      }
      else {
         Maxwell::buildTimeHistoryNames(numProbes,
                                        numTrackingParticles,
                                        a_species_names,
                                        time_hist_names);
      }
   }
   else {
      if (a_isCollision) {
         Poisson::buildCollTimeHistoryNames(a_species_names, time_hist_names);
      }
      else {
         Poisson::buildTimeHistoryNames(numProbes,
                                        numTrackingParticles,
                                        a_species_names,
                                        time_hist_names);
      }
   }

   // Read from this time history file each time sequence and the times at
   // which they were taken.
   int num_time_hists = static_cast<int>(time_hist_names.size());
   time_hists.resize(num_time_hists);
   for (int th = 0; th < num_time_hists; ++th) {
      Loki_Utilities::printF("Processing time history %s\n",
         time_hist_names[th].c_str());
      th_reader.readTimeHistory(time_hist_names[th], time_hists[th]);
   }
   th_reader.readTimeHistory("sequence_times", time_hist_times);

   // Create the post processed time history file.
   ostringstream pp_time_hist_file_name;
   if (a_isCollision) {
      pp_time_hist_file_name << a_prefix << "_coll_timeSeries.hdf";
   }
   else {
      pp_time_hist_file_name << a_prefix << "_timeSeries.hdf";
   }
   TimeHistWriter th_writer(pp_time_hist_file_name.str(), time_hist_times);

   // Write each time history to the post_processed time history file.
   for (int th = 0; th < num_time_hists; ++th) {
      th_writer.writeTimeHistory(time_hist_names[th], time_hists[th]);
   }

   return;
}


void
writeFields(
   const string& a_prefix,
   int nx,
   int ny,
   int a_nGhost,
   int a_isMaxwell,
   bool a_isCollision,
   const vector<string>& a_species_names,
   int a_plot_ke_vel_bdy_flux)
{
   Loki_Utilities::printF("\n\n\n\n STARTING FIELDS*******\n\n\n\n");
   int nx_w_ghost = nx+2*a_nGhost;

   // Try to open the first field file and get the total number of time slices
   // from it.
   ostringstream file_name;
   if (a_isCollision) {
      file_name << a_prefix << "_coll.fields_" << 0 << ".hdf";
   }
   else {
      file_name << a_prefix << ".fields_" << 0 << ".hdf";
   }
   int total_num_time_slices;
   {
      FieldReader field_reader(file_name.str());
      field_reader.readTotalNumTimeSlices(total_num_time_slices);
   }
   file_name.str("");

   // See how many field files exist.
   int num_field_files = 0;
   while(true) {
      struct stat stat_buf;
      if (a_isCollision) {
         file_name << a_prefix << "_coll.fields_" << num_field_files << ".hdf";
      }
      else {
         file_name << a_prefix << ".fields_" << num_field_files << ".hdf";
      }
      if (stat(file_name.str().c_str(), &stat_buf) != 0) {
         break;
      }
      ++num_field_files;
      file_name.str("");
   }
   file_name.str("");

   // Generate plot names depending on if it is a Maxwell or a Poisson system.
   vector<string> plot_names;
   if (a_isMaxwell == 1) {
      if (a_isCollision) {
         Maxwell::buildCollPlotNames(a_species_names, plot_names);
      }
      else {
         Maxwell::buildPlotNames(a_plot_ke_vel_bdy_flux == 1,
                                 a_species_names,
                                 plot_names);
      }
   }
   else {
      if (a_isCollision) {
         Poisson::buildCollPlotNames(a_species_names, plot_names);
      }
      else {
         Poisson::buildPlotNames(a_plot_ke_vel_bdy_flux == 1,
                                 a_species_names,
                                 plot_names);
      }
   }
   int num_plots = static_cast<int>(plot_names.size());

   // Create the writer for the post-processed field file.
   if (a_isCollision) {
      file_name << a_prefix << "_coll_fields.hdf";
   }
   else {
      file_name << a_prefix << "_fields.hdf";
   }
   FieldWriter pp_writer(file_name.str(),
      plot_names,
      nx,
      ny,
      total_num_time_slices);
   file_name.str("");

   // Now we know how many plots to read from each field file.  Loop through the
   // field files, read each plot and plot time and save them to the post
   // processed field file.
   vector<double> plot_time_slice;
   double* plot_time_slice_no_ghost = new double[nx*ny];
   vector<double> plot_times;
   int which_time_slice = 0;
   for (int file = 0; file < num_field_files; ++file) {
      // Get a reader for a field file.
      if (a_isCollision) {
         file_name << a_prefix << "_coll.fields_" << file << ".hdf";
      }
      else {
         file_name << a_prefix << ".fields_" << file << ".hdf";
      }
      {
         FieldReader field_reader(file_name.str());
         file_name.str("");

         // We only need to get the x and y coordinates from the 1st field file.
         if (file == 0) {
            vector<double> x_coords;
            vector<double> y_coords;
            field_reader.readCoords(x_coords, y_coords);
            pp_writer.writeCoords(x_coords, y_coords);
         }

         // Read the number of time slices in this field file.
         int num_time_slices_in_this_file;
         field_reader.readNumTimeSlicesInFile(num_time_slices_in_this_file);

         // Now loop over all the time slices in this field file.
         // Read the info for each time slice and write it to the post processed
         // field file.
         for (int time_slice = 0;
              time_slice < num_time_slices_in_this_file;
              ++time_slice) {
            double time;
            ostringstream dset_name;
            dset_name << "time_slice_" << which_time_slice << "_time";
            field_reader.readTime(dset_name.str(), time);
            dset_name.str("");
            plot_times.push_back(time);

            // Loop through the plots for this time slice.
            for (int plot = 0; plot < num_plots; ++plot) {
               Loki_Utilities::printF("Processing time slice %d of field %s\n",
                  which_time_slice,
                  plot_names[plot].c_str());

               // Read plot data.
               dset_name << "time_slice_" << which_time_slice << "_"
                         << plot_names[plot];
               field_reader.readField(dset_name.str(), plot_time_slice);
               dset_name.str("");

               // Strip out the ghosts from the plot data.
               for (int iy = 0; iy < ny; ++iy) {
                  for (int ix = 0; ix < nx; ++ix) {
                     plot_time_slice_no_ghost[iy*nx+ix] =
                        plot_time_slice[(iy+a_nGhost)*nx_w_ghost + ix+a_nGhost];
                  }
               }

               // Now write this plot's time slice to the post processed file.
               pp_writer.writeFieldTimeSlice(plot_names[plot],
                  nx,
                  ny,
                  which_time_slice,
                  plot_time_slice_no_ghost);
            }

            ++which_time_slice;
            if (which_time_slice == total_num_time_slices) {
               goto bailout;
            }
         }
      }
   }
   bailout:
   delete [] plot_time_slice_no_ghost;

   // Write the plot times.
   pp_writer.writePlotTimes(plot_times);
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
   // initialize MPI
   MPI_Init(&a_argc, &a_argv);
   Loki_Utilities::initialize();

   int isMaxwell;
   int species_list_size;
   vector<string> species_names;
   // These will be read in so their initializations probably are not needed.
   int nGhost = 3;
   int spatial_solution_order = 6;
   int temporal_solution_order = 6;

   int major_version;
   int minor_version;
   int patch_level;
   int nx, ny;

   // For backward compatibility initilize this to true.  If it is not in the
   // restart file then we're processing an older file which will contain these
   // plots.
   int plot_ke_vel_bdy_flux = 1;

   // Figure out what we're reading and how to process it.
   string prefix;
   bool skip_dists = false;
   deque<int> process_dists(0);
   bool coll = false;
   parseCommandLine(prefix, skip_dists, process_dists, coll, a_argc, a_argv);

   if (!skip_dists) {
      // Loop over either all the checkpoint files or, if specific checkpoints
      // have been requested, only those requested.
      Loki_Utilities::printF("starting with distribution functions\n");
      bool first_dist = true;
      bool processing_specific_dists = process_dists.size() != 0;

      // Set the index of the first checkpoint file to be processed.
      int idx;
      if (processing_specific_dists) {
         idx = process_dists.front();
         process_dists.pop_front();
      }
      else {
         idx = 0;
      }
      while (true) {

         // Check if this checkpoint file exists.  If it does, then process it.
         // If not then we're done with the checkpoint files that have been
         // written.
         ostringstream base_file_name;
         base_file_name << prefix << "/dist_" << idx << ".hdf";
         struct stat stat_buf;
         if (stat(base_file_name.str().c_str(), &stat_buf) != 0) {
            break;
         }
         else {
            // Find out how many bulkddata files exist.
            ostringstream bulkdata_file_name;
            int bulkdata_idx = 0;
            while (true) {
               bulkdata_file_name << base_file_name.str();
               bulkdata_file_name << ".g" << bulkdata_idx;
               if (stat(bulkdata_file_name.str().c_str(), &stat_buf) != 0) {
                  break;
               }
               ++bulkdata_idx;
               bulkdata_file_name.str("");
            }
            Loki_Utilities::printF("opening dist_%i\n", idx);
            RestartReader reader(base_file_name.str(), bulkdata_idx, true);
            base_file_name.str("");

            // We only need to check that the post-processor is compatible with
            // the version of Loki that generated this data one time.
            if (first_dist) {
               // This version of the post-processor will not work with a
               // version of Loki earlier than 3.0.0.
               reader.readIntegerValue("major version", major_version);
               reader.readIntegerValue("minor version", minor_version);
               reader.readIntegerValue("patch level", patch_level);
               if (major_version < 3 || minor_version < 0 || patch_level < 0) {
                  cout << "ERROR: vp4DPostProcess.C: Data generated with incompatible version of Loki."
                       << endl;
                  exit(-1);
               }
            }

            // Get the time stamp of this checkpoint and the number of ghosts.
            // Then figure out the solution orders.
            double t(0.0);
            reader.readDoubleValue("time", t);

            reader.readIntegerValue("nGhost", nGhost);

            reader.readIntegerValue("plot_ke_vel_bdy_flux",
               plot_ke_vel_bdy_flux);

            if (nGhost == 2) {
               spatial_solution_order = 4;
               temporal_solution_order = 4;
            }
            else {
               spatial_solution_order = 6;
               temporal_solution_order = 6;
            }

            // See if this is a Maxwell run.
            reader.readIntegerValue("isMaxwell", isMaxwell);

            // Get the list of species names.
            reader.readIntegerValue("species_list_size", species_list_size);
            cout << "Obtained species_list_size of "
                 << species_list_size << endl;
            reader.pushSubDir("species_list");
            cout << "Located species_list subdir..." << endl;

            // We only need to get the species names one time.
            if (first_dist) {
               for (int s(0); s < species_list_size; ++s) {
                  string tmp;
                  stringstream tag;
                  tag << "species." << s + 1;
                  reader.readString(tag.str(), tmp);
                  species_names.push_back(tmp);
               }
            }
            reader.popSubDir();

            // Construct each species.
            KineticSpeciesPtrVect kinetic_species;
            kinetic_species.resize(species_list_size);
            for (int s(0); s < species_list_size; ++s) {
               kinetic_species[s] =
                  new KineticSpecies(reader,
                     s+1,
                     species_list_size,
                     spatial_solution_order,
                     temporal_solution_order,
                     plot_ke_vel_bdy_flux,
                     false,
                     species_names[s]);
               if (first_dist && s == 0) {
                  nx = kinetic_species[s]->domain().numberOfCells(0);
                  ny = kinetic_species[s]->domain().numberOfCells(1);
               }
            }
            cout << "**Created vector of new KineticSpecies" << endl;

            // Open the post-processed distribution function file for this
            // checkpoint and write the distribution functions to it.
            base_file_name << prefix << "_dist_" << idx << ".hdf";
            RestartWriter writer(base_file_name.str(), bulkdata_idx, true);
            base_file_name.str("");

            // Save the time stamp, the number of ghosts, and whether the KE
            // velocity boundary fluxes were written.
            writer.writeDoubleValue("time", t, true);
            writer.writeIntegerValue("nGhost", nGhost, true);
            writer.writeIntegerValue("plot_ke_vel_bdy_flux",
               plot_ke_vel_bdy_flux,
               true);

            // Write the number of species.
            writer.writeIntegerValue("species_list_size",
               species_list_size,
               true);

            // Write the names of each species.
            writer.pushSubDir("species_list");
            for (int s(0); s < species_list_size; ++s) {
               stringstream tag;
               tag << "species." << s + 1;
               writer.writeString(tag.str().c_str(), species_names[s], true);
            }
            writer.popSubDir();

            // Have each species write it distribution function and other
            // basic parameters.  We don't care about Krook layers.
            for (int s(0); s < species_list_size; ++s) {
               kinetic_species[s]->putToRestart_SkipKrook(writer, 0.0);
            }
            cout << "**Processed KineticSpecies for restart " << idx << endl;
         }

        // Update the index of the checkpoint file to be processed next.  If
        // specific checkpoint files have been specified then get the next
        // index from that list.  If the list is exhausted then we're done.
        // Otherwise just increment the index and try to process it.
         if (processing_specific_dists) {
            if (process_dists.empty()) {
               break;
            }
            idx = process_dists.front();
            process_dists.pop_front();
         }
         else {
            ++idx;
         }
         first_dist = false;
      }
      Loki_Utilities::printF("\nfinished with distribution functions\n");
   }
   else {
      // If we're not processing the distribution functions all we need to do
      // is to get some basic info from the first checkpoint file.

      // Open the first checkpoint file and get the number of ghosts.
      ostringstream base_file_name;
      base_file_name << prefix << "/dist_0.hdf";

      RestartReader reader(base_file_name.str(), 1, true);
      reader.readIntegerValue("nGhost", nGhost);

      reader.readIntegerValue("plot_ke_vel_bdy_flux", plot_ke_vel_bdy_flux);

      // See if this is a Maxwell run.
      reader.readIntegerValue("isMaxwell", isMaxwell);

      // Get list of species names.
      reader.readIntegerValue("species_list_size", species_list_size);
      cout << "Obtained species_list_size of " << species_list_size << endl;
      reader.pushSubDir("species_list");
      cout << "Located species_list sub_db..." << endl;

      for (int s(0); s < species_list_size; ++s) {
         string tmp;
         stringstream tag;
         tag << "species." << s + 1;
         reader.readString(tag.str(), tmp);
         species_names.push_back(tmp);
      }
      reader.popSubDir();

      // This version of the post-processor will not work with a version of Loki
      // earlier than 3.0.0.
      reader.readIntegerValue("major version", major_version);
      reader.readIntegerValue("minor version", minor_version);
      reader.readIntegerValue("patch level", patch_level);
      if (major_version < 3 || minor_version < 0 || patch_level < 0) {
         cout << "ERROR: vp4DPostProcess.C: Data generated with incompatible version of Loki."
              << endl;
         exit(-1);
      }

      // Get Nx and Ny from the first species.
      reader.pushSubDir(species_names[0]);
      tbox::Dimension dim(4);
      tbox::IntVector n_cells(dim);
      n_cells.getFromDatabase(reader, "N");
      nx = n_cells[0];
      ny = n_cells[1];
      reader.popSubDir();
   }

   // Finished with distribution functions so process the time histories and
   // fields.
   writeTimeHistories(prefix, isMaxwell, false, species_names);

   writeFields(prefix,
      nx,
      ny,
      nGhost,
      isMaxwell,
      false,
      species_names,
      plot_ke_vel_bdy_flux);

   if (coll) {
      writeTimeHistories(prefix, isMaxwell, true, species_names);
      writeFields(prefix,
         nx,
         ny,
         nGhost,
         isMaxwell,
         true,
         species_names,
         plot_ke_vel_bdy_flux);
   }

   MPI_Finalize();
   return 0;
}
