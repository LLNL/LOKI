/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "Loki_Utilities.H"
#include "RestartReader.H"
#include "FieldReader.H"
#include "TimeHistReader.H"
#include "hdf5.h"

using namespace std;
using namespace Loki;


int
GetSpeciesDistributionData(
   hid_t rootGID,
   const string& species,
   const int printLevel,
   double*& data_out)
{
   string dist_name(species);
   dist_name += "\\distribution.p0";
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dsid = H5Dopen(rootGID, dist_name.c_str());
#else
   hid_t dsid = H5Dopen(rootGID, dist_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dsid = H5Dopen(rootGID, dist_name.c_str());
#endif
   if (dsid < 0) {
      cout << "FAILED: Can not open distribution dataset." << endl;
      return -1;
   }

   // Now read the requested distribution from the dataset.
   herr_t readStatus =
      H5Dread(dsid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
   if (readStatus < 0) {
      cout << "FAILED: Could not read distribution data for species " << species
           << endl;
      return -2;
   }

   // Close the dataset.
   herr_t errf = H5Dclose(dsid);
   if (errf < 0) {
      cout << "Can not close timeseries dataset" << endl;
      return -3;
   }

   // Successful return
   return 0;
}


int
CheckFieldData(
   const stringstream& baseline_field_file,
   const stringstream& test_field_file,
   const string& field_tol_file_name)
{
   // Get a reader for the baseline fields.
   FieldReader baseline_reader(baseline_field_file.str());

   // Get a reader for the test fields.
   FieldReader test_reader(test_field_file.str());

   // Loop over all the fields specified in the field tolerance file and check
   // that the maximium relative difference between the baseline and test of
   // each is less than the specified relative tolerance.
   vector<double> field_data;
   vector<double> tfield_data;
   ifstream field_tol_file;
   field_tol_file.open(field_tol_file_name.c_str());
   while (true) {
      // Get the name and relative tolerance of this field.
      string field_name, field_tol_str;
      getline(field_tol_file, field_name);
      // We've run out of fields to check.
      if (field_tol_file.eof()) {
         break;
      }
      getline(field_tol_file, field_tol_str);
      double field_tol = atof(field_tol_str.c_str());

      // Allocate storge and then read this field's data from the baseline and
      // test field files.
      baseline_reader.readField(field_name, field_data);
      test_reader.readField(field_name, tfield_data);
      int baseline_field_size = static_cast<int>(field_data.size());
      int test_field_size = static_cast<int>(tfield_data.size());

      if (baseline_field_size != test_field_size) {
         cout << "FAILED. Test and baseline field sizes differ." << endl;
         return -1;
      }

      // Compare field data and check relative difference between baseline and
      // test.
      double max_rel_diff = 0.0;
      for (int i = 0; i < baseline_field_size; ++i) {
         double baseline = field_data[i];
         double test = tfield_data[i];
         double rel_diff;
         if (test == 0.0) {
            if (baseline == 0.0) {
               rel_diff = 0.0;
            }
            else {
               rel_diff = fabs((test-baseline)/baseline);
            }
         }
         else {
            rel_diff = fabs((test-baseline)/test);
         }
         if (rel_diff > max_rel_diff) {
            max_rel_diff = rel_diff;
         }
      }
      if (max_rel_diff > field_tol) {
         cout << "Maximum relative difference for field " << field_name << ": "
              << max_rel_diff << " exceeds tolerance " << field_tol << endl;
      }
   }

   // Close the field tolerance file.
   field_tol_file.close();

   return 0;
}


int
CheckTimeSeriesData(
   const stringstream& baseline_timeseries_file,
   const stringstream& test_timeseries_file,
   const string& ts_tol_file_name)
{
   // Get a reader for the baseline time histories.
   TimeHistReader baseline_reader(baseline_timeseries_file.str());

   // Get a reader for the test time histories.
   TimeHistReader test_reader(test_timeseries_file.str());

   // Loop over all the timeseries specified in the timeseries tolerance file
   // and check that the maximium relative difference between the baseline and
   // test of each is less than the specified relative tolerance.
   vector<double> ts_data;
   vector<double> tts_data;
   ifstream ts_tol_file;
   ts_tol_file.open(ts_tol_file_name.c_str());
   while (true) {
      // Get the name and relative tolerance of this time series.
      string ts_name, ts_tol_str;
      getline(ts_tol_file, ts_name);
      // We've run out of time series to check.
      if (ts_tol_file.eof()) {
         break;
      }
      getline(ts_tol_file, ts_tol_str);
      double ts_tol = atof(ts_tol_str.c_str());

      // Read this series' data from the baseline and test time history files.
      baseline_reader.readTimeHistory(ts_name, ts_data);
      test_reader.readTimeHistory(ts_name, tts_data);

      if (ts_data.size() != tts_data.size()) {
         cout << "FAILED. Test and baseline timeseries sizes differ." << endl;
         return -1;
      }

      // Compare timeseries data and check relative difference between baseline
      // and test.
      double max_rel_diff = 0.0;
      int ts_size = static_cast<int>(ts_data.size());
      for (int i = 0; i < ts_size; ++i) {
         double baseline = ts_data[i];
         double test = tts_data[i];
         double rel_diff;
         if (test == 0.0) {
            if (baseline == 0.0) {
               rel_diff = 0.0;
            }
            else {
               rel_diff = fabs((test-baseline)/baseline);
            }
         }
         else {
            rel_diff = fabs((test-baseline)/test);
         }
         if (rel_diff > max_rel_diff) {
            max_rel_diff = rel_diff;
         }
      }
      if (max_rel_diff > ts_tol) {
         cout << "Maximum relative difference for time series " << ts_name
              << ": " << max_rel_diff << " exceeds tolerance " << ts_tol
              << endl;
      }
   }

   // Close the time series tolerance file.
   ts_tol_file.close();

   return 0;
}


int
CheckDistributionData(
   const stringstream& baseline_dist_bulkdata_file,
   const stringstream& test_dist_bulkdata_file,
   const string& dist_tol_file_name,
   int** species_N,
   map<string, int>& species_names_to_index,
   const int printLevel,
   const int ng)
{
   // Open the baseline bulkdata file.
   hid_t baseline_bulkdata_file =
      H5Fopen(baseline_dist_bulkdata_file.str().c_str(),
         H5F_ACC_RDONLY,
         H5P_DEFAULT);
   if (baseline_bulkdata_file < 0) {
      cout << "FAILED: Can not open baseline bulkdata file." << endl;
      return 1;
   }

   // Get the root group from the baseline bulkdata file.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Gopen_vers) && H5Gopen_vers == 1
   hid_t baseline_rootGID = H5Gopen(baseline_bulkdata_file, "/root");
#else
   hid_t baseline_rootGID =
      H5Gopen(baseline_bulkdata_file, "/root", H5P_DEFAULT);
#endif
#else
   hid_t baseline_rootGID = H5Gopen(baseline_bulkdata_file, "/root");
#endif
   if (baseline_rootGID < 0) {
      cout << "FAILED: Can not open root group of baseline distribution bulkdata file."
           << endl;
      return 2;
   }

   // Open the test bulkdata file.
   hid_t test_bulkdata_file =
      H5Fopen(test_dist_bulkdata_file.str().c_str(),
         H5F_ACC_RDONLY,
         H5P_DEFAULT);
   if (test_bulkdata_file < 0) {
      cout << "FAILED: Can not open test bulkdata file." << endl;
      return 3;
   }

   // Get the root group from the test bulkdata file.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Gopen_vers) && H5Gopen_vers == 1
   hid_t test_rootGID = H5Gopen(test_bulkdata_file, "/root");
#else
   hid_t test_rootGID = H5Gopen(test_bulkdata_file, "/root", H5P_DEFAULT);
#endif
#else
   hid_t test_rootGID = H5Gopen(test_bulkdata_file, "/root");
#endif
   if (test_rootGID < 0) {
      cout << "FAILED: Can not open root group of test distribution bulkdata file."
           << endl;
      return 4;
   }

   // For each species read and compare the distributions from the baseline and
   // test bulkdata files.

   // Read the species of interest and their tolerances from the distribution
   // tolerance file and for each compare the baseline and test distributions.
   ifstream dist_tols;
   dist_tols.open(dist_tol_file_name.c_str());
   while (true) {
      // Get the name and relative tolerance of this species.
      string species_name, species_tol_str;
      getline(dist_tols, species_name);
      // We've run out of species to check.
      if (dist_tols.eof()) {
         break;
      }
      getline(dist_tols, species_tol_str);
      double species_tol = atof(species_tol_str.c_str());
      if (printLevel > 0) {
         cout << "Reading distribution for species: " << species_name << endl;
      }

      // Allocate storage for this species' distribution.
      const int* this_species_N =
         species_N[species_names_to_index[species_name]];
      int n1 = this_species_N[0]+2*ng;
      int n2 = this_species_N[1]+2*ng;
      int n3 = this_species_N[2]+2*ng;
      int n4 = this_species_N[3]+2*ng;
      int data_size = n1*n2*n3*n4;
      double *data_out = new double [data_size];
      double *tdata_out = new double [data_size];

      int baseline_dist_status =
         GetSpeciesDistributionData(baseline_rootGID,
            species_name,
            printLevel,
            data_out);
      int test_dist_status =
         GetSpeciesDistributionData(test_rootGID,
            species_name,
            printLevel,
            tdata_out);
      if (baseline_dist_status != 0 || test_dist_status != 0) {
         cout << "FAILED. Species distribution not read." << endl;
         return 5;
      }

      // Compare data
      double max_rel_diff = 0.0;
      for (int i4 = ng; i4 < this_species_N[3]+ng; ++i4) {
         for (int i3 = ng; i3 < this_species_N[2]+ng; ++i3) {
            for (int i2 = ng; i2 < this_species_N[1]+ng; ++i2) {
               for (int i1 = ng; i1 < this_species_N[0]+ng; ++i1) {
                  int idx = i1 + i2*n1 + i3*n1*n2 + i4*n1*n2*n3;
                  double baseline = data_out[idx];
                  double test = tdata_out[idx];
                  double rel_diff;
                  if (test == 0.0) {
                     if (baseline == 0.0) {
                        rel_diff = 0.0;
                     }
                     else {
                        rel_diff = fabs((test-baseline)/baseline);
                     }
                  }
                  else {
                     rel_diff = fabs((test-baseline)/test);
                  }
                  if (rel_diff > max_rel_diff) {
                     max_rel_diff = rel_diff;
                  }
               }
            }
         }
      }
      if (max_rel_diff > species_tol) {
         cout << "Maximum relative difference for species " << species_name
              << ": " << max_rel_diff << " exceeds tolerance " << species_tol
              << endl;
      }

      delete [] data_out;
      delete [] tdata_out;
   }

   // Close the baseline bulkdata file.
   herr_t status;
   status = H5Fclose(baseline_bulkdata_file);
   if (status < 0) {
      cout << "FAILED: Can not close baseline bulkdata file." << endl;
      return 6;
   }

   // Close the test bulkdata file.
   status = H5Fclose(test_bulkdata_file);
   if (status < 0) {
      cout << "FAILED: Can not close test bulkdata file." << endl;
      return 7;
   }

   return 0;
}


int
main(
   int argc,
   char* argv[])
{
   // Check that the user has supplied the input file name.
   if (argc != 2) {
      cout << "Run with one argument for the input file name." << endl;
      return 1;
   }

   MPI_Init(&argc, &argv);
   Loki_Utilities::initialize();

   // Read general info from the input file.
   ifstream input;
   input.open(argv[1]);

   string dist_tol_file, ts_tol_file, field_tol_file, test_dir_name,
      test_file_prefix;
   int fileIndex, do_coll;
   input >> dist_tol_file >> ts_tol_file >> field_tol_file >> test_dir_name
         >> test_file_prefix >> fileIndex >> do_coll;

   const int printLevel = 0;

   // Get the baseline meta and bulk data file names.
   stringstream baseline_dist_metadata_file;
   baseline_dist_metadata_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
      test_dir_name << "/" << test_file_prefix << "_dist_" << fileIndex <<
      ".hdf";
   stringstream baseline_dist_bulkdata_file;
   baseline_dist_bulkdata_file << baseline_dist_metadata_file.str() << ".g0";

   // Get the test meta and bulk data file names.
   stringstream test_dist_metadata_file;
   test_dist_metadata_file << test_file_prefix << "_dist_" << fileIndex <<
      ".hdf";
   stringstream test_dist_bulkdata_file;
   test_dist_bulkdata_file << test_dist_metadata_file.str() << ".g0";

   // Open the baseline and test metadata files.
   if (printLevel > 0) {
      cout << "Attempting to load files " << baseline_dist_metadata_file.str()
         << endl << " and " << endl << test_dist_metadata_file.str() << endl;
   }
   RestartReader* baseline_metadata_db =
      new RestartReader(baseline_dist_metadata_file.str(), 1, true);
   RestartReader* test_metadata_db =
      new RestartReader(test_dist_metadata_file.str(), 1, true);
   if (printLevel > 0) {
      cout << "Loaded meta files " << baseline_dist_metadata_file.str() << endl
         << " and " << endl << test_dist_metadata_file.str() << endl;
   }

   // Now read needed metadata from baseline and test metadata files and ensure
   // consistency,

   // Get the number of ghosts.
   int nGhost;
   baseline_metadata_db->readIntegerValue("nGhost", nGhost);
   int nGhostTest;
   test_metadata_db->readIntegerValue("nGhost", nGhostTest);
   if (nGhost != nGhostTest) {
      cout <<
         "FAILED: Baseline and test have different number of ghosts." << endl;
      return 4;
   }
   if (printLevel > 0) {
      cout << "Number of ghosts: " << nGhost << endl;
   }

   // Get the number of species.  For each species get its name and
   // computational domain.
   // Check number of species.
   int baseline_num_species;
   baseline_metadata_db->readIntegerValue("species_list_size",
      baseline_num_species);
   int test_num_species;
   test_metadata_db->readIntegerValue("species_list_size", test_num_species);
   if (baseline_num_species != test_num_species || baseline_num_species < 1) {
      cout <<
         "FAILED: Baseline and test have different number of species." << endl;
       return 5;
   }
   if (printLevel > 0) {
      cout << "Number of species: " << baseline_num_species << endl;
   }

   // Check species names.
   baseline_metadata_db->pushSubDir("species_list");
   test_metadata_db->pushSubDir("species_list");
   vector<string> species_names(baseline_num_species);
   map<string, int> species_names_to_index;
   for (int s = 1; s <= baseline_num_species; ++s) {
      string tmp1, tmp2;
      stringstream tag;
      tag << "species." << s;
      baseline_metadata_db->readString(tag.str(), tmp1);
      test_metadata_db->readString(tag.str(), tmp2);
      if (tmp1.compare(tmp2) != 0) {
         cout << "FAILED: Baseline and test species lists differ." << endl;
         return 8;
      }
      species_names[s-1] = tmp1;
      species_names_to_index[tmp1] = s-1;
      if (printLevel > 0) {
         cout << "Species " << s << ": " << tag.str() << ", name " << tmp1
            << endl;
      }
   }
   baseline_metadata_db->popSubDir();
   test_metadata_db->popSubDir();

   // Check the computational domains.
   int field_dims[2];
   int **species_N = new int*[baseline_num_species];
   for (int s = 1; s <= baseline_num_species; ++s) {
      // Get the size of the computational domain of this species.
      string& species_name = species_names[s-1];
      baseline_metadata_db->pushSubDir(species_name);
      test_metadata_db->pushSubDir(species_name);

      species_N[s-1] = new int[4];
      baseline_metadata_db->readIntegerArray("N", species_N[s-1], 4);
      if (s == 1) {
         field_dims[0] = species_N[s-1][0];
         field_dims[1] = species_N[s-1][1];
      }
      int Ntest[4];
      const int *N = species_N[s-1];
      test_metadata_db->readIntegerArray("N", Ntest, 4);
      if (N[0] != Ntest[0] || N[1] != Ntest[1] ||
          N[2] != Ntest[2] || N[3] != Ntest[3]) {
         cout <<
            "FAILED: Baseline and test have different computational domains." <<
            endl;
         return 11;
      }
      if (N[0] != field_dims[0] || N[1] != field_dims[1]) {
         cout << "FAILED: Inconsistent configuration space dims." << endl;
         return 12;
      }

      if (printLevel > 0) {
         cout << "size ";
         for (int i=0; i<4; ++i) {
            cout << N[i] << " ";
         }
         cout << endl;
      }
      baseline_metadata_db->popSubDir();
      test_metadata_db->popSubDir();
   }

   // Done with the metadata files.
   delete baseline_metadata_db;
   delete test_metadata_db;

   // Get the baseline field file name.
   stringstream baseline_field_file;
   baseline_field_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
      test_dir_name << "/" << test_file_prefix << "_fields.hdf";

   // Get the test field file names.
   stringstream test_field_file;
   test_field_file << test_file_prefix << "_fields.hdf";

   // Check the field data.
   int err = CheckFieldData(baseline_field_file,
      test_field_file,
      field_tol_file);
   if (err != 0) {
     return err;
   }

   // Get the baseline time series file names.
   stringstream baseline_timeseries_file;
   baseline_timeseries_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
      test_dir_name << "/" << test_file_prefix << "_timeSeries.hdf";

   // Get the test time series file names.
   stringstream test_timeseries_file;
   test_timeseries_file << test_file_prefix << "_timeSeries.hdf";

   // Check time series data.
   err = CheckTimeSeriesData(baseline_timeseries_file,
      test_timeseries_file,
      ts_tol_file);
   if (err != 0) {
      return err;
   }

   // If checking the collision diagnostics is requested check those fields and
   // time series.
   if (do_coll) {
      // Get the names of the collision tolerance files.
      string coll_ts_tol_file, coll_field_tol_file;
      input >> coll_ts_tol_file >> coll_field_tol_file;

      // Get the baseline collision diagnostics field file name.
      stringstream baseline_coll_field_file;
      baseline_coll_field_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
         test_dir_name << "/" << test_file_prefix << "_coll_fields.hdf";

      // Get the test collision diagnostics field file name.
      stringstream test_coll_field_file;
      test_coll_field_file << test_file_prefix << "_coll_fields.hdf";

      // Check the field data.
      int err = CheckFieldData(baseline_coll_field_file,
         test_coll_field_file,
         coll_field_tol_file);
      if (err != 0) {
        return err;
      }

      // Get the baseline collision diagnostics time series file name.
      stringstream baseline_coll_timeseries_file;
      baseline_coll_timeseries_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
         test_dir_name << "/" << test_file_prefix << "_coll_timeSeries.hdf";

      // Get the test collision diagnostics time series file name.
      stringstream test_coll_timeseries_file;
      test_coll_timeseries_file << test_file_prefix << "_coll_timeSeries.hdf";

      // Check time series data.
      err = CheckTimeSeriesData(baseline_coll_timeseries_file,
         test_coll_timeseries_file,
         coll_ts_tol_file);
      if (err != 0) {
         return err;
      }
   }
   input.close();

   // Check the distribution data.
   err = CheckDistributionData(baseline_dist_bulkdata_file,
      test_dist_bulkdata_file,
      dist_tol_file,
      species_N,
      species_names_to_index,
      printLevel,
      nGhost);
   if (err != 0) {
      return err;
   }

   // Delete allocated storage.
   for (int i = 0; i < baseline_num_species; ++i) {
      delete [] species_N[i];
   }
   delete [] species_N;

   MPI_Finalize();

   return 0;
}

