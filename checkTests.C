#include <iostream>
#include <sstream>
#include <vector>

#include "HDF5_DataBase.h"
#include "hdf5.h"

#define MAX_NAME 1024
#define NDIM 6  // Hard-coding this as 6 assumes 4D data.

using namespace std;

double* GetDataSet(const hid_t dsid, hsize_t *dims)
{
  hid_t filespace = H5Dget_space(dsid);

  herr_t status_n = H5Sget_simple_extent_dims(filespace, dims, NULL);
  hid_t memspace = H5Screate_simple(NDIM,dims,NULL);

  const int printLevel = 0;
  
  if (printLevel > 0)
    {
      for (int i=0; i<NDIM; ++i)
	cout << "dim[" << i << "] " << dims[i] << endl;
    }

  const long int dataSize = dims[0]*dims[1]*dims[2]*dims[3]*dims[4]*dims[5];
  if (dataSize < 1)
    return NULL;
  
  double data[dims[0]][dims[1]][dims[2]][dims[3]][dims[4]][dims[5]];
  /*
  // Why doesn't the following dynamic allocation work? It results in a segmentation fault in H5Dread().
  double ******data = new double*****[dims[0]];

  for (int i0=0; i0<dims[0]; ++i0)
    {
      data[i0] = new double****[dims[1]];
      for (int i1=0; i1<dims[1]; ++i1)
	{
	  data[i0][i1] = new double***[dims[2]];
	  for (int i2=0; i2<dims[2]; ++i2)
	    {
	      data[i0][i1][i2] = new double**[dims[3]];
	      for (int i3=0; i3<dims[3]; ++i3)
		{
		  data[i0][i1][i2][i3] = new double*[dims[4]];
		  for (int i4=0; i4<dims[4]; ++i4)
		    {
		      data[i0][i1][i2][i3][i4] = new double[dims[5]];
		    }
		}
	    }
	}
    }
  */
  
  herr_t readStatus = H5Dread(dsid, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, data);

  H5Dclose(dsid);
  H5Sclose(memspace);
  
  double *data_out = new double[dataSize];

  if (data_out == NULL)
    {
      cout << "FAILED. Data array not allocated." << endl;
      return NULL;
    }
  
  unsigned long int os = 0;
  for (int i0=0; i0<dims[0]; ++i0)
    {
      for (int i1=0; i1<dims[1]; ++i1)
	{
	  for (int i2=0; i2<dims[2]; ++i2)
	    {
	      for (int i3=0; i3<dims[3]; ++i3)
		{
		  for (int i4=0; i4<dims[4]; ++i4)
		    {
		      for (int i5=0; i5<dims[5]; ++i5)
			data_out[os + i5] = data[i0][i1][i2][i3][i4][i5];

		      os += dims[5];
		    }
		}
	    }
	}
    }
  
  return data_out;
}

int GetTimeSeriesData(const stringstream& filename, const int nTimeHistory, hsize_t *dims, double **data_out)
{
  char filename_char[filename.str().length() + 1];
  strcpy(filename_char, filename.str().c_str());
  
  hid_t file = H5Fopen(filename_char, H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t rootGID = H5Gopen(file,"/root", H5P_DEFAULT);  // Assumption: /root will always be root group under RootGroup.
  
  herr_t err;
  hsize_t nobj;
  char group_name[MAX_NAME];
  char memb_name[MAX_NAME];
  
  const int printLevel = 0;
  
  err = H5Gget_num_objs(rootGID, &nobj);
  ssize_t len = H5Iget_name(rootGID, group_name, MAX_NAME);

  if (printLevel > 0)
    cout << "Root group number: " << rootGID << ", root group name " << group_name << ", len " << len
	 << ", number of objects in root group: " << nobj << endl;

  if (strcmp("/root", group_name) != 0)
    {
      cout << "FAILED: H5 root group is not named /root" << endl;
      return 20;
    }

  /*
  bool dataFound = false;
  
  for (int i = 0; i < nobj; ++i)
    {
      len = H5Gget_objname_by_idx(rootGID, (hsize_t)i, memb_name, (size_t)MAX_NAME);
      if (strcmp("series_data", memb_name) == 0)
	{
	  if (dataFound)
	    return 1;
	}
    }
  
  if (!dataFound)
    return 1;
  */

  hid_t dsid = H5Dopen(rootGID, "series_data", H5P_DEFAULT);
  (*data_out) = GetDataSet(dsid, dims);

  if (printLevel > 0)
    {
      for (int i=0; i<NDIM; ++i)
	cout << "TS dim[" << i << "] " << dims[i] << endl;
    }
  
  if (dims[0] != 1 || dims[1] != 1 || dims[2] != 1 || dims[3] != 1 || nTimeHistory != dims[4]
      || dims[4] < 1 || dims[5] < 1)
    {
      cout << "FAILED. Time series data dimensions are invalid." << endl;
      return 19;
    }
  
  herr_t status = H5Fclose(file);

  return 0;
}

int GetFieldData(const stringstream& field_filename, const vector<string>& species_names, bool& Maxwell, hsize_t **fieldDims,
		 double ***data_out, bool **fieldKeFluxVy)
{
  char field_filename_char[field_filename.str().length() + 1];
  strcpy(field_filename_char, field_filename.str().c_str());
  
  hid_t field_file = H5Fopen(field_filename_char, H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t rootGID = H5Gopen(field_file,"/root", H5P_DEFAULT);  // Assumption: /root will always be root group under RootGroup.

  herr_t err;
  hsize_t nobj;
  char group_name[MAX_NAME];
  char memb_name[MAX_NAME];
  
  const int printLevel = 0;
  
  err = H5Gget_num_objs(rootGID, &nobj);
  ssize_t len = H5Iget_name(rootGID, group_name, MAX_NAME);

  if (printLevel > 0)
    cout << "Root group number: " << rootGID << ", root group name " << group_name << ", len " << len
	 << ", number of objects in root group: " << nobj << endl;

  if (strcmp("/root", group_name) != 0)
    {
      cout << "FAILED: H5 root group is not named /root" << endl;
      return -1;
    }

  Maxwell = false;
  for (int i = 0; i < nobj; ++i)
    {
      len = H5Gget_objname_by_idx(rootGID, (hsize_t)i, memb_name, (size_t)MAX_NAME);
      if (strcmp("BZ", memb_name) == 0)
	Maxwell = true;

      /*
      int otype =  H5Gget_objtype_by_idx(rootGID, (size_t)i );
      if (printLevel > 10)
	cout << "Member " << i << " name: " << memb_name << endl; // ", type " << otype << ", dataset type " << H5G_DATASET << endl;
      */
    }

  vector<string> fields;

  fields.push_back("EX");

  // Note that EY is omitted.
  
  if (Maxwell)
    {
      fields.push_back("EZ");
      fields.push_back("BX");
      fields.push_back("BY");
      fields.push_back("BZ");

      for (vector<string>::const_iterator it = species_names.begin(); it != species_names.end(); ++it)
	fields.push_back(*it + " VZ");
    }
  
  for (vector<string>::const_iterator it = species_names.begin(); it != species_names.end(); ++it)
    {
      fields.push_back(*it + " ke flux vx hi");
      fields.push_back(*it + " ke flux vx lo");
      fields.push_back(*it + " ke flux vy hi");
      fields.push_back(*it + " ke flux vy lo");
    }

  *fieldDims = new hsize_t[fields.size() * NDIM];
  hsize_t dims[NDIM];
  
  *data_out = new double*[fields.size()];

  *fieldKeFluxVy = new bool[fields.size()];
  
  int cnt = 0;
  for (vector<string>::const_iterator it = fields.begin(); it != fields.end(); ++it, ++cnt)
    {
      // cout << "Field " << cnt << " is " << *it << endl;
      
      hid_t dsid = H5Dopen(rootGID, it->c_str(), H5P_DEFAULT);

      (*data_out)[cnt] = GetDataSet(dsid, dims);

      if ((*data_out)[cnt] == NULL)
	return -2;
      
      for (int i=0; i<NDIM; ++i)
	(*fieldDims)[(cnt*NDIM) + i] = dims[i];

      (*fieldKeFluxVy)[cnt] = false;
    }

  herr_t status = H5Fclose(field_file);

  int os = fields.size() - (4*species_names.size());
  for (cnt=0; cnt < species_names.size(); ++cnt)
    {
      (*fieldKeFluxVy)[os + (4*cnt) + 2] = true;
      (*fieldKeFluxVy)[os + (4*cnt) + 3] = true;
    }
  
  return fields.size();
}

int GetSpeciesDistributionData(const stringstream& bulkdata_filename, const string& species, const int printLevel,
			       hsize_t *dims, double **data_out)
{
  char bulkdata_filename_char[bulkdata_filename.str().length() + 1];
  strcpy(bulkdata_filename_char, bulkdata_filename.str().c_str());
  
  hid_t bulkdata_file = H5Fopen(bulkdata_filename_char, H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t rootGID = H5Gopen(bulkdata_file,"/root", H5P_DEFAULT);  // Assumption: /root will always be root group under RootGroup.

  herr_t err;
  hsize_t nobj;
  char group_name[MAX_NAME];
  char memb_name[MAX_NAME];
  
  err = H5Gget_num_objs(rootGID, &nobj);
  ssize_t len = H5Iget_name(rootGID, group_name, MAX_NAME);

  if (printLevel > 0)
    cout << "Root group number: " << rootGID << ", root group name " << group_name << ", len " << len
	 << ", number of objects in root group: " << nobj << endl;

  if (strcmp("/root", group_name) != 0)
    {
      cout << "FAILED: H5 root group is not named /root" << endl;
      return 8;
    }

  /*
  for (int i = 0; i < nobj; ++i)
    {
      len = H5Gget_objname_by_idx(rootGID, (hsize_t)i, memb_name, (size_t)MAX_NAME);
      int otype =  H5Gget_objtype_by_idx(rootGID, (size_t)i );
      if (printLevel > 10)
	cout << "Member " << i << " name: " << memb_name << ", type " << otype << ", dataset type " << H5G_DATASET << endl;
    }
  */

  bool distFound = false;
  hid_t dsid = -1;
      
  for (int i = 0; i < nobj; ++i)
    {
      len = H5Gget_objname_by_idx(rootGID, (hsize_t)i, memb_name, (size_t)MAX_NAME);
      string truncated_memb_name(memb_name, species.length());
      if (species.compare(0, species.length(), truncated_memb_name) == 0)
	{
	  if (printLevel > 0)
	    cout << "Found distribution: " << memb_name << " for species " << species << endl;
	      
	  if (distFound)
	    {
	      cout << "FAILED. Non-unique object name." << endl;
	      return 10;
	    }
	      
	  distFound = true;

	  int otype =  H5Gget_objtype_by_idx(rootGID, (size_t)i);
	  dsid = H5Dopen(rootGID, memb_name, H5P_DEFAULT);

	  if (otype != H5G_DATASET)
	    {
	      cout << "FAILED. Distribution must be of type DataSet." << endl;
	      return 9;
	    }
	}
    }
      
  if (!distFound)
    {
      cout << "FAILED. Species distribution not found." << endl;
      return 12;
    }

  *data_out = GetDataSet(dsid, dims);

  herr_t status = H5Fclose(bulkdata_file);

  if (*data_out == NULL)
    {
      cout << "FAILED. Data not allocated." << endl;
      return 13;
    }
    
  return 0;  // Successful return
}

int CheckFieldData(const stringstream& baseline_field_file, const stringstream& test_field_file, const vector<string>& species_names,
		   const double absTol, const double tol, const double keFluxVyTol, bool& Maxwell)
{
  hsize_t *fieldDims = NULL;
  hsize_t *tfieldDims = NULL;
  
  double **fieldData = NULL;
  double **tfieldData = NULL;

  Maxwell = false;
  bool tMaxwell = false;

  bool *fieldKeFluxVy = NULL;
  bool *tFieldKeFluxVy = NULL;
  
  const int numFields = GetFieldData(baseline_field_file, species_names, Maxwell, &fieldDims, &fieldData, &fieldKeFluxVy);
  const int tnumFields = GetFieldData(test_field_file, species_names, tMaxwell, &tfieldDims, &tfieldData, &tFieldKeFluxVy);
  
  if (numFields < 1 || tnumFields != numFields || Maxwell != tMaxwell)
    {
      cout << "FAILED. Could not read field data." << endl;
      return 17;
    }

  // Compare field data
  {
    for (int i=0; i<numFields; ++i)
      {
	double maxDiff = 0.0;
	double maxRelDiff = 0.0;

	for (int j=0; j<6; ++j)
	  {
	    if (fieldDims[(i*NDIM)+j] != tfieldDims[(i*NDIM)+j])
	      {
		cout << "FAILED. Field dimensions are different in baseline and test files." << endl;
		return 18;
	      }
	  }

	hsize_t *dims = fieldDims + (i*NDIM);	
	
	unsigned long int os = 0;
	double maxVal = 0.0;
	for (int i0=0; i0<dims[0]; ++i0)
	  {
	    for (int i1=0; i1<dims[1]; ++i1)
	      {
		for (int i2=0; i2<dims[2]; ++i2)
		  {
		    for (int i3=0; i3<dims[3]; ++i3)
		      {
			for (int i4=0; i4<dims[4]; ++i4)
			  {
			    for (int i5=0; i5<dims[5]; ++i5)
			      {
				maxVal = std::max(maxVal, fabs(fieldData[i][os + i5]));
			      }
			    
			    os += dims[5];
			  }
		      }
		  }
	      }
	  }

	const double atol = absTol * maxVal;
	
	os = 0;
	for (int i0=0; i0<dims[0]; ++i0)
	  {
	    for (int i1=0; i1<dims[1]; ++i1)
	      {
		for (int i2=0; i2<dims[2]; ++i2)
		  {
		    for (int i3=0; i3<dims[3]; ++i3)
		      {
			for (int i4=0; i4<dims[4]; ++i4)
			  {
			    for (int i5=0; i5<dims[5]; ++i5)
			      {
				const double diff = fabs(tfieldData[i][os + i5] - fieldData[i][os + i5]);
				maxDiff = std::max(maxDiff, diff);

				const double relDiff = diff / std::max(fabs(tfieldData[i][os + i5]), fabs(fieldData[i][os + i5]));

				if (diff > atol)
				  maxRelDiff = std::max(maxRelDiff, relDiff);
			      }
			    
			    os += dims[5];
			  }
		      }
		  }
	      }
	  }

	const double tol_i = fieldKeFluxVy[i] ? keFluxVyTol : tol;
	
	if (maxRelDiff > tol_i)
	  cout << "Maximum relative difference for field " << i << ": " << maxRelDiff << " exceeds tolerance " << tol_i << endl;
	
	delete [] fieldData[i];
	delete [] tfieldData[i];
      }

    delete [] fieldDims;
    delete [] tfieldDims;

    delete [] fieldData;
    delete [] tfieldData;

    delete [] fieldKeFluxVy;
    delete [] tFieldKeFluxVy;
  }
  
  return 0;
}

int CheckDistributionData(const stringstream& baseline_dist_bulkdata_file, const stringstream& test_dist_bulkdata_file, const string& species,
			  const int *species_N, const int printLevel, const int ng, const double absTol, const double tol)
{
  hsize_t dims[NDIM];
  double *data_out = NULL;

  hsize_t tdims[NDIM];
  double *tdata_out = NULL;

  if (GetSpeciesDistributionData(baseline_dist_bulkdata_file, species, printLevel, dims, &data_out) != 0 ||
      GetSpeciesDistributionData(test_dist_bulkdata_file, species, printLevel, tdims, &tdata_out) != 0)
    {
      cout << "FAILED. Species distribution not read." << endl;
      return 14;
    }

  if (printLevel > 100)
    {
      unsigned long int os = 0;
      for (int i2=0; i2<dims[2]; ++i2)
	for (int i3=0; i3<dims[3]; ++i3)
	  for (int i4=0; i4<dims[4]; ++i4)
	    {
	      for (int i5=0; i5<dims[5]; ++i5)
		cout << data_out[os + i5] << endl;
	      // cout << data_out[0][0][i2][i3][i4][i5] << endl;
		  
	      os += dims[5];
	    }
    }
      
  for (int i=0; i<4; ++i)
    {
      if (species_N[i] + (2*ng) != dims[NDIM-1-i])
	{
	  cout << "FAILED. Dimensions or ghost counts are incorrect." << endl;
	  return 13;
	}
    }

  for (int i=0; i<NDIM; ++i)
    if (dims[i] != tdims[i])
      {
	cout << "FAILED. Dimensions do not agree." << endl;
	return 16;
      }

  double maxVal = 0.0;
  unsigned long int os;
  for (int i2=ng; i2<dims[2]-ng; ++i2)
    for (int i3=ng; i3<dims[3]-ng; ++i3)
      for (int i4=ng; i4<dims[4]-ng; ++i4)
	{
	  for (int i5=ng; i5<dims[5]-ng; ++i5)
            {
              os = i2*(dims[3]*dims[4]*dims[5]) + i3*(dims[4]*dims[5]) + i4*dims[5];
	      maxVal = std::max(maxVal, fabs(data_out[os + i5]));
            }
	}

  const double atol = absTol * maxVal;
  
  // Compare data
  double maxDiff = 0.0;
  double maxRelDiff = 0.0;
  for (int i2=ng; i2<dims[2]-ng; ++i2)
    for (int i3=ng; i3<dims[3]-ng; ++i3)
      for (int i4=ng; i4<dims[4]-ng; ++i4)
	{
	  for (int i5=ng; i5<dims[5]-ng; ++i5)
	    {
              os = i2*(dims[3]*dims[4]*dims[5]) + i3*(dims[4]*dims[5]) + i4*dims[5];
	      const double diff = fabs(tdata_out[os + i5] - data_out[os + i5]);
	      maxDiff = std::max(maxDiff, diff);
	      
	      //if (diff / std::max(fabs(tdata_out[os + i5]), fabs(data_out[os + i5])) > 1.0)
	      //cout << tdata_out[os + i5] << " " << data_out[os + i5] << endl;
	      
	      if (diff > atol)
		maxRelDiff = std::max(maxRelDiff, diff / std::max(fabs(tdata_out[os + i5]), fabs(data_out[os + i5])));
	    }
	      
	}

  if (printLevel > 0)
    cout << "Species " << species << " maximum absolute difference on non-ghost nodes: " << maxDiff << endl
	 << "maximum relative difference on non-ghost nodes: " << maxRelDiff << endl;

  if (maxRelDiff > tol)
    cout << "Maximum relative difference for species " << species << " distribution on non-ghost nodes: " << maxRelDiff
	 << " exceeds tolerance " << tol << endl;
      
  delete [] data_out;
  delete [] tdata_out;

  return 0;
}

int CheckTimeSeriesData(const stringstream& baseline_file, const stringstream& test_file, const bool Maxwell,
			const vector<string>& species, const int nspecies,
			const int nprobes, const int ntracking, const char *tstolFile)
{
  double *data = NULL;
  double *tdata = NULL;

  hsize_t dims[NDIM];
  hsize_t tdims[NDIM];

  const int nTimeHistory = Maxwell ? 12 + (nprobes * (6 + nspecies)) + (4*ntracking) + (5*nspecies) :
    5 + (2*nprobes) + (4*ntracking) + (5*nspecies);
  
  int err = GetTimeSeriesData(baseline_file, nTimeHistory, dims, &data);
  if (err != 0)
    return err;

  err = GetTimeSeriesData(test_file, nTimeHistory, tdims, &tdata);
  if (err != 0)
    return err;

  if (dims[5] != tdims[5])
    {
      cout << "FAILED. Baseline and test dimensions differ in time series file." << endl;
      return 21;
    }

  vector<string> q;

  int qId[dims[4]];

  for (int i=0; i<dims[4]; ++i)
    qId[i] = -1;

  unsigned long int os = 0;

  if (Maxwell)
    {
      q.push_back("maxEnorm");
      q.push_back("intEnorm");
      q.push_back("maxEx");
      q.push_back("maxEy");
      q.push_back("maxEz");
      q.push_back("intEnorm2");
      q.push_back("maxBnorm");
      q.push_back("intBnorm");
      q.push_back("maxBx");
      q.push_back("maxBy");
      q.push_back("maxBz");
      q.push_back("intBnorm2");

      q.push_back("probeEx");
      q.push_back("probeEy");
      q.push_back("probeEz");
      q.push_back("probeBx");
      q.push_back("probeBy");
      q.push_back("probeBz");
      
      for (int i=0; i<nspecies; ++i)
	q.push_back(species[i] + ".probeVz");

      for (int i=0; i<12; ++i)
	qId[i] = i;

      os = 12;
      for (int i=0; i<nprobes; ++i)
	{
	  for (int j=0; j<6; ++j)
	    qId[os+j] = 12+j;

	  os += 6;

	  for (int k=0; k<nspecies; ++k, ++os)
	    qId[os+k] = 18+k;  // vz at the probe, for each kinetic species
	}
    }
  else  // Poisson
    {
      q.push_back("maxEnorm");
      q.push_back("intEnorm");
      q.push_back("maxEx");
      q.push_back("maxEy");
      q.push_back("intEnorm2");

      q.push_back("probeEx");
      q.push_back("probeEy");

      for (int i=0; i<5; ++i)
	qId[i] = i;

      os = 5;
      for (int i=0; i<nprobes; ++i)
	{
	  for (int j=0; j<2; ++j)
	    qId[os+j] = 5+j;

	  os += 2;
	}
    }

  // Quantities in common to VP and VM
  q.push_back("particleX");
  q.push_back("particleY");
  q.push_back("particleVx");
  q.push_back("particleVy");
      
  for (int i=0; i<nspecies; ++i)
    {
      q.push_back(species[i] + ".intKE");
      q.push_back(species[i] + ".KEfluxXlo");
      q.push_back(species[i] + ".KEfluxXhi");
      q.push_back(species[i] + ".KEfluxYlo");
      q.push_back(species[i] + ".KEfluxYhi");
    }

  for (int i=0; i<ntracking; ++i)
    {
      for (int j=0; j<4; ++j)
	{
	  if (Maxwell)
	    qId[os+j] = 18 + nspecies + j;
	  else
	    qId[os+j] = 7 + j;
	}

      os += 4;
    }

  for (int i=0; i<nspecies; ++i)
    {
      for (int j=0; j<5; ++j)
	{
	  if (Maxwell)
	    qId[os+j] = 22 + nspecies + (5*i) + j;
	  else
	    qId[os+j] = 11 + (5*i) + j;
	}
      
      os += 5;
    }

  if (os != dims[4])
    {
      cout << "BUG in qId." << endl;
      return 23;
    }
  
  for (int i=0; i<dims[4]; ++i)
    {
      if (qId[i] < 0)
	{
	  cout << "BUG in qId." << endl;
	  return 23;
	}
    }

  const int ntol = q.size();
  double tol[ntol];

  string quantity;
  ifstream tolfile;

  tolfile.open(tstolFile);

  for (int i=0; i<ntol; ++i)
    {
      tolfile >> quantity;
      // cout << quantity << " == " << q[i] << endl;
      if (strcmp(quantity.c_str(), q[i].c_str()) != 0)
	{
	  cout << "FAILED. Invalid timeseries tolerance file." << endl;
	  return 24;
	}

      tolfile >> tol[i];
    }
  
  tolfile.close();
  
  bool diff[dims[4]];
  for (int i4=0; i4<dims[4]; ++i4)
    diff[i4] = true;

  os = 0;
  if (!Maxwell)  // VP case
    {
      // Do not diff the following
      diff[3] = false;  // max Ey

      os = 5;
      
      for (int i=0; i<nprobes; ++i)
	{
	  diff[os + 1] = false;  // Ey for each probe
	  os += 2;
	}

      os += 4 * ntracking;
      
      for (int i=0; i<nspecies; ++i)
	{
	  diff[os + 3] = false;  // KE flux through y low physical boundary for each species
	  diff[os + 4] = false;  // KE flux through y high physical boundary for each species
	  
	  os += 5;
	}
    }

  os = 0;
  int qnt = 0;
  for (int i4=0; i4<dims[4]; ++i4)  // Loop over timeseries quantities.
    {
      if (diff[i4])
	{
	  double maxDiff = 0.0;
	  double maxRelDiff = 0.0;
	  
	  for (int i5=1; i5<dims[5]; ++i5)  // Loop over timesteps, skipping t=0. 
	    {
	      const double diff = fabs(data[os + i5] - tdata[os + i5]);
	      maxDiff = std::max(maxDiff, diff);

	      if (diff > 0.0)
		maxRelDiff = std::max(maxRelDiff, diff / std::max(fabs(data[os + i5]), fabs(tdata[os + i5])));
	    }

	  if (maxRelDiff > tol[qId[i4]])
	    cout << "Maximum relative difference for timeseries quantity " << q[qId[i4]] << ": " << maxRelDiff
		 << " exceeds tolerance " << tol[qId[i4]] << endl;
	}
      
      os += dims[5];
    }

  delete [] data;
  delete [] tdata;

  return 0;
}

int main(int argc, char* argv[])
{
  if (argc != 2)
    {
      cout << "Run with one argument for the input file name." << endl;

      //cout << "Run with arguments: absolute_tolerance distribution_tolerance field_tolerance "
      //"species_kefluxvy_field_tolerance tstolFile filename_root index nprobes ntracking" << endl;
      
      return 1;
    }

  /*
  const int argId_absTol = 1;
  const int argId_distTol = 2;
  const int argId_fieldTol = 3;
  const int argId_fieldKeFluxVyTol = 4;
  const int argId_tstolFile = 5;
  const int argId_file = 6;
  const int argId_index = 7;
  const int argId_nprobes = 8;
  const int argId_ntracking = 9;

  const int nprobes = atoi(argv[argId_nprobes]);
  const int ntracking = atoi(argv[argId_ntracking]);
  
  const double atol = atof(argv[argId_absTol]);

  const double dtol = atof(argv[argId_distTol]);
  const double ftol = atof(argv[argId_fieldTol]);
  const double fftol = atof(argv[argId_fieldKeFluxVyTol]);
  char* tstolFile = argv[argId_tstolFile];
  */

  ifstream input;

  input.open(argv[1]);

  int fileIndex, nprobes, ntracking;
  double atol, dtol, ftol, fftol;
  string tstolFile, testDir, testFile;
  input >> atol >> dtol >> ftol >> fftol >> tstolFile >> testDir >> testFile >> fileIndex >> nprobes >> ntracking;
  
  input.close();
     
  if (nprobes < 0 || ntracking < 0)
    {
      cout << "Number of probes and number of tracking particles must be non-negative" << endl;
      return 1;
    }
  
  if (dtol < 0.0 || ftol < 0.0)
    {
      cout << "Tolerance must be non-negative" << endl;
      return 1;
    }
  
  MPI_Init(&argc, &argv);

  const int printLevel = 0;
  
  // Get the baseline meta and bulk data file names.
  std::stringstream baseline_dist_metadata_file;
  baseline_dist_metadata_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
    testDir << "/" << testFile << "_dist_" << fileIndex << ".hdf";
  std::stringstream baseline_dist_bulkdata_file;
  baseline_dist_bulkdata_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
    testDir << "/" << testFile << "_dist_" << fileIndex << ".hdf.g0";

  // Get the baseline and test field file names.
  std::stringstream baseline_field_file;
  baseline_field_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
    testDir << "/" << testFile << "_fields.hdf";

  std::stringstream test_field_file;
  test_field_file << testFile << "_fields.hdf";

  // Get the baseline and test time series file names.
  std::stringstream baseline_timeseries_file;
  baseline_timeseries_file << "/usr/gapps/valhalla/LOKI/BASELINES/" <<
    testDir << "/" << testFile << "_timeSeries.hdf";

  std::stringstream test_timeseries_file;
  test_timeseries_file << testFile << "_timeSeries.hdf";
  
  // Get the test meta and bulk data file names.
  std::stringstream test_dist_metadata_file;
  test_dist_metadata_file << testFile << "_dist_" << fileIndex << ".hdf";
  std::stringstream test_dist_bulkdata_file;
  test_dist_bulkdata_file << testFile << "_dist_" << fileIndex << ".hdf.g0";

  if (printLevel > 0)
    cout << "Attempting to load files " << baseline_dist_metadata_file.str() << endl
	 << " and " << endl << test_dist_metadata_file.str() << endl;
   
  // Open the baseline and test metadata files.
  HDF_DataBase baseline_metadata_db;
  if (baseline_metadata_db.mount(baseline_dist_metadata_file.str(), "R") != 0) {
    std::cout << "FAILED: Can not open baseline metadata file." << std::endl;
    return 1;
  }
  HDF_DataBase test_metadata_db;
  if (test_metadata_db.mount(test_dist_metadata_file.str(), "R") != 0) {
    std::cout << "FAILED: Can not open test metadata file." << std::endl;
    return 1;
  }

  if (printLevel > 0)
    cout << "Loaded meta files " << baseline_dist_metadata_file.str() << endl
	 << " and " << endl << test_dist_metadata_file.str() << endl;

  // Get the number of ghosts.
  int nGhost;
  baseline_metadata_db.get(nGhost, "nGhost");
  int nGhostTest;
  test_metadata_db.get(nGhostTest, "nGhost");
  if (nGhost != nGhostTest) {
    std::cout <<
      "FAILED: Baseline and test have different number of ghosts." <<
      std::endl;
  }

  if (printLevel > 0)
    cout << "Number of ghosts: " << nGhost << endl;
   
  // Get the number of species and their names.
  int baseline_num_species;
  baseline_metadata_db.get(baseline_num_species, "species_list_size");
  int test_num_species;
  test_metadata_db.get(test_num_species, "species_list_size");
  if (baseline_num_species != test_num_species || baseline_num_species < 1) {
    std::cout <<
      "FAILED: Baseline and test have different number of species." <<
      std::endl;
    return 2;
  }

  if (printLevel > 0)
    cout << "Number of species: " << baseline_num_species << endl;

  HDF_DataBase baseline_metadata_sub_db;
  if (baseline_metadata_db.locate(baseline_metadata_sub_db, "species_list") != 0) {
    std::cout << "FAILED: Baseline missing species list." << std::endl;
    return 3;
  }
  HDF_DataBase test_metadata_sub_db;
  if (test_metadata_db.locate(test_metadata_sub_db, "species_list") != 0) {
    std::cout << "FAILED: Test missing species list." << std::endl;
    return 3;
  }
  
  std::vector<string> species_names;
  int **species_N = new int*[baseline_num_species];
  for (int s = 1; s <= baseline_num_species; ++s) {
    aString tmp1, tmp2;
    std::stringstream tag;
    tag << "species." << s;
    baseline_metadata_sub_db.get(tmp1, (aString)tag.str());
    test_metadata_sub_db.get(tmp2, (aString)tag.str());
    if (tmp1.compare(tmp2) != 0) {
      std::cout <<
	"FAILED: Baseline and test species lists differ." << std::endl;
      return 4;
    }
    species_names.push_back(tmp1);

    // Get the size of the computational domain.
    if (printLevel > 0)
      cout << "Species " << s << ": " << tag.str() << ", name " << tmp1 << endl;

    HDF_DataBase baseline_metadata_species_db;
    if (baseline_metadata_db.locate(baseline_metadata_species_db, tmp1) != 0) {
      std::cout << "FAILED: Baseline missing species entry." << std::endl;
      return 5;
    }
    HDF_DataBase test_metadata_species_db;
    if (test_metadata_db.locate(test_metadata_species_db, tmp1) != 0) {
      std::cout << "FAILED: Test missing species entry." << std::endl;
      return 5;
    }

    species_N[s-1] = new int[4];
    baseline_metadata_species_db.get(species_N[s-1], "N", 4);
    int Ntest[4];
    const int *N = species_N[s-1];
    test_metadata_species_db.get(Ntest, "N", 4);
    if (N[0] != Ntest[0] || N[1] != Ntest[1] ||
	N[2] != Ntest[2] || N[3] != Ntest[3]) {
      std::cout << "FAILED: Baseline and test have different computational domains." << std::endl;
    }

    if (printLevel > 0)
      {
	cout << "size ";
	for (int i=0; i<4; ++i)
	  cout << N[i] << " ";
	cout << endl;
      }
  }

  // Done with the metadata files.
  baseline_metadata_db.unmount();
  test_metadata_db.unmount();

  {
    // TODO: it seems these DataBase files are not used.
    // Now open the baseline and test bulkdata files.
    HDF_DataBase baseline_bulkdata_db;
    if (baseline_bulkdata_db.mount(baseline_dist_bulkdata_file.str(), "R") != 0) {
      std::cout << "FAILED: Can not open baseline bulkdata file." << std::endl;
      return 1;
    }

    HDF_DataBase test_bulkdata_db;
    if (test_bulkdata_db.mount(test_dist_bulkdata_file.str(), "R") != 0) {
      std::cout << "FAILED: Can not open test bulkdata file." << std::endl;
      return 1;
    }

    // Done with the metadata files.
    baseline_bulkdata_db.unmount();
    test_bulkdata_db.unmount();
  }
  
  bool Maxwell = false;
  {
    int err = CheckFieldData(baseline_field_file, test_field_file, species_names, atol, ftol, fftol, Maxwell);

    if (err != 0)
      return err;
  }
  
  // Read the distributions from each bulkdata file, for each species, and compare baseline and test.
  
  int scnt = 0;
  for (std::vector<string>::const_iterator it = species_names.begin(); it != species_names.end(); ++it, ++scnt)
    {
      if (printLevel > 0)
	cout << "Reading distribution for species: " << *it << endl;

      int err = CheckDistributionData(baseline_dist_bulkdata_file, test_dist_bulkdata_file, *it, species_N[scnt], printLevel, nGhost, atol, dtol);
      if (err != 0)
	return err;
    }

  // Check time series data.
  {
    int err = CheckTimeSeriesData(baseline_timeseries_file, test_timeseries_file, Maxwell, species_names, baseline_num_species, nprobes, ntracking,
				  tstolFile.c_str());
    if (err != 0)
      return err;
  }
  
  MPI_Finalize();

  return 0;
}

