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
#include "External2DIC.H"
#include "Array.H"
#include "BoxOps.H"
#include "External2DICF.H"
#include "Loki_Utilities.H"
#include "hdf5.h"

namespace Loki {

bool
External2DIC::isType(
   const aString& a_name)
{
   if (a_name.matches("External 2D")) {
      return true;
   }
   return false;
}


External2DIC::External2DIC(
   ParmParse& a_pp,
   real a_vflowinitx,
   real a_vflowinity)
{
   // Set some hopefully sane input parameter defaults.
   m_parameters.resize(NUM_PARAMS);
   m_parameters(ALPHA)         = 1.0;
   m_parameters(BETA)          = 1.0;
   m_parameters(VX0)           = 0.0;
   m_parameters(VY0)           = 0.0;
   m_parameters(VFLOWINITX)    = a_vflowinitx;
   m_parameters(VFLOWINITY)    = a_vflowinity; 
   m_parameters(X_WAVE_NUMBER) = 0.0;
   m_parameters(Y_WAVE_NUMBER) = 0.0;
   m_parameters(FLOW_VEL_PHI)  = 0.0;

   // Now read what the user really wants.
   parseParameters(a_pp);
}


External2DIC::~External2DIC()
{
}


void
External2DIC::set(
   RealArray& a_u,
   const ProblemDomain& a_domain,
   const tbox::Box& a_grown_global_box,
   real a_time) const
{
   NULL_USE(a_time);
   tbox::Box local_box(BoxOps::getLocalBox(a_u));

   // Open 2D initial condition distribution file.
   hid_t file_id = H5Fopen(m_ext_2d_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file_id < 0) {
      OV_ABORT("Unable to open external 2D distribution file.");
   }

   // Open the dataset containing the 2D distribution.
   hid_t dataset_id;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   dataset_id = H5Dopen(file_id, "2D_dist");
#else
   dataset_id = H5Dopen(file_id, "2D dist", H5P_DEFAULT);
#endif
#else
   dataset_id = H5Dopen(file_id, "2D_dist");
#endif
   if (dataset_id > 0) {
      // Check that the dataset is the right size.
      int data_size =
         (a_grown_global_box.upper(0)-a_grown_global_box.lower(0)+1)*
         (a_grown_global_box.upper(1)-a_grown_global_box.lower(1)+1);
      hid_t dspace = H5Dget_space(dataset_id);
      if (dspace < 0) {
         OV_ABORT("Can not get dataspace for dataset \"2D dist\".");
      }
      hsize_t nsel = H5Sget_select_npoints(dspace);
      if (static_cast<int>(nsel) != data_size) {
         OV_ABORT("Distribution size does not match configuration space.");
      }

      // Now that everything checks out, read the dataset.
      double* data = new double [data_size];
      herr_t errf = H5Dread(dataset_id,
                            H5T_NATIVE_DOUBLE,
                            H5S_ALL,
                            H5S_ALL,
                            H5P_DEFAULT,
                            data);
      if (errf < 0) {
         OV_ABORT("Can not read dataset \"2D dist\".");
      }

      // Close the dataspace.
      errf = H5Sclose(dspace);
      if (errf < 0) {
         OV_ABORT("Can not close dataspace for dataset \"2D dist\".");
      }

      // Close the dataset containing the 2D distribution.
      errf = H5Dclose(dataset_id);
      if (errf < 0) {
         OV_ABORT("Can not close dataset \"2D dist\".");
      }

      // Now pass the initial condition data and input parameters to fortran
      // to set the distribution function.
      setIC_external_2D(BOX4D_TO_FORT(local_box),
         BOX4D_TO_FORT(a_grown_global_box),
         *a_u.getDataPointer(),
         *data,
         *m_parameters.getDataPointer(),
         PROBLEMDOMAIN_TO_FORT(a_domain),
         INTVECTOR4D_TO_FORT(a_domain.numberOfCells()));

      delete [] data;
   }
   else {
      OV_ABORT("Can not open dataset \"2D dist\".");
   }
}


void
External2DIC::printParameters() const
{
   // Print all the user input parameters.
   printF("\n  Using external 2D distribution:\n");
   printF("    external 2D file    = %s\n", m_ext_2d_file.c_str());
   printF("    alpha               = %e\n", m_parameters(ALPHA));
   printF("    beta                = %e\n", m_parameters(BETA));
   printF("    vx0                 = %e\n", m_parameters(VX0));
   printF("    vy0                 = %e\n", m_parameters(VY0));
   printF("    x wave number       = %e\n", m_parameters(X_WAVE_NUMBER));
   printF("    y wave number       = %e\n", m_parameters(Y_WAVE_NUMBER));
   printF("    flow velocity phase = %e\n", m_parameters(FLOW_VEL_PHI));
}


void
External2DIC::parseParameters(
   ParmParse& a_pp)
{
   // The file name is required.
   if (!a_pp.query("file_name", m_ext_2d_file)) {
      OV_ABORT("Must supply name of external 2D distribution file.");
   }

   // Currently none of these inputs is required so the defaults better be
   // sane.
   a_pp.query("alpha",         m_parameters(ALPHA));
   a_pp.query("beta",          m_parameters(BETA));
   a_pp.query("vx0",           m_parameters(VX0));
   a_pp.query("vy0",           m_parameters(VY0));
   a_pp.query("x_wave_number", m_parameters(X_WAVE_NUMBER));
   a_pp.query("y_wave_number", m_parameters(Y_WAVE_NUMBER));
   a_pp.query("phase",         m_parameters(FLOW_VEL_PHI));
}

} // end namespace Loki
