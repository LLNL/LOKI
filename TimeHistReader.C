/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "TimeHistReader.H"

namespace Loki {


TimeHistReader::TimeHistReader(
   const string& a_name)
   : ReaderWriterBase()
{
   // Open this time history file and the root group.
   openFileAndRoot(a_name, m_file, m_root);
}


TimeHistReader::~TimeHistReader()
{
   // Close the file and root group.
   closeFileAndRoot(m_file, m_root);
}


void
TimeHistReader::readTimeHistory(
   const string& a_name,
   vector<double>& a_vals)
{
   readDoubleArray(a_name, m_root, a_vals);
}


void
TimeHistReader::readNumProbes(
   int& a_num_probes)
{
   readIntegerValue("numProbes", m_root, a_num_probes);
}


void
TimeHistReader::readNumTrackingParticles(
   int& a_num_tracking_particles)
{
   readIntegerValue("numTrackingParticles", m_root, a_num_tracking_particles);
}
} // end namespace Loki
