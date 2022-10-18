/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ElectronTrigTZSource.H"
#include "ElectronTZSourceF.H"
#include "Loki_Utilities.H"

namespace Loki {

const string ElectronTrigTZSource::s_CLASS_NAME("ElectronTrigTZSource");

bool
ElectronTrigTZSource::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}

ElectronTrigTZSource::ElectronTrigTZSource(
   LokiInputParser& a_pp)
   : m_dparameters(NUM_DPARAMS)
{
   if (!a_pp.query("amp", m_dparameters[AMP])) {
      LOKI_ABORT("Must supply amp");
   }
}


ElectronTrigTZSource::~ElectronTrigTZSource()
{
}


void
ElectronTrigTZSource::set(
   ParallelArray& a_dist_func,
   const ProblemDomain& a_domain,
   double a_time,
   const ParallelArray& a_velocities) const
{
   FORT_SET_ELECTRON_TRIG_TZ_SOURCE(*a_dist_func.getData(),
      BOX4D_TO_FORT(a_dist_func.dataBox()),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      *a_velocities.getData(),
      m_dparameters[0]);
}


void
ElectronTrigTZSource::computeError(
   ParallelArray& a_tz_error,
   const ParallelArray& a_dist_func,
   const ProblemDomain& a_domain,
   double a_time,
   const ParallelArray& a_velocities) const
{
   FORT_COMPUTE_ELECTRON_TRIG_TZ_SOURCE_ERROR(*a_tz_error.getData(),
      *a_dist_func.getData(),
      BOX4D_TO_FORT(a_tz_error.dataBox()),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      *a_velocities.getData(),
      m_dparameters[0]);
}


void
ElectronTrigTZSource::printParameters() const
{
   Loki_Utilities::printF("  Using twilight zone source:\n");
   Loki_Utilities::printF("    amplitude     = %f\n", m_dparameters[AMP]);
}

} // end namespace Loki
