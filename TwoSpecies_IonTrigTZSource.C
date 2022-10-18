/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "TwoSpecies_IonTrigTZSource.H"
#include "TwoSpecies_IonTZSourceF.H"
#include "Loki_Utilities.H"

namespace Loki {

const string TwoSpecies_IonTrigTZSource::s_CLASS_NAME("TwoSpecies_IonTrigTZSource");

bool
TwoSpecies_IonTrigTZSource::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}

TwoSpecies_IonTrigTZSource::TwoSpecies_IonTrigTZSource(
   LokiInputParser& a_pp)
   : m_dparameters(NUM_DPARAMS)
{
   if (!a_pp.query("amp", m_dparameters[AMP])) {
      LOKI_ABORT("Must supply amp");
   }
   if (!a_pp.query("electron_mass", m_dparameters[EMASS])) {
      LOKI_ABORT("Must supply electron mass");
   }
   if (!a_pp.query("ion_mass", m_dparameters[IMASS])) {
      LOKI_ABORT("Must supply ion mass");
   }
}


TwoSpecies_IonTrigTZSource::~TwoSpecies_IonTrigTZSource()
{
}


void
TwoSpecies_IonTrigTZSource::set(
   ParallelArray& a_dist_func,
   const ProblemDomain& a_domain,
   double a_time,
   const ParallelArray& a_velocities) const
{
   FORT_SET_TWOION_TRIG_TZ_SOURCE(*a_dist_func.getData(),
      BOX4D_TO_FORT(a_dist_func.dataBox()),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      *a_velocities.getData(),
      m_dparameters[0]);
}


void
TwoSpecies_IonTrigTZSource::computeError(
   ParallelArray& a_tz_error,
   const ParallelArray& a_dist_func,
   const ProblemDomain& a_domain,
   double a_time,
   const ParallelArray& a_velocities) const
{
   FORT_COMPUTE_TWOION_TRIG_TZ_SOURCE_ERROR(*a_tz_error.getData(),
      *a_dist_func.getData(),
      BOX4D_TO_FORT(a_tz_error.dataBox()),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      *a_velocities.getData(),
      m_dparameters[0]);


}


void
TwoSpecies_IonTrigTZSource::printParameters() const
{
   Loki_Utilities::printF("  Using twilight zone source:\n");
   Loki_Utilities::printF("    amplitude     = %f\n", m_dparameters[AMP]);
   Loki_Utilities::printF("    electron mass = %f\n", m_dparameters[EMASS]);
   Loki_Utilities::printF("    ion mass      = %f\n", m_dparameters[IMASS]);
}

} // end namespace Loki
