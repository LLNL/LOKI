/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "LokiInputParser.H"
#include "ParallelArray.H"
#include "ProblemDomain.H"
#include "tbox/Pointer.H"
#include <string>
#include <vector>

using namespace std;

namespace Loki {

class KineticSpeciesPtrVect;

/**
 * Abstract base class that defines the interface of the system we are solving
 * either Vlasov-Poisson or Vlasov-Maxwell.
 */
class System
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp The input database.
    * @param[in] a_spatial_solution_order Spatial order of solution.
    * @param[in] a_temporal_solution_order Temporal order of solution.
    */
   System(
      LokiInputParser& a_pp,
      int a_spatial_solution_order,
      int a_temporal_solution_order);

   /**
    * @brief Destructor
    */
   virtual
   ~System();

   /**
    * @brief Computes and returns a stable time step estimate.
    */
   virtual
   double
   stableDt() = 0;

   /**
    * @brief Advance the solution by a single time step.  Returns the
    *        simulation time advanced to.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    */
   virtual
   double
   advance(
      double a_time,
      double a_dt) = 0;

   /**
    * @brief Initialize the system given the current simulation time.
    *
    * @param[in] a_is_from_restart If true simulation is being run from
    *                              restart.
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_coll_diag_on If true, collision diagnostics are on.
    */
   virtual
   void
   initialize(
      bool a_is_from_restart,
      double a_time,
      double a_dt,
      bool a_coll_diag_on) = 0;

   /**
    * @brief Write the field and accumulated sequences to the field file(s).
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_probes Coordinates of probes.
    * @param[in] a_num_probes Number of probes.
    * @param[in] a_saved_seq Number of sequence time steps saved.
    * @param[in] a_saved_save Number of plots saved.
    * @param[in] a_time_hist_file_name Base name for time history files.
    */
   virtual
   void
   plot(
      double a_time,
      double a_dt,
      const vector<vector<double> >& a_probes, 
      int a_num_probes,
      int a_saved_seq,
      int& a_saved_save,
      string& a_time_hist_file_name) = 0;

   virtual
   void plotColl(
      double a_time,
      double a_dt,
      int a_saved_seq,
      int& a_saved_save,
      string& a_time_hist_file_name) = 0;

   /**
    * @brief Accumulate another time point for the sequences.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_probes Coordinates of probes.
    * @param[in] a_num_probes Number of probes.
    * @param[out] a_saved_seq Number of sequence time steps saved.
    */
   virtual
   void
   accumulateSequences(
      double a_time,
      const vector<vector<double> >& a_probes, 
      int a_num_probes,
      int& a_saved_seq) = 0;

   virtual
   void
   accumulateCollisionSequences(
      double a_time,
      double a_dt,
      int a_coll_op_idx,
      int& a_saved_seq) = 0;
  
   /**
    * @brief Returns the number of cells in the simulation.
    */
   virtual
   long int
   problemSize() const = 0;

   /**
    * @brief Updates ghost cells and boundary conditions in KineticSpecies
    *        distributions and fields.  Sets particle boundary conditions.
    *
    * @param[in] a_particles_only If true only updates particle boundary
    *                             conditions.
    */
   virtual
   void
   updateGhosts(
      bool a_particles_only = false) = 0;

   /**
    * @brief Print the parameters of the system being run .
    */
   virtual
   void
   printParameters() const = 0;

protected:
   // Unimplemented default constructor.
   System();

   // Parse input parameters from top level of input database.
   void
   parseParameters(
      LokiInputParser& a_pp)
   {
      string test_str = "true";
      a_pp.query("do_new_algorithm", test_str);
      m_do_new_algorithm = test_str.compare("true") == 0 ? true : false;

      test_str = "false";
      a_pp.query("plot_ke_fluxes", test_str);
      m_plot_ke_vel_bdy_flux = test_str.compare("true") == 0 ? true : false;

      // Check for constant Bz field. Added by IEO.
      if (a_pp.contains("bzconst")) {
         a_pp.get("bzconst", m_bz_const);
      }
      else {
         m_bz_const = 0.0;
      }
   }

   void
   PopulateInterspeciesData(
      KineticSpeciesPtrVect& ks);

   // Configuration space dimension.
   tbox::Dimension m_cdim;

   // Phase space dimension.
   tbox::Dimension m_pdim;

   // The global problem domain.
   tbox::Pointer<ProblemDomain> m_cfg_domain;

   // Names of all KineticSpecies.
   vector<string> m_species_names;

   // The charge density of each KineticSpecies.
   vector<ParallelArray> m_charge_density;

   // The electric field generated by all electric field drivers from each
   // KineticSpecies.
   vector<ParallelArray> m_ext_efield;

   // The total electric field generated by all electric field drivers from all
   // KineticSpecies.
   ParallelArray m_net_ext_efield;

   // The time sequences.
   vector<vector<double> > m_sequences, m_sequences_coll;

   // The times at which each time sequence has a data point.
   vector<double> m_time_seq, m_time_seq_coll;

   // The number of time sequences.
   int m_num_seq, m_num_seq_coll;

   // The number of data points for each time sequence.
   int m_length_seq, m_length_seq_coll;

   // If true use new, faster, algorithm.
   bool m_do_new_algorithm;

   // Order of spatial solution.
   int m_spatial_solution_order;

   // Order of temporal solution.
   int m_temporal_solution_order;

   // Background magnetic field
   double m_bz_const; //IEO

   // If true, plot species' velocity boundary Kinetic Energy fluxes.
   bool m_plot_ke_vel_bdy_flux;
};

}  // end namespace Loki

#endif
