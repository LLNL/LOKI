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
#ifndef  _VPSYSTEM_H_
#define  _VPSYSTEM_H_

#include "RKIntegrator.H"
#include "ODEInterface.H"
#include "System.H"
#include "VPState.H"

#include "Overture.h"
#include "Ogshow.h"
#include "ParmParse.H"
#include "ProblemDomain.H"
#include "Loki_Utilities.H"

namespace Loki {

/**
 * Implemention of the System and ODEInterface interfaces for the solution of
 * the coupled Vlasov-Poisson system.  The VPSystem time integrates a VPState
 * using the appropriate RKIntegrator.  The RKIntegrator makes use of the
 * ODEInterface methods implemented here as callbacks to perform the time
 * integration.
 */
class VPSystem : public System, public ODEInterface<VPState>, public Serializable
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp The input database.
    * @param[in] a_num_probes Number of probes for which the system will write
    *                         data.
    * @param[in] a_spatial_solution_order Spatial order of solution.
    * @param[in] a_temporal_solution_order Temporal order of solution.
    */
   VPSystem(
      ParmParse& a_pp,
      int a_num_probes,
      int a_spatial_solution_order,
      int a_temporal_solution_order);

   /**
    * @brief Destructor.
    */
   virtual
   ~VPSystem();

   //@{

   //! @name Implementation of ODEInterface<VPState> interface.

   /**
    * @brief Defines the VPState data holder for data without ghost cells.
    *
    * @param[out] a_rhs Data holder for rhs (without ghost cells).
    * @param[in] a_prototype Prototype data holder.
    */
   virtual
   void
   defineRHSData(
      VPState&       a_rhs,
      const VPState& a_prototype);

   /**
    * @brief Defines the VPState data holder for data with ghost cells.
    *
    * @param[out] a_soln Data holder for solution with ghost cells.
    * @param[in] a_prototype Prototype data holder.
    */
   virtual
   void
   defineSolnData(
      VPState&       a_soln,
      const VPState& a_prototype);

   /**
    * @brief Copies all data from srcSoln to dstSoln.
    *
    * @param[out] a_dstSoln Destination data holder.
    * @param[in] a_srcSoln Source data holder.
    */
   virtual
   void
   copySolnData(
      VPState&       a_dstSoln,
      const VPState& a_srcSoln);

   /**
    * @brief Sets the VPState data to zero.
    *
    * @param[out] a_soln Data holder.
    */
   virtual
   void
   zeroSolnData(
      VPState& a_soln);

   /**
    * @brief Adds (scale * increment) to the VPState data holder.
    *
    * @param[out] a_soln Solution data to which is added.
    * @param[in] a_increment Amount to add.
    * @param[in] a_scale Multiplier by which to scale increment.
    * @param[in] a_final_rk If true, this has been called as part of the final
    *                       RK update and tracking particle positions need to
    *                       be limited.
    */
   virtual
   void
   addSolnData(
      VPState&       a_soln,
      const VPState& a_increment,
      real           a_scale,
      bool           a_final_rk = false);

   /**
    * @brief Checks for conformity between soln and protoSoln ignoring ghost
    * cells.
    *
    * @param[in] a_rhs Data holder.
    * @param[in] a_protoSoln Prototype soln data holder.
    */
   virtual
   bool
   validRHSData(
      const VPState& a_rhs,
      const VPState& a_protoSoln);

   /**
    * @brief Checks for conformity between soln and protoSoln.
    *
    * @param[in] a_soln Data holder.
    * @param[in] a_protoSoln Prototype soln data holder.
    */
   virtual
   bool
   validSolnData(
      const VPState& a_soln,
      const VPState& a_protoSoln);

   /**
    * @brief Given the current time and solution, evaluates the RHS of the ODE.
    *
    * @param[out] a_rhs Data holder for rhs (without ghost cells).
    * @param[in] a_soln Current solution with which to evaluate RHS.
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_stage Current RK stage number.
    */
   virtual
   void
   evalRHS(
      VPState& a_rhs,
      VPState& a_soln,
      real     a_time,
      real     a_dt,
      int      a_stage);

   /**
    * @brief Callback to allow for additional operations after stage advance.
    *
    * @param[in,out] a_soln Data holder.
    * @param[in] a_stage Which solution stage.
    */
   virtual
   void
   postStageAdvance(
      VPState& a_soln,
      int a_stage);

   //@}

   //@{

   //! @name Implementation of System interface.

   /**
    * @brief Computes and returns a stable time step estimate.
    */
   virtual
   Real
   stableDt();

   /**
    * @brief Advance the solution by a single time step.  Returns the
    * simulation time advanced to.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    */
   virtual
   Real
   advance(
      real a_time,
      real a_dt);

   /**
    * @brief Initialize the system given the current simulation time.
    *
    * @param[in] a_is_from_restart If true simulation is being run from
    *                              restart.
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    */
   virtual
   void
   initialize(
      bool a_is_from_restart,
      real a_time,
      real a_dt);

   /**
    * @brief Returns number of fields to be plotted.
    */
   virtual
   int
   getNumFrameSeries();

   /**
    * @brief Write the field and accumulated sequences to the Overture show
    * file.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_probes Coordinates of probes.
    * @param[in] a_num_probes Number of probes.
    * @param[in] a_saved_seq Number of sequence time steps saved.
    * @param[in] a_saved_save Number of plots saved.
    * @param[out] a_show The Overture show file.
    */
   virtual
   void
   plot(
      real a_time,
      real a_dt,
      const RealArray& a_probes,
      int a_num_probes,
      int a_saved_seq,
      int& a_saved_save,
      Ogshow& a_show);

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
      real a_time,
      const RealArray& a_probes,
      int a_num_probes,
      int& a_saved_seq);

   /**
    * @brief Returns the number of cells in the simulation.
    */
   virtual
   long int
   problemSize() const;

   /**
    * @brief Updates ghost cells in KineticSpecies distributions.
    */
   virtual
   void
   updateGhosts();

   //@}

   //@{

   //! @name Implementation of Serializable interface.

   /**
    * @brief Read the KineticSpecies names from the restart database.
    *
    * @param[in] a_db The restart database.
    */
   virtual
   void
   getFromRestart(
      const HDF_DataBase& a_db);

   /**
    * @brief Write the KineticSpecies names to the restart database.
    *
    * @param[out] a_db The restart database.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      HDF_DataBase& a_db,
      real a_time);

   //@}

private:
   // Unimplemented default constructor.
   VPSystem();

   // Unimplemented copy constructor.
   VPSystem(
      const VPSystem& other);

   // Unimplemented assignment operator.
   VPSystem&
   operator = (
      const VPSystem& rhs);

   // Creates the state variable, m_state.
   void
   createVPState(
      ParmParse& a_pp);

   // Load balance the contents of m_state.
   void
   loadBalance();

   // Print class' parameters.
   void
   printParameters() const
   {
   }

   // Parse input parameters from top level of input database.
   void
   parseParameters(
      ParmParse& a_pp)
   {
      aString test_str = "true";
      a_pp.query("do_new_algorithm", test_str);
      m_do_new_algorithm = test_str.matches("true") ? true : false;
   }

   // Print the decomposition of m_state.
   void
   printDecomposition() const;

   // Configuration space dimension.
   tbox::Dimension m_cdim;

   // Phase space dimension.
   tbox::Dimension m_pdim;

   // The global problem domain.
   tbox::Pointer<ProblemDomain> m_cfg_domain;

   // Updated state.
   VPState m_state;

   // State at last time step.
   VPState m_state_old;

   // Names of all KineticSpecies.
   std::vector<string> m_species_names;

   // RKIntegrator to solve ODE system.
   tbox::Pointer<RKIntegrator<VPState> > m_integrator;

   // The net charge density.
   realArray m_net_charge_density;

   // The contribution to m_net_charge_density from each KineticSpecies.
   std::vector<realArray> m_charge_density;

   // The 2D species kinetic energy.
   realArray m_species_ke;

   // The 2D species physical boundary energy flux.
   realArray m_species_phys_bdry_flux;

   // The 2D species velocity boundary energy flux.
   realArray m_species_vel_bdry_flux;

   // The electric field solved by m_poisson.
   realArray m_efield;

   // The electric field generated by all electric field drivers from each
   // KineticSpecies.
   std::vector<realArray> m_ext_efield;

   // The total electric field generated by all electric field drivers from all
   // KineticSpecies.
   realArray m_net_ext_efield;

   // The time sequences.
   RealArray m_sequences;

   // The times at which each time sequence has a data point.
   RealArray m_time_seq;

   // The number of time sequences.
   int m_num_seq;

   // The number of data points for each time sequence.
   int m_length_seq;

   // If true use new, faster, algorithm.
   bool m_do_new_algorithm;

   // Order of spatial solution.
   int m_spatial_solution_order;

   // Order of temporal solution.
   int m_temporal_solution_order;
};

} // end namespace Loki

#endif
