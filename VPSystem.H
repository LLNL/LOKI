/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef  _VPSYSTEM_H_
#define  _VPSYSTEM_H_

#include "RKIntegrator.H"
#include "ODEInterface.H"
#include "RestartReader.H"
#include "RestartWriter.H"
#include "System.H"
#include "VPState.H"
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
    * @param[in] a_coll_diag_on If true then collision diagnostics are
    *                           requested.
    */
   VPSystem(
      LokiInputParser& a_pp,
      int a_num_probes,
      int a_spatial_solution_order,
      int a_temporal_solution_order,
      bool a_coll_diag_on);

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
    * @param[in] a_sum_reduce_inc If true, sum reduce the particle increment.
    *                             This is necessary at different points in the
    *                             RK schemes.
    */
   virtual
   void
   addSolnData(
      VPState&       a_soln,
      const VPState& a_increment,
      double         a_scale,
      bool           a_sum_reduce_inc);

   /**
    * @brief Checks for conformity between soln and protoSoln ignoring ghost
    *        cells.
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
    * @param[in] a_first_rk_stage True if called from the first RK stage.
    * @param[in] a_last_rk_stage True if called from the last RK stage.
    */
   virtual
   void
   evalRHS(
      VPState& a_rhs,
      VPState& a_soln,
      double   a_time,
      double   a_dt,
      bool     a_first_rk_stage,
      bool     a_last_rk_stage);

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
   double
   stableDt();

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
      double a_dt);

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
      bool a_coll_diag_on);

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
      string& a_time_hist_file_name);

   virtual
   void
   plotColl(
      double a_time,
      double a_dt,
      int a_saved_seq,
      int& a_saved_save,
      string& a_time_hist_file_name);

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
      int& a_saved_seq);

   virtual
   void
   accumulateCollisionSequences(
      double a_time,
      double a_dt,
      int a_coll_op_idx,
      int& a_saved_seq);

   /**
    * @brief Returns the number of cells in the simulation.
    */
   virtual
   long int
   problemSize() const;

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
      bool a_particles_only = false);

   /**
    * @brief Print the parameters of the system being run .
    */
   virtual
   void
   printParameters() const;

   //@}

   //@{

   //! @name Implementation of Serializable interface.

   /**
    * @brief Read the KineticSpecies names from the restart database.
    *
    * @param[in] a_reader The restart file reader.
    */
   virtual
   void
   getFromRestart(
      RestartReader& a_reader);

   /**
    * @brief Write the KineticSpecies names to the restart database.
    *
    * @param[out] a_writer The restart file writer.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      RestartWriter& a_writer,
      double a_time);

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
      LokiInputParser& a_pp,
      bool a_coll_diag_on);

   // Load balance the contents of m_state.
   void
   loadBalance();

   // Print the decomposition of m_state.
   void
   printDecomposition() const;

   // Updated state.
   VPState m_state;

   // State at last time step.
   VPState m_state_old;

   // RKIntegrator to solve ODE system.
   tbox::Pointer<RKIntegrator<VPState> > m_integrator;

   // The net charge density.
   ParallelArray m_net_charge_density;

   // The 2D species x momentum
   ParallelArray m_species_momx;

   // The 2D species y momentum
   ParallelArray m_species_momy;

   // The 2D species kinetic energy.
   ParallelArray m_species_ke;

   // The 2D species entropy
   ParallelArray m_species_ent;
};

} // end namespace Loki

#endif
