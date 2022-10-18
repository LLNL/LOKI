/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef  _VMSYSTEM_H_
#define  _VMSYSTEM_H_

#include "RKIntegrator.H"
#include "ODEInterface.H"
#include "RestartReader.H"
#include "RestartWriter.H"
#include "System.H"
#include "VMState.H"
#include "Loki_Utilities.H"

namespace Loki {

/**
 * Implemention of the System and ODEInterface interfaces for the solution of
 * the coupled Vlasov-Maxwell system.  The VMSystem time integrates a VMState
 * using the appropriate RKIntegrator.  The RKIntegrator makes use of the
 * ODEInterface methods implemented here as callbacks to perform the time
 * integration.
 */
class VMSystem : public System, public ODEInterface<VMState>, public Serializable
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
   VMSystem(
      LokiInputParser& a_pp,
      int a_num_probes,
      int a_spatial_solution_order,
      int a_temporal_solution_order,
      bool a_coll_diag_on);

   /**
    * @brief Destructor.
    */
   virtual
   ~VMSystem();

   //@{

   //! @name Implementation of ODEInterface<VMState> interface.

   /**
    * @brief Defines the VMState data holder for data without ghost cells.
    *
    * @param[out] a_rhs Data holder for rhs (without ghost cells).
    * @param[in] a_prototype Prototype data holder.
    */
   virtual
   void
   defineRHSData(
      VMState&       a_rhs,
      const VMState& a_prototype);

   /**
    * @brief Copies all data from srcSoln to dstSoln.
    *
    * @param[out] a_dstSoln Destination data holder.
    * @param[in] a_srcSoln Source data holder.
    */
   virtual
   void
   copySolnData(
      VMState&       a_dstSoln,
      const VMState& a_srcSoln);

   /**
    * @brief Sets the VMState data to zero.
    *
    * @param[out] a_soln Data holder.
    */
   virtual
   void
   zeroSolnData(
      VMState& a_soln);

   /**
    * @brief Adds (scale * increment) to the VMState data holder.
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
      VMState&       a_soln,
      const VMState& a_increment,
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
      const VMState& a_rhs,
      const VMState& a_protoSoln);

   /**
    * @brief Checks for conformity between soln and protoSoln.
    *
    * @param[in] a_soln Data holder.
    * @param[in] a_protoSoln Prototype soln data holder.
    */
   virtual
   bool
   validSolnData(
      const VMState& a_soln,
      const VMState& a_protoSoln);

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
      VMState& a_rhs,
      VMState& a_soln,
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
      VMState& a_soln,
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
      int& a_save_save,
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
   VMSystem();

   // Unimplemented copy constructor.
   VMSystem(
      const VMSystem& other);

   // Unimplemented assignment operator.
   VMSystem&
   operator = (
      const VMSystem& rhs);

   // Creates the state variable, m_state.
   void
   createVMState(
      LokiInputParser& a_pp,
      bool a_coll_diag_on);

   // Load balance the contents of m_state.
   void
   loadBalance();

   // Print the decomposition of m_state.
   void
   printDecomposition() const;

   // Updated state.
   VMState m_state;

   // State at last time step.
   VMState m_state_old;

   // RKIntegrator to solve ODE system.
   tbox::Pointer<RKIntegrator<VMState> > m_integrator;

   // The net X, Y, and Z component of the current density.
   ParallelArray m_net_x_current_density;
   ParallelArray m_net_y_current_density;
   ParallelArray m_net_z_current_density;

   // The contribution to m_net_(xyz)_current_density from each KineticSpecies.
   vector<ParallelArray> m_x_current_density;
   vector<ParallelArray> m_y_current_density;
   vector<ParallelArray> m_z_current_density;

   // The square of the effective plasma frequency.
   // This is sum(Z^2*N/m)=sum(Z*charge_density/m) over all species.
   ParallelArray m_omega_eff2;
};

} // end namespace Loki

#endif
