/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ODE_INTEGRATOR_H_
#define _ODE_INTEGRATOR_H_

namespace Loki {

/**
 * Abstract base class that defines the interface to an ODE operator.  The
 * template parameter is the entity which must be time integrated.
 */
template <class TSoln>
class ODEInterface
{
public:
   /**
    * @brief Default constructor.
    */
   ODEInterface()
   {
   }

   /**
    * @brief Destructor.
    */
   virtual
   ~ODEInterface()
   {
   }

   /**
    * @brief Defines the TSoln data holder for data without ghost cells.
    *
    * @param[out] a_rhs Data holder for rhs (without ghost cells).
    * @param[in] a_prototype Prototype data holder.
    */
   virtual
   void
   defineRHSData(
      TSoln& a_rhs,
      const TSoln& a_prototype) = 0;

   /**
    * @brief Copies all data from srcSoln to dstSoln.
    *
    * @param[out] a_dstSoln Destination data holder.
    * @param[in] a_srcSoln Source data holder.
    */
   virtual
   void
   copySolnData(
      TSoln& a_dstSoln,
      const TSoln& a_srcSoln) = 0;

   // FIXME: soln should be const for an ODE but we are solving a PDE which
   // requires that its boundary conditions be set.  Instead of making a
   // copy of the soln inside evalRHS we are making soln non-const.
   /**
    * @brief Given the current time and solution, evaluates the RHS of the ODE.
    *
    * @param[out] a_rhs Data holder for rhs (without ghost cells).
    * @param[in] a_soln Current solution with which to evaluate RHS.
    * @param[in] a_time Current time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_first_rk_stage True if called from the first RK stage.
    * @param[in] a_last_rk_stage True if called from the last RK stage.
    */
   virtual
   void
   evalRHS(
      TSoln& a_rhs,
      TSoln& a_soln,
      double a_time,
      double a_dt,
      bool   a_first_rk_stage,
      bool   a_last_rk_stage) = 0;

   /**
    * @brief Checks for conformity between soln and protoSoln.
    *
    * @param[in] a_soln Data holder.
    * @param[in] a_protoSoln Prototype soln data holder.
    */
   virtual
   bool
   validSolnData(
      const TSoln& a_soln,
      const TSoln& a_protoSoln) = 0;

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
      const TSoln& a_rhs,
      const TSoln& a_protoSoln) = 0;

   /**
    * @brief Sets the TSoln data to zero.
    *
    * @param[out] a_soln Data holder.
    */
   virtual
   void
   zeroSolnData(
      TSoln& a_soln) = 0;

   /**
    * @brief Adds (scale * increment) to the TSoln data holder.
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
      TSoln&       a_soln,
      const TSoln& a_increment,
      double       a_scale,
      bool         a_sum_reduce_inc) = 0;

   /**
    * @brief Callback to allow for additional operations after stage advance.
    *
    * @param[in,out] a_soln Data holder.
    * @param[in] a_stage Which solution stage.
    */
   virtual
   void
   postStageAdvance(
      TSoln& a_soln,
      int a_stage) = 0;

private:
   // Unimplemented copy constructor.
   ODEInterface(
      const ODEInterface<TSoln>& other);

   // Unimplemented assignment operator.
   ODEInterface&
   operator = (
      const ODEInterface<TSoln>& rhs);

};

} // end namespace Loki

#endif
