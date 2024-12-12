/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_QUASIDYNAMICEQRK32_HPP
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_QUASIDYNAMICEQRK32_HPP

#include "QDRateAndStateBase.hpp"
#include "kernels/ExplicitRateAndStateKernels.hpp"

namespace geos
{

class ExplicitQDRateAndState : public QDRateAndStateBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  ExplicitQDRateAndState() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  ExplicitQDRateAndState( const string & name,
                          Group * const parent );

  /// Destructor
  virtual ~ExplicitQDRateAndState() override;

  static string derivedSolverPrefix() { return "Explicit";};

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  struct viewKeyStruct : public QDRateAndStateBase::viewKeyStruct
  {
    /// target slip increment
    constexpr static char const * timeStepTol() { return "timeStepTol"; }
  };

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override final;


  virtual real64 setNextDt( real64 const & currentDt,
                            DomainPartition & domain ) override final;

  /**
   * @brief Evaluates whether an adaptive time step was successful
   * @param domain
   */
  void evalTimestep(DomainPartition & domain );

  /**
   * @brief Computes stage rates for the initial Runge-Kutta substage and updates slip and state
   * @param dt
   * @param domain
   */
  void stepRateStateODEInitialSubstage( real64 const dt, DomainPartition & domain ) const;

  /**
   * @brief Computes stage rates at the Runge-Kutta substage specified by stageIndex and updates slip and state
   * @param stageIndex
   * @param dt
   * @param domain
   */
  void stepRateStateODESubstage( integer const stageIndex,
                                 real64 const dt,
                                 DomainPartition & domain ) const;

  /**
   * @brief Updates slip and state to t + dt and approximates the error
   * @param dt
   * @param domain
   */
  void stepRateStateODEAndComputeError( real64 const dt, DomainPartition & domain ) const;

  /**
   * @brief Updates rate-and-state slip velocity
   * @param domain
   */
  void updateSlipVelocity( real64 const & time_n,
                           real64 const & dt,
                           DomainPartition & domain ) const;

protected:

  /// Runge-Kutta Butcher table (specifies the embedded RK method)
  // TODO: The specific type should not be hardcoded!
  // Should be possible to change RK-method based on the table.
  rateAndStateKernels::BogackiShampine32Table m_butcherTable;

  bool m_successfulStep; // Flag indicating if the adative time step was accepted

  real64 m_stepUpdateFactor; // Factor to update timestep with

  /**
   * @brief Proportional-integral-derivative controller used for updating time step
   * based error estimate in the current and previous time steps.
   */
  class PIDController
  {
public:

    GEOS_HOST_DEVICE
    PIDController( std::array< const real64, 3 > const & cparams,
                   const real64 atol,
                   const real64 rtol,
                   const real64 safety ):
      controlParameters{ cparams },
      absTol( atol ),
      relTol( rtol ),
      acceptSafety( safety ),
      errors{ {0.0, 0.0, 0.0} }
    {}

    /// Default copy constructor
    PIDController( PIDController const & ) = default;

    /// Default move constructor
    PIDController( PIDController && ) = default;

    /// Deleted default constructor
    PIDController() = delete;

    /// Deleted copy assignment operator
    PIDController & operator=( PIDController const & ) = delete;

    /// Deleted move assignment operator
    PIDController & operator=( PIDController && ) =  delete;

    /// Parameters for the PID error controller
    const std::array< const real64, 3 > controlParameters; // Controller parameters

    real64 const absTol; // absolut tolerence

    real64 const relTol; // relative tolerence

    real64 const acceptSafety; // Acceptance safety

    std::array< real64, 3 > errors; // Errors for current and two previous updates
                                    // stored as [n+1, n, n-1]

    real64 computeUpdateFactor( integer const algHighOrder, integer const algLowOrder )
    {
      // PID error controller + limiter
      real64 const k = LvArray::math::min( algHighOrder, algLowOrder ) + 1.0;
      real64 const eps0 = 1.0/(errors[0] + std::numeric_limits< real64 >::epsilon()); // n + 1
      real64 const eps1 = 1.0/(errors[1] + std::numeric_limits< real64 >::epsilon()); // n
      real64 const eps2 = 1.0/(errors[2] + std::numeric_limits< real64 >::epsilon()); // n-1
      // Compute update factor eps0^(beta0/k)*eps1^(beta1/k)*eps2^(beta2/k) where
      // beta0 - beta2 are the control parameters. Also apply limiter to smoothen changes.
      // Limiter is 1.0 + atan(x - 1.0). Here use atan(x) = atan2(x, 1.0).
      return 1.0 + LvArray::math::atan2( pow( eps0, controlParameters[0] / k ) *
                                         pow( eps1, controlParameters[1] / k ) *
                                         pow( eps2, controlParameters[2] / k ) - 1.0, 1.0 );
    }
  };

  PIDController m_controller;

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_QUASIDYNAMICEQRK32_HPP */
