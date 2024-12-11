/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_IMPLICITRATEANDSTATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_IMPLICITRATEANDSTATEKERNELS_HPP_

#include "RateAndStateKernelsBase.hpp"
#include "denseLinearAlgebra/denseLASolvers.hpp"

namespace geos
{

namespace rateAndStateKernels
{

/**
 * @class ImplicitFixedStressRateAndStateKernel
 *
 * @brief
 *
 * @details
 */
class ImplicitFixedStressRateAndStateKernel
{
public:

  ImplicitFixedStressRateAndStateKernel( SurfaceElementSubRegion & subRegion,
                                         constitutive::RateAndStateFriction const & frictionLaw,
                                         real64 const shearImpedance ):
    m_slipRate( subRegion.getField< fields::rateAndState::slipRate >() ),
    m_stateVariable( subRegion.getField< fields::rateAndState::stateVariable >() ),
    m_stateVariable_n( subRegion.getField< fields::rateAndState::stateVariable_n >() ),
    m_slipRate_n( subRegion.getField< fields::rateAndState::slipRate_n >() ),
    m_traction( subRegion.getField< fields::contact::traction >() ),
    m_slipVelocity( subRegion.getField< fields::rateAndState::slipVelocity >() ),
    m_shearImpedance( shearImpedance ),
    m_frictionLaw( frictionLaw.createKernelUpdates()  )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( )
    {}

    real64 jacobian[2][2]{};

    real64 rhs[2]{};

  };

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              real64 const dt,
              StackVariables & stack ) const
  {
    real64 const normalTraction = m_traction[k][0];
    real64 const shearTractionMagnitude = LvArray::math::sqrt( m_traction[k][1] * m_traction[k][1] + m_traction[k][2] * m_traction[k][2] );

    // Eq 1: Scalar force balance for slipRate and shear traction magnitude
    stack.rhs[0] = shearTractionMagnitude - m_shearImpedance * m_slipRate[k]
                   - normalTraction * m_frictionLaw.frictionCoefficient( k, m_slipRate[k], m_stateVariable[k] );
    real64 const dFriction[2] = { -normalTraction * m_frictionLaw.dFrictionCoefficient_dStateVariable( k, m_slipRate[k], m_stateVariable[k] ),
                                  -m_shearImpedance - normalTraction * m_frictionLaw.dFrictionCoefficient_dSlipRate( k, m_slipRate[k], m_stateVariable[k] ) };

    // Eq 2: slip law
    stack.rhs[1] = (m_stateVariable[k] - m_stateVariable_n[k]) / dt - m_frictionLaw.stateEvolution( k, m_slipRate[k], m_stateVariable[k] );
    real64 const dStateEvolutionLaw[2] = { 1.0 / dt - m_frictionLaw.dStateEvolution_dStateVariable( k, m_slipRate[k], m_stateVariable[k] ),
                                           -m_frictionLaw.dStateEvolution_dSlipRate( k, m_slipRate[k], m_stateVariable[k] ) };

    // Assemble Jacobian matrix
    stack.jacobian[0][0] = dFriction[0];          // derivative of Eq 1 w.r.t. stateVariable
    stack.jacobian[0][1] = dFriction[1];          // derivative of Eq 1 w.r.t. slipRate
    stack.jacobian[1][0] = dStateEvolutionLaw[0]; // derivative of Eq 2 w.r.t. stateVariable
    stack.jacobian[1][1] = dStateEvolutionLaw[1]; // derivative of Eq 2 w.r.t. slipRate
  }

  GEOS_HOST_DEVICE
  void solve( localIndex const k,
              StackVariables & stack ) const
  {
    /// Solve 2x2 system
    real64 solution[2] = {0.0, 0.0};
    LvArray::tensorOps::scale< 2 >( stack.rhs, -1.0 );
    denseLinearAlgebra::solve< 2 >( stack.jacobian, stack.rhs, solution );

    // Slip rate is bracketed between [0, shear traction magnitude / shear impedance]
    // Check that the update did not end outside of the bracket.
    real64 const shearTractionMagnitude = LvArray::math::sqrt( m_traction[k][1] * m_traction[k][1] + m_traction[k][2] * m_traction[k][2] );
    real64 const upperBound = shearTractionMagnitude / m_shearImpedance - m_slipRate[k];
    real64 const lowerBound = -m_slipRate[k];

    real64 scalingFactor = 1.0;
    if( solution[1] > upperBound )
    {
      scalingFactor = 0.5 * upperBound / solution[1];
    }
    else if( solution[1] < lowerBound )
    {
      scalingFactor = 0.5 * lowerBound / solution[1];
    }

    LvArray::tensorOps::scale< 2 >( solution, scalingFactor );

    // Update variables
    m_stateVariable[k]  +=  solution[0];
    m_slipRate[k]       +=  solution[1];
  }

  GEOS_HOST_DEVICE
  void projectSlipRate( localIndex const k ) const
  {
    real64 const frictionCoefficient = m_frictionLaw.frictionCoefficient( k, m_slipRate[k], m_stateVariable[k] );
    projectSlipRateBase( k, frictionCoefficient, m_shearImpedance, m_traction, m_slipRate, m_slipVelocity );
  }

  GEOS_HOST_DEVICE
  void udpateVariables( localIndex const k ) const
  {
    projectSlipRate( k );
    m_stateVariable_n[k] = m_stateVariable[k];
    m_slipRate_n[k] = m_slipRate[k];
  }

  GEOS_HOST_DEVICE
  camp::tuple< int, real64 > checkConvergence( StackVariables const & stack,
                                               real64 const tol ) const
  {
    real64 const residualNorm = LvArray::tensorOps::l2Norm< 2 >( stack.rhs );
    int const converged = residualNorm < tol ? 1 : 0;
    camp::tuple< int, real64 > result { converged, residualNorm };
    return result;
  }

  GEOS_HOST_DEVICE
  void resetState( localIndex const k ) const
  {
    m_stateVariable[k] = m_stateVariable_n[k];
    m_slipRate[k] = m_slipRate_n[k];
  }

template< typename POLICY >
static real64 solveRateAndStateEquation( SurfaceElementSubRegion & subRegion,
                                         ImplicitFixedStressRateAndStateKernel & kernel,
                                         real64 dt,
                                         integer const maxNewtonIter,
                                         real64 const newtonTol )
{
  bool converged = false;
  for( integer attempt = 0; attempt < 5; attempt++ )
  {
    if( attempt > 0 )
    {
      forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        kernel.resetState( k );
      } );
    }
    GEOS_LOG_RANK_0( GEOS_FMT( "  Attempt {} ", attempt ) );
    converged = newtonSolve< POLICY >( subRegion, kernel, dt, maxNewtonIter, newtonTol );
    if( converged )
    {
      forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        kernel.udpateVariables( k );
      } );
      return dt;
    }
    else
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "  Attempt {} failed. Halving dt and retrying.", attempt ) );
      dt *= 0.5;
    }
  }
  if( !converged )
  {
    GEOS_ERROR( "Maximum number of attempts reached without convergence." );
  }
  return dt;
}

private:

  arrayView1d< real64 > const m_slipRate;

  arrayView1d< real64 > const m_stateVariable;

  arrayView1d< real64 > const m_stateVariable_n;

  arrayView1d< real64 > const m_slipRate_n;

  arrayView2d< real64 const > const m_traction;

  arrayView2d< real64 > const m_slipVelocity;

  real64 const m_shearImpedance;

  constitutive::RateAndStateFriction::KernelWrapper m_frictionLaw;

};

} /* namespace rateAndStateKernels */

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_IMPLICITRATEANDSTATEKERNELS_HPP_ */
