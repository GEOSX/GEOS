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

#ifndef GEOS_PHYSICSSOLVERS_RATEANDSTATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_RATEANDSTATEKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "physicsSolvers/inducedSeismicity/rateAndStateFields.hpp"
#include "denseLinearAlgebra/denseLASolvers.hpp"

namespace geos
{

namespace rateAndStateKernels
{

// TBD: Pass the kernel and add getters for relevant fields to make this function general purpose and avoid
// wrappers?
GEOS_HOST_DEVICE
static void projectSlipRateBase( localIndex const k,
                          real64 const frictionCoefficient,
                          real64 const shearImpedance,
                          arrayView2d< real64 const > const traction,
                          arrayView1d< real64 const > const slipRate,
                          arrayView2d< real64 > const slipVelocity)
{
  // Project slip rate onto shear traction to get slip velocity components
  real64 const frictionForce = traction[k][0] * frictionCoefficient;
  real64 const projectionScaling = 1.0 / ( shearImpedance +  frictionForce / slipRate[k] );
  slipVelocity[k][0] = projectionScaling * traction[k][1];
  slipVelocity[k][1] = projectionScaling * traction[k][2];
}

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
    real64 const dStateEvolutionLaw[2] = { 1 / dt - m_frictionLaw.dStateEvolution_dStateVariable( k, m_slipRate[k], m_stateVariable[k] ),
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
    denseLinearAlgebra::solve< 2 >( stack.jacobian, stack.rhs, solution );

    // Update variables
    m_stateVariable[k]  -=  solution[0];
    m_slipRate[k]       -=  solution[1];
  }

  GEOS_HOST_DEVICE
  void projectSlipRate( localIndex const k ) const
  {
    real64 const frictionCoefficient = m_frictionLaw.frictionCoefficient( k, m_slipRate[k], m_stateVariable[k] );
    projectSlipRateBase(k, frictionCoefficient, m_shearImpedance, m_traction, m_slipRate, m_slipVelocity);
  }

  GEOS_HOST_DEVICE
  std::pair< int, real64 > checkConvergence( StackVariables const & stack,
                                             real64 const tol ) const
  {
    real64 const residualNorm = LvArray::tensorOps::l2Norm< 2 >( stack.rhs );
    int const converged = residualNorm < tol ? 1 : 0;
    return std::make_pair( converged, residualNorm );
  }

private:

  arrayView1d< real64 > const m_slipRate;

  arrayView1d< real64 > const m_stateVariable;

  arrayView1d< real64 const > const m_stateVariable_n;

  arrayView2d< real64 const > const m_traction;

  arrayView2d< real64 > const m_slipVelocity;

  real64 const m_shearImpedance;

  constitutive::RateAndStateFriction::KernelWrapper m_frictionLaw;

};

/**
 * @class ExplicitRateAndStateKernel
 *
 * @brief
 *
 * @details
 */
class ExplicitRateAndStateKernel
{
public:

  ExplicitRateAndStateKernel( SurfaceElementSubRegion & subRegion,
                              constitutive::RateAndStateFriction const & frictionLaw,
                              real64 const shearImpedance ):
    m_slipRate( subRegion.getField< fields::rateAndState::slipRate >() ),
    m_stateVariable( subRegion.getField< fields::rateAndState::stateVariable >() ),
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

    real64 jacobian;
    real64 rhs;

  };

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              real64 const dt,
              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( dt );
    real64 const normalTraction = m_traction[k][0];
    real64 const shearTractionMagnitude = LvArray::math::sqrt( m_traction[k][1] * m_traction[k][1] + m_traction[k][2] * m_traction[k][2] );

    // Slip rate is bracketed between [0, shear traction magnitude / shear impedance] 
    // If slip rate is outside the bracket, re-initialize to the middle value
    real64 const upperBound = shearTractionMagnitude/m_shearImpedance;
    real64 const bracketedSlipRate = m_slipRate[k] > upperBound ? upperBound/2 : m_slipRate[k];
  
    stack.rhs = shearTractionMagnitude - m_shearImpedance *bracketedSlipRate - normalTraction * m_frictionLaw.frictionCoefficient( k,bracketedSlipRate, m_stateVariable[k] );
    stack.jacobian = -m_shearImpedance - normalTraction * m_frictionLaw.dFrictionCoefficient_dSlipRate( k, bracketedSlipRate, m_stateVariable[k] );
  }

  GEOS_HOST_DEVICE
  void solve( localIndex const k,
              StackVariables & stack) const
  {   
    m_slipRate[k] -= stack.rhs/stack.jacobian;
  }

  
  GEOS_HOST_DEVICE
  std::pair< int, real64 > checkConvergence( StackVariables const & stack,
                                             real64 const tol ) const
  {
    real64 const residualNorm = LvArray::math::abs(stack.rhs);
    int const converged = residualNorm < tol ? 1 : 0;
    return std::make_pair( converged, residualNorm );
  }

  GEOS_HOST_DEVICE
  void projectSlipRate( localIndex const k ) const
  {
    real64 const frictionCoefficient = m_frictionLaw.frictionCoefficient( k, m_slipRate[k], m_stateVariable[k] );
    projectSlipRateBase(k, frictionCoefficient, m_shearImpedance, m_traction, m_slipRate, m_slipVelocity);
  }

private:

  arrayView1d< real64 > const m_slipRate;

  arrayView1d< real64 > const m_stateVariable;

  arrayView2d< real64 const > const m_traction;

  arrayView2d< real64 > const m_slipVelocity;

  real64 const m_shearImpedance;

  constitutive::RateAndStateFriction::KernelWrapper m_frictionLaw;

};


/**
 * @brief Performs the kernel launch
 * @tparam POLICY the policy used in the RAJA kernels
 */
template< typename KernelType, typename POLICY >
static void
createAndLaunch( SurfaceElementSubRegion & subRegion,
                 string const & frictionLawNameKey,
                 real64 const shearImpedance,
                 integer const maxIterNewton,
                 real64 const newtonTol,
                 real64 const time_n,
                 real64 const dt )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n );

  string const & frictionaLawName = subRegion.getReference< string >( frictionLawNameKey );
  constitutive::RateAndStateFriction const & frictionLaw = subRegion.getConstitutiveModel< constitutive::RateAndStateFriction >( frictionaLawName );
  KernelType kernel( subRegion, frictionLaw, shearImpedance );

  // Newton loop (outside of the kernel launch)
  bool allConverged = false;
  for( integer iter = 0; iter < maxIterNewton; iter++ )
  {
    RAJA::ReduceMin< parallelDeviceReduce, int > converged( 1 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > residualNorm( 0.0 );
    forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KernelType::StackVariables stack;
      kernel.setup( k, dt, stack );
      kernel.solve( k, stack );
      auto result = kernel.checkConvergence( stack, newtonTol );
      converged.min( std::get< 0 >( result ) );
      residualNorm.max( std::get< 1 >( result ) );
    } );

    real64 const maxResidualNorm = MpiWrapper::max( residualNorm.get() );
    GEOS_LOG_RANK_0( GEOS_FMT( "-----iter {} : residual = {:.10e} ", iter, maxResidualNorm ) );

    if( converged.get() )
    {
      allConverged = true;
      break;
    }
  }
  if( !allConverged )
  {
    GEOS_ERROR( " Failed to converge" );
  }
  forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    kernel.projectSlipRate( k );
  } );
}

} /* namespace rateAndStateKernels */

}/* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_RATEANDSTATEKERNELS_HPP_ */
