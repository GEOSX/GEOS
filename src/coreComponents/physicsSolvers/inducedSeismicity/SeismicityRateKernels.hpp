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

#ifndef GEOS_PHYSICSSOLVERS_SEISMICITYRATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_SEISMICITYRATEKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "physicsSolvers/inducedSeismicity/inducedSeismicityFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geos
{

namespace seismicityRateKernels
{
/**
 * @class DieterichSeismicityRate
 *
 * @brief Solving the ODE for seismicity rate from Dieterich, 1994
 *
 * @details This solver finds a solution R(x, t) - the seismicity rate - to the ordinary differential equation (ODE)
 * formulated by Dieterich, 1994 given a certain stressing history. The stressing history can consist
 * of mechanical stresses and pore pressure. The solver class includes a member variable
 * pointing to the stress solver that is specified in the XML file. SolverStep for the
 * stress solver is then called in the SolverStep function for the seismicity rate, to take
 * the updated stress history as the input.
 *
 * Solving the ODE is currently implemented by computing the closed-form integral solution
 * to the ODE which involves numerical calculation of an integral of a stress functional.
 * We initially solve for the log of the seismicity rate in order to avoid overflow that
 * typically occurs in the exponential of the stress history.
 */
class SeismicityRateKernel
{
public:

  SeismicityRateKernel( ElementSubRegionBase & subRegion ):
    m_R( subRegion.getField< fields::inducedSeismicity::seismicityRate >() ),
    m_logDenom( subRegion.getField< fields::inducedSeismicity::logDenom >() ),
    m_sigma_0( subRegion.getField< fields::inducedSeismicity::initialProjectedNormalTraction >() ),
    m_sigma_n( subRegion.getField< fields::inducedSeismicity::projectedNormalTraction_n >() ),
    m_sigma( subRegion.getField< fields::inducedSeismicity::projectedNormalTraction >() ),
    m_tau_0( subRegion.getField< fields::inducedSeismicity::initialProjectedShearTraction >() ),
    m_tau_n( subRegion.getField< fields::inducedSeismicity::projectedShearTraction_n >() ),
    m_tau( subRegion.getField< fields::inducedSeismicity::projectedShearTraction >() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( real64 const directEffect,
                    real64 const backgroundStressingRate ):
      directEffectValue( directEffect ),
      backgroundStressingRateValue( backgroundStressingRate ),
      effectiveNormalTraction_0( 0.0 ),
      effectiveNormalTraction_n( 0.0 ),
      effectiveNormalTraction( 0.0 )
    {}

    real64 const directEffectValue;

    real64 const backgroundStressingRateValue;

    real64 effectiveNormalTraction_0;

    real64 effectiveNormalTraction_n;

    real64 effectiveNormalTraction;
  };

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    stack.effectiveNormalTraction_0 = -m_sigma_0[k];
    stack.effectiveNormalTraction_n = -m_sigma_n[k];
    stack.effectiveNormalTraction   = -m_sigma[k];
  }

  GEOS_HOST_DEVICE
  void computeSeismicityRate( localIndex const k,
                              real64 const & time_n,
                              real64 const & dt,
                              StackVariables & stack ) const
  {

    // arguments of stress exponential at current and previous time step
    real64 const g = ( m_tau[k] + stack.backgroundStressingRateValue*(time_n+dt))/(stack.directEffectValue*stack.effectiveNormalTraction )
                     - m_tau_0[k]/(stack.directEffectValue * stack.effectiveNormalTraction_0 );

    real64 const g_n = ( m_tau_n[k] + stack.backgroundStressingRateValue*time_n)/(stack.directEffectValue*stack.effectiveNormalTraction_n )
                       - m_tau_0[k]/(stack.directEffectValue*stack.effectiveNormalTraction_0);

    // Compute the difference of the log of the denominator of closed for integral solution.
    // This avoids directly computing the exponential of the current stress state which is more prone to overflow.
    m_logDenom[k] += std::log( 1 + dt/(2*(stack.directEffectValue*stack.effectiveNormalTraction_0/stack.backgroundStressingRateValue))
                               *(std::exp( g - m_logDenom[k] ) + std::exp( g_n - m_logDenom[k] ) ));

    // Convert log seismicity rate to raw value
    m_R[k] = LvArray::math::exp( g - m_logDenom[k] );
  }

protected:

  arrayView1d< real64 > m_R;

  arrayView1d< real64 > m_logDenom;

  arrayView1d< real64 const > m_sigma_0;

  arrayView1d< real64 const > m_sigma_n;

  arrayView1d< real64 const > m_sigma;

  arrayView1d< real64 const > m_tau_0;

  arrayView1d< real64 const > m_tau_n;

  arrayView1d< real64 const > m_tau;
};


class SeismicityRateKernelPoroelastic : public SeismicityRateKernel
{

public:

  SeismicityRateKernelPoroelastic( ElementSubRegionBase & subRegion ):
    SeismicityRateKernel( subRegion ),
    m_pressure_0( subRegion.getField< fields::flow::initialPressure >() ),
    m_pressure_n( subRegion.getField< fields::flow::pressure_n >() ),
    m_pressure( subRegion.getField< fields::flow::pressure >() )
  {}

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    stack.effectiveNormalTraction_0 = -m_sigma_0[k] - m_pressure_0[k];
    stack.effectiveNormalTraction_n = -m_sigma_n[k] - m_pressure_n[k];
    stack.effectiveNormalTraction   = -m_sigma[k] - m_pressure[k];
  }

private:

  arrayView1d< real64 const > m_pressure_0;

  arrayView1d< real64 const > m_pressure_n;

  arrayView1d< real64 const > m_pressure;
};


/**
 * @brief Performs the kernel launch
 * @tparam POLICY the policy used in the RAJA kernels
 */
template< typename POLICY, bool ISPORO >
static void
createAndLaunch( ElementSubRegionBase & subRegion,
                 real64 const time_n,
                 real64 const dt,
                 real64 const directEffectValue,
                 real64 const backgroundStressingRateValue )
{
  GEOS_MARK_FUNCTION;

  using kernelType = std::conditional_t< ISPORO, SeismicityRateKernelPoroelastic, SeismicityRateKernel >;
  kernelType kernel( subRegion );

  forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    typename kernelType::StackVariables stack( directEffectValue, backgroundStressingRateValue );
    kernel.setup( k, stack );
    kernel.computeSeismicityRate( k, time_n, dt, stack );
  } );
}

} /* namespace seismicityRateKernels */

}/* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SEISMICITYRATEKERNELS_HPP_ */
