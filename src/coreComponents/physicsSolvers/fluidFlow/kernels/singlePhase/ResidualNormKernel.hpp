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

/**
 * @file ResidualNormKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_RESIDUALNORMKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_RESIDUALNORMKERNEL_HPP

#include "physicsSolvers/SolverBaseKernels.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** ResidualNormKernel ********************************/

/**
 * @class IsothermalResidualNormKernel
 */
class IsothermalResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = solverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  IsothermalResidualNormKernel( globalIndex const rankOffset,
                                arrayView1d< real64 const > const & localResidual,
                                arrayView1d< globalIndex const > const & dofNumber,
                                arrayView1d< localIndex const > const & ghostRank,
                                ElementSubRegionBase const & subRegion,
                                real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_mass_n( subRegion.template getField< fields::flow::mass_n >() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );
    real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow] ) / massNormalizer;
    if( valMass > stack.localValue[0] )
    {
      stack.localValue[0] = valMass;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );
    stack.localValue[0] += m_localResidual[stack.localRow] * m_localResidual[stack.localRow];
    stack.localNormalizer[0] += massNormalizer;
  }


protected:

  /// View on mass at the previous converged time step
  arrayView1d< real64 const > const m_mass_n;

};

/**
 * @class ThermalResidualNormKernel
 */
class ThermalResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 2 >
{
public:

  using Base = solverBaseKernels::ResidualNormKernelBase< 2 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ThermalResidualNormKernel( globalIndex const rankOffset,
                             arrayView1d< real64 const > const & localResidual,
                             arrayView1d< globalIndex const > const & dofNumber,
                             arrayView1d< localIndex const > const & ghostRank,
                             ElementSubRegionBase const & subRegion,
                             real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_mass_n( subRegion.template getField< fields::flow::mass_n >() ),
    m_energy_n( subRegion.template getField< fields::flow::energy_n >() )
  {}

  GEOS_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const ei,
                                     real64 & massNormalizer,
                                     real64 & energyNormalizer ) const
  {
    massNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );
    energyNormalizer = LvArray::math::max( m_minNormalizer, LvArray::math::abs( m_energy_n[ei] ) ); // energy can be negative
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );

    // step 1: mass residual

    real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow] ) / massNormalizer;
    if( valMass > stack.localValue[0] )
    {
      stack.localValue[0] = valMass;
    }

    // step 2: energy residual
    real64 const valEnergy = LvArray::math::abs( m_localResidual[stack.localRow + 1] ) / energyNormalizer;
    if( valEnergy > stack.localValue[1] )
    {
      stack.localValue[1] = valEnergy;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );

    // step 1: mass residual

    stack.localValue[0] += m_localResidual[stack.localRow] * m_localResidual[stack.localRow];
    stack.localNormalizer[0] += massNormalizer;

    // step 2: energy residual

    stack.localValue[1] += m_localResidual[stack.localRow + 1] * m_localResidual[stack.localRow + 1];
    stack.localNormalizer[1] += energyNormalizer;
  }


protected:

  /// View on mass at the previous converged time step
  arrayView1d< real64 const > const m_mass_n;

  /// View on energy at the previous converged time step
  arrayView1d< real64 const > const m_energy_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch (isothermal version)
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( solverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    IsothermalResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, subRegion, minNormalizer );
    if( normType == solverBaseKernels::NormType::Linf )
    {
      IsothermalResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      IsothermalResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

  /**
   * @brief Create a new kernel and launch (thermal version)
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[in] solidInternalEnergy the solid internal energy model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( solverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[2],
                   real64 (& residualNormalizer)[2] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ThermalResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, subRegion, minNormalizer );
    if( normType == solverBaseKernels::NormType::Linf )
    {
      ThermalResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ThermalResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_RESIDUALNORMKERNEL_HPP
