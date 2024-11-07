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
 * @file ThermalResidualNormKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALRESIDUALNORMKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALRESIDUALNORMKERNEL_HPP

#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 3 >
{
public:

  using Base = ResidualNormKernelBase< 3 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numComponents,
                      integer const numPhases,
                      ElementSubRegionBase const & subRegion,
                      constitutive::MultiFluidBase const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      constitutive::SolidInternalEnergy const & solidInternalEnergy,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numComponents( numComponents ),
    m_numPhases( numPhases ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_phaseVolFrac_n( subRegion.getField< fields::flow::phaseVolumeFraction_n >() ),
    m_totalDens_n( fluid.totalDensity_n() ),
    m_phaseDens_n( fluid.phaseDensity_n() ),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n() ),
    m_solidInternalEnergy_n( solidInternalEnergy.getInternalEnergy_n() )
  {}

  GEOS_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const ei,
                                     real64 & massNormalizer,
                                     real64 & energyNormalizer ) const
  {
    massNormalizer = LvArray::math::max( m_minNormalizer, m_totalDens_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );
    real64 const poreVolume = m_porosity_n[ei][0] * m_volume[ei];
    energyNormalizer = m_solidInternalEnergy_n[ei][0] * ( 1.0 - m_porosity_n[ei][0] ) * m_volume[ei];
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      energyNormalizer += m_phaseInternalEnergy_n[ei][0][ip] * m_phaseDens_n[ei][0][ip] * m_phaseVolFrac_n[ei][ip] * poreVolume;
    }
    // warning: internal energy can be negative
    energyNormalizer = LvArray::math::max( m_minNormalizer, LvArray::math::abs( energyNormalizer ) );
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );
    real64 const volumeNormalizer = LvArray::math::max( m_minNormalizer, m_porosity_n[ei][0] * m_volume[ei] );

    // step 1: mass residual

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
      if( valMass > stack.localValue[0] )
      {
        stack.localValue[0] = valMass;
      }
    }

    // step 2: volume residual

    real64 const valVol = LvArray::math::abs( m_localResidual[stack.localRow + m_numComponents] ) / volumeNormalizer;
    if( valVol > stack.localValue[1] )
    {
      stack.localValue[1] = valVol;
    }

    // step 3: energy residual

    real64 const valEnergy = LvArray::math::abs( m_localResidual[stack.localRow + m_numComponents + 1] ) / energyNormalizer;
    if( valEnergy > stack.localValue[2] )
    {
      stack.localValue[2] = valEnergy;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    // note: for the L2 norm, we bundle the volume and mass residuals/normalizers
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );

    // step 1: mass residual

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      stack.localValue[0] += m_localResidual[stack.localRow + idof] * m_localResidual[stack.localRow + idof];
      stack.localNormalizer[0] += massNormalizer;
    }

    // step 2: volume residual

    real64 const valVol = m_localResidual[stack.localRow + m_numComponents] * m_totalDens_n[ei][0]; // we need a mass here, hence the
                                                                                                    // multiplication
    stack.localValue[1] += valVol * valVol;
    stack.localNormalizer[1] += massNormalizer;

    // step 3: energy residual

    stack.localValue[2] += m_localResidual[stack.localRow + m_numComponents + 1] * m_localResidual[stack.localRow + m_numComponents + 1];
    stack.localNormalizer[2] += energyNormalizer;
  }

protected:

  /// Number of fluid components
  integer const m_numComponents;

  /// Number of fluid phases
  integer const m_numPhases;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  arrayView2d< real64 const > const m_porosity_n;

  /// View on phase properties at the previous converged time step
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac_n;
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const m_totalDens_n;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens_n;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseInternalEnergy_n;

  /// View on solid properties at the previous converged time step
  arrayView2d< real64 const > const m_solidInternalEnergy_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
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
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   constitutive::SolidInternalEnergy const & solidInternalEnergy,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[3],
                   real64 (& residualNormalizer)[3] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               numComps, numPhases, subRegion, fluid, solid, solidInternalEnergy, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }

  }
};


} // namespace thermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALRESIDUALNORMKERNEL_HPP
