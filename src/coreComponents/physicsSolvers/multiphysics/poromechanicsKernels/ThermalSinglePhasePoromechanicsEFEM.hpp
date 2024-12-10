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
 * @file ThermalSinglePhasePoromechanicsEFEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_HPP_

#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM.hpp"

namespace geos
{

namespace thermoPoromechanicsEFEMKernels
{
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ThermalSinglePhasePoromechanicsEFEM :
  public poromechanicsEFEMKernels::SinglePhasePoromechanicsEFEM< SUBREGION_TYPE,
                                                                 CONSTITUTIVE_TYPE,
                                                                 FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = poromechanicsEFEMKernels::SinglePhasePoromechanicsEFEM< SUBREGION_TYPE,
                                                                       CONSTITUTIVE_TYPE,
                                                                       FE_TYPE >;

  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_fracturePresDofNumber;
  using Base::m_matrixPresDofNumber;
  using Base::m_wDofNumber;
  using Base::m_fluidDensity;
  using Base::m_fluidDensity_n;
  using Base::m_dFluidDensity_dPressure;
  using Base::m_porosity_n;
  using Base::m_surfaceArea;
  using Base::m_elementVolumeFrac;
  using Base::m_deltaVolume;
  using Base::m_cellsToEmbeddedSurfaces;
  using Base::m_dt;



  ThermalSinglePhasePoromechanicsEFEM( NodeManager const & nodeManager,
                                       EdgeManager const & edgeManager,
                                       FaceManager const & faceManager,
                                       localIndex const targetRegionIndex,
                                       SUBREGION_TYPE const & elementSubRegion,
                                       FE_TYPE const & finiteElementSpace,
                                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                                       EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion,
                                       arrayView1d< globalIndex const > const dispDofNumber,
                                       arrayView1d< globalIndex const > const jumpDofNumber,
                                       string const inputFlowDofKey,
                                       globalIndex const rankOffset,
                                       CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                       arrayView1d< real64 > const inputRhs,
                                       real64 const inputDt,
                                       real64 const (&inputGravityVector)[3],
                                       string const fluidModelKey );

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geos::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;


    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dFluidMassIncrement_dTemperature( 0.0 ),
      energyIncrement( 0.0 ),
      dEnergyIncrement_dJump( 0.0 ),
      dEnergyIncrement_dPressure( 0.0 ),
      dEnergyIncrement_dTemperature( 0.0 ),
      localKwTm{ 0.0 }
    {}

    /// Derivative of fluid mass accumulation wrt temperature
    real64 dFluidMassIncrement_dTemperature{};
    /// Energy accumulation
    real64 energyIncrement{};
    /// Derivative of energy accumulation wrt normal jump
    real64 dEnergyIncrement_dJump{};
    /// Derivative of energy accumulation wrt pressure
    real64 dEnergyIncrement_dPressure{};
    /// Derivative of energy accumulation wrt temperature
    real64 dEnergyIncrement_dTemperature{};
    /// C-array storage for the element local KwTm matrix.
    real64 localKwTm[numWdofs]{};
  };
  //*****************************************************************************


  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );
  //END_kernelLauncher


  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

private:

  /// Views on fluid density derivative wrt temperature
  arrayView2d< real64 const > const m_dFluidDensity_dTemperature;

  /// Views on fluid internal energy
  arrayView2d< real64 const > const m_fluidInternalEnergy_n;
  arrayView2d< real64 const > const m_fluidInternalEnergy;
  arrayView2d< real64 const > const m_dFluidInternalEnergy_dPressure;
  arrayView2d< real64 const > const m_dFluidInternalEnergy_dTemperature;

  /// Views on temperature
  arrayView1d< real64 const > const m_temperature_n;
  arrayView1d< real64 const > const m_temperature;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_matrixTemperature;
};


using ThermalSinglePhasePoromechanicsEFEMKernelFactory = finiteElement::KernelFactory< ThermalSinglePhasePoromechanicsEFEM,
                                                                                       EmbeddedSurfaceSubRegion const &,
                                                                                       arrayView1d< globalIndex const > const,
                                                                                       arrayView1d< globalIndex const > const,
                                                                                       string const,
                                                                                       globalIndex const,
                                                                                       CRSMatrixView< real64, globalIndex const > const,
                                                                                       arrayView1d< real64 > const,
                                                                                       real64 const,
                                                                                       real64 const (&)[3],
                                                                                       string const >;

} // namespace thermoPoromechanicsEFEMKernels

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_HPP_
