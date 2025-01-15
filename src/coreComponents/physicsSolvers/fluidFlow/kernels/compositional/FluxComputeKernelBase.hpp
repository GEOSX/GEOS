/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FluxComputeKernelBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_FLUXCOMPUTEKERNELBASE_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_FLUXCOMPUTEKERNELBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

enum class KernelFlags
{
  /// Flag to specify whether capillary pressure is used or not
  CapPressure = 1 << 0, // 1
  /// Flag to specify whether diffusion is used or not
  Diffusion = 1 << 1, // 2
  /// Flag to specify whether dispersion is used or not
  Dispersion = 1 << 2, // 4
  /// Flag indicating whether total mass equation is formed or not
  TotalMassEquation = 1 << 3, // 8
  /// Flag indicating whether gravity treatment is checking phase presence or not
  CheckPhasePresenceInGravity = 1 << 4, // 16
  /// Flag indicating whether C1-PPU is used or not
  C1PPU = 1 << 5, // 32
  /// Flag indicating whether IHU is used or not
  IHU = 1 << 6, // 64
  /// Flag indicating whether HU 2-phase simplified version is used or not
  HU2PH = 1 << 7 // 128
};

/******************************** FluxComputeKernelBase ********************************/

/**
 * @brief Base class for FluxComputeKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of components/dofs).
 */
class FluxComputeKernelBase
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using CompFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::gravityCoefficient,
                      fields::flow::pressure,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity,
                      fields::flow::phaseVolumeFraction,
                      fields::flow::dPhaseVolumeFraction,
                      fields::flow::phaseMobility,
                      fields::flow::dPhaseMobility >;
  using MultiFluidAccessors =
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::dPhaseMassDensity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  using CapPressureAccessors =
    StencilMaterialAccessors< constitutive::CapillaryPressureBase,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] compFlowAccessors accessor for wrappers registered by the solver
   * @param[in] multiFluidAccessors accessor for wrappers registered by the multifluid model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed all together
   */
  FluxComputeKernelBase( integer const numPhases,
                         globalIndex const rankOffset,
                         DofNumberAccessor const & dofNumberAccessor,
                         CompFlowAccessors const & compFlowAccessors,
                         MultiFluidAccessors const & multiFluidAccessors,
                         real64 const dt,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         BitFlags< KernelFlags > kernelFlags );

protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary and secondary variables

  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  /// Views on phase volume fractions
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseVolFrac;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseVolFrac;

  /// Views on derivatives of comp fractions
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens;

  /// Views on phase component fractions
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const m_phaseCompFrac;
  ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const m_dPhaseCompFrac;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  BitFlags< KernelFlags > const m_kernelFlags;
};

namespace helpers
{
template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

template< localIndex numComp, localIndex numFluxSupportPoints >
GEOS_HOST_DEVICE
static void calculateMeanDensity( localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  integer const checkPhasePresenceInGravity,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  real64 & densMean, real64 (& dDensMean_dPres)[numFluxSupportPoints], real64 (& dDensMean_dComp)[numFluxSupportPoints][numComp] )
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  densMean = 0;
  integer denom = 0;
  real64 dDens_dC[numComp]{};
  for( localIndex i = 0; i < numFluxSupportPoints; ++i )
  {
    localIndex const er = seri[i];
    localIndex const esr = sesri[i];
    localIndex const ei = sei[i];

    bool const phaseExists = (phaseVolFrac[er][esr][ei][ip] > 0);
    if( checkPhasePresenceInGravity && !phaseExists )
    {
      dDensMean_dPres[i] = 0.0;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dComp[i][jc] = 0.0;
      }
      continue;
    }

    // density
    real64 const density = phaseMassDens[er][esr][ei][0][ip];
    real64 const dDens_dPres = dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];

    applyChainRule( numComp,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseMassDens[er][esr][ei][0][ip],
                    dDens_dC,
                    Deriv::dC );

    // average density and derivatives
    densMean += density;
    dDensMean_dPres[i] = dDens_dPres;
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dDensMean_dComp[i][jc] = dDens_dC[jc];
    }
    denom++;
  }
  if( denom > 1 )
  {
    densMean /= denom;
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      dDensMean_dPres[i] /= denom;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dComp[i][jc] /= denom;
      }
    }
  }
}

}

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_FLUXCOMPUTEKERNELBASE_HPP
