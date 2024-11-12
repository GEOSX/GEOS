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
 * @file FluxComputeKernelBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUXCOMPUTEKERNELBASE_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUXCOMPUTEKERNELBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxKernelsHelper.hpp"

namespace geos
{

namespace singlePhaseFVMKernels
{

/******************************** FluxComputeKernelBase ********************************/

/**
 * @brief Base class for FluxComputeKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of dofs).
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

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::pressure_n,
                      fields::flow::gravityCoefficient,
                      fields::flow::mobility,
                      fields::flow::dMobility_dPressure >;

  using SinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure >;

  using SlurryFluidAccessors =
    StencilMaterialAccessors< constitutive::SlurryFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  using ProppantPermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure,
                              fields::permeability::dPerm_dDispJump,
                              fields::permeability::permeabilityMultiplier >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] singleFlowAccessors accessor for wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors accessor for wrappers registered by the singlefluid model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FluxComputeKernelBase( globalIndex const rankOffset,
                         DofNumberAccessor const & dofNumberAccessor,
                         SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                         SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                         PermeabilityAccessors const & permeabilityAccessors,
                         real64 const & dt,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs )
    : m_rankOffset( rankOffset ),
    m_dt( dt ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
    m_ghostRank( singlePhaseFlowAccessors.get( fields::ghostRank {} ) ),
    m_gravCoef( singlePhaseFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
    m_pres( singlePhaseFlowAccessors.get( fields::flow::pressure {} ) ),
    m_mob( singlePhaseFlowAccessors.get( fields::flow::mobility {} ) ),
    m_dMob_dPres( singlePhaseFlowAccessors.get( fields::flow::dMobility_dPressure {} ) ),
    m_dens( singlePhaseFluidAccessors.get( fields::singlefluid::density {} ) ),
    m_dDens_dPres( singlePhaseFluidAccessors.get( fields::singlefluid::dDensity_dPressure {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary and secondary variables
  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  /// Views on fluid mobility
  ElementViewConst< arrayView1d< real64 const > > const m_mob;
  ElementViewConst< arrayView1d< real64 const > > const m_dMob_dPres;

  /// Views on fluid density
  ElementViewConst< arrayView2d< real64 const > > const m_dens;
  ElementViewConst< arrayView2d< real64 const > > const m_dDens_dPres;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;
};

} // namespace singlePhaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUXCOMPUTEKERNELBASE_HPP
