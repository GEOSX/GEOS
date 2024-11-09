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
 * @file FluxComputeKernelBase.cpp
 */

#include "physicsSolvers/fluidFlow/kernels/compositional/FluxComputeKernelBase.hpp"

#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/SurfaceElementStencil.hpp"
#include "finiteVolume/EmbeddedSurfaceToCellStencil.hpp"
#include "finiteVolume/FaceElementToCellStencil.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"

namespace geos
{
using namespace constitutive;

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** FluxComputeKernelBase ********************************/

FluxComputeKernelBase::FluxComputeKernelBase( integer const numPhases,
                                              globalIndex const rankOffset,
                                              DofNumberAccessor const & dofNumberAccessor,
                                              CompFlowAccessors const & compFlowAccessors,
                                              MultiFluidAccessors const & multiFluidAccessors,
                                              real64 const dt,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs,
                                              BitFlags< KernelFlags > kernelFlags )
  : m_numPhases( numPhases ),
  m_rankOffset( rankOffset ),
  m_dt( dt ),
  m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
  m_ghostRank( compFlowAccessors.get( fields::ghostRank {} ) ),
  m_gravCoef( compFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
  m_pres( compFlowAccessors.get( fields::flow::pressure {} ) ),
  m_phaseVolFrac( compFlowAccessors.get( fields::flow::phaseVolumeFraction {} ) ),
  m_dPhaseVolFrac( compFlowAccessors.get( fields::flow::dPhaseVolumeFraction {} ) ),
  m_dCompFrac_dCompDens( compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity {} ) ),
  m_phaseCompFrac( multiFluidAccessors.get( fields::multifluid::phaseCompFraction {} ) ),
  m_dPhaseCompFrac( multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction {} ) ),
  m_localMatrix( localMatrix ),
  m_localRhs( localRhs ),
  m_kernelFlags( kernelFlags )
{}

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos
