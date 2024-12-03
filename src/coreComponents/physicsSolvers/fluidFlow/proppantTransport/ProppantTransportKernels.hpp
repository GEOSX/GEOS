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
 * @file ProppantTransportKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORTKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORTKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/ParticleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/ParticleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"

namespace geos
{

namespace proppantTransportKernels
{
/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void launch( FLUID_WRAPPER const & fluidWrapper,
                      arrayView1d< real64 const > const & pres,
                      arrayView2d< real64 const > const & componentConcentration )
  {
    forAll< parallelDevicePolicy<> >( fluidWrapper.numElems(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const NC = fluidWrapper.numFluidComponents();
      stackArray1d< real64, constitutive::SlurryFluidBase::MAX_NUM_COMPONENTS > compConc( NC );

      for( localIndex c = 0; c < NC; ++c )
      {
        compConc[c] = componentConcentration[a][c];
      }

      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.updateFluidProperty( a, q,
                                          pres[a],
                                          compConc,
                                          0.0 );
      }
    } );
  }
};

/******************************** ComponentDensityUpdateKernel ********************************/

struct ComponentDensityUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void launch( FLUID_WRAPPER const & fluidWrapper,
                      arrayView1d< real64 const > const & pres,
                      arrayView2d< real64 const > const & componentConcentration )
  {
    forAll< parallelDevicePolicy<> >( fluidWrapper.numElems(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const NC = fluidWrapper.numFluidComponents();
      stackArray1d< real64, constitutive::SlurryFluidBase::MAX_NUM_COMPONENTS > compConc( NC );

      for( localIndex c = 0; c < NC; ++c )
      {
        compConc[c] = componentConcentration[a][c];
      }

      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.updateComponentDensity( a, q,
                                             pres[a],
                                             compConc );
      }
    } );
  }
};

/******************************** ProppantUpdateKernel ********************************/

struct ProppantUpdateKernel
{
  template< typename PROPPANT_WRAPPER >
  static void launch( PROPPANT_WRAPPER const & proppantWrapper,
                      arrayView1d< real64 const > const & proppantConc,
                      arrayView2d< real64 const > const & fluidDens,
                      arrayView2d< real64 const > const & dFluidDens_dPres,
                      arrayView3d< real64 const > const & dFluidDens_dCompConc,
                      arrayView2d< real64 const > const & fluidVisc,
                      arrayView2d< real64 const > const & dFluidVisc_dPres,
                      arrayView3d< real64 const > const & dFluidVisc_dCompConc )
  {

    forAll< parallelDevicePolicy<> >( proppantWrapper.numElems(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      proppantWrapper.update( a,
                              proppantConc[a],
                              fluidDens[a][0],
                              dFluidDens_dPres[a][0],
                              dFluidDens_dCompConc[a][0],
                              fluidVisc[a][0],
                              dFluidVisc_dPres[a][0],
                              dFluidVisc_dCompConc[a][0] );
    } );
  }
};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{
  GEOS_HOST_DEVICE
  inline
  static
  void
  compute( localIndex const NC,
           real64 const proppantConc_n,
           real64 const proppantConcNew,
           arraySlice1d< real64 const > const & componentDens_n,
           arraySlice1d< real64 const > const & componentDensNew,
           arraySlice1d< real64 const > const & GEOS_UNUSED_PARAM( dCompDens_dPres ),
           arraySlice2d< real64 const > const & dCompDens_dCompConc,
           real64 const volume,
           real64 const packPoreVolume,
           real64 const proppantLiftVolume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian );

  static void
  launch( localIndex const size,
          localIndex const NC,
          localIndex const NDOF,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & proppantConc_n,
          arrayView1d< real64 const > const & proppantConc,
          arrayView2d< real64 const > const & componentDens_n,
          arrayView3d< real64 const > const & componentDens,
          arrayView3d< real64 const > const & dCompDens_dPres,
          arrayView4d< real64 const > const & dCompDens_dCompConc,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & proppantPackVolFrac,
          arrayView1d< real64 const > const & proppantLiftFlux,
          real64 const dt,
          real64 const maxProppantConcentration,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

/******************************** FluxKernel ********************************/

struct FluxKernel
{
  using FlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::elementAperture,
                      fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::proppant::proppantConcentration,
                      fields::proppant::isProppantMobile >;

  using CellBasedFluxFlowAccessors =
    StencilAccessors< fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::elementAperture >;

  using ParticleFluidAccessors =
    StencilMaterialAccessors< constitutive::ParticleFluidBase,
                              fields::particlefluid::settlingFactor,
                              fields::particlefluid::dSettlingFactor_dPressure,
                              fields::particlefluid::dSettlingFactor_dProppantConcentration,
                              fields::particlefluid::dSettlingFactor_dComponentConcentration,
                              fields::particlefluid::collisionFactor,
                              fields::particlefluid::dCollisionFactor_dProppantConcentration >;

  using SlurryFluidAccessors =
    StencilMaterialAccessors< constitutive::SlurryFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure,
                              fields::slurryfluid::dDensity_dProppantConcentration,
                              fields::slurryfluid::dDensity_dComponentConcentration,
                              fields::singlefluid::viscosity,
                              fields::singlefluid::dViscosity_dPressure,
                              fields::slurryfluid::dViscosity_dProppantConcentration,
                              fields::slurryfluid::dViscosity_dComponentConcentration,
                              fields::slurryfluid::componentDensity,
                              fields::slurryfluid::dComponentDensity_dPressure,
                              fields::slurryfluid::dComponentDensity_dComponentConcentration,
                              fields::slurryfluid::fluidDensity,
                              fields::slurryfluid::dFluidDensity_dPressure,
                              fields::slurryfluid::dFluidDensity_dComponentConcentration >;

  using CellBasedFluxSlurryFluidAccessors =
    StencilMaterialAccessors< constitutive::SlurryFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::viscosity >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::permeabilityMultiplier >;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = typename ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementView< VIEWTYPE >;

  static void
  launch( SurfaceElementStencilWrapper const & stencilWrapper,
          localIndex const numDofPerCell,
          real64 const dt,
          globalIndex const rankOffset,
          integer const updateProppantPacking,
          R1Tensor const & unitGravityVector,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & proppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & componentDens,
          ElementViewConst< arrayView3d< real64 const > > const & dComponentDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dComponentDens_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dDens_dComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & visc,
          ElementViewConst< arrayView2d< real64 const > > const & dVisc_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & dVisc_dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dVisc_dComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & fluidDensity,
          ElementViewConst< arrayView2d< real64 const > > const & dFluidDens_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dFluidDens_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & settlingFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dProppantConc,
          ElementViewConst< arrayView2d< real64 const > > const & dSettlingFactor_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & collisionFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dCollisionFactor_dProppantConc,
          ElementViewConst< arrayView1d< integer const > > const & isProppantMobile,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
          ElementViewConst< arrayView1d< real64 const > > const & aperture,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

  static void
  launchCellBasedFluxCalculation( SurfaceElementStencilWrapper const & stencilWrapper,
                                  R1Tensor const & unitGravityVector,
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                  ElementViewConst< arrayView2d< real64 const > > const & dens,
                                  ElementViewConst< arrayView2d< real64 const > > const & visc,
                                  ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                  ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
                                  ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                  ElementView< arrayView2d< real64 > > const & cellBasedFlux );

  /**
   * @brief Compute flux and its derivatives for a given multi-element connector.
   *
   * This is a specialized version that flux in a single region, and uses
   * element pairing instead of a proper junction.
   */
  template< localIndex MAX_NUM_FLUX_ELEMS >
  GEOS_HOST_DEVICE
  static void
  computeJunction( localIndex const numElems,
                   localIndex const numDofPerCell,
                   arraySlice1d< localIndex const > const & stencilElementIndices,
                   arrayView1d< real64 const > const & pres,
                   arrayView1d< real64 const > const & proppantConc,
                   arrayView3d< real64 const > const & componentDens,
                   arrayView3d< real64 const > const & dComponentDens_dPres,
                   arrayView4d< real64 const > const & dComponentDens_dComponentConc,
                   arrayView1d< real64 const > const & gravDepth,
                   arrayView2d< real64 const > const & dens,
                   arrayView2d< real64 const > const & dDens_dPres,
                   arrayView2d< real64 const > const & dDens_dProppantConc,
                   arrayView3d< real64 const > const & dDens_dComponentConc,
                   arrayView2d< real64 const > const & visc,
                   arrayView2d< real64 const > const & dVisc_dPres,
                   arrayView2d< real64 const > const & dVisc_dProppantConc,
                   arrayView3d< real64 const > const & dVisc_dComponentConc,
                   arrayView2d< real64 const > const & fluidDensity,
                   arrayView2d< real64 const > const & dFluidDens_dPres,
                   arrayView3d< real64 const > const & dFluidDens_dComponentConc,
                   arrayView1d< real64 const > const & settlingFactor,
                   arrayView1d< real64 const > const & dSettlingFactor_dPres,
                   arrayView1d< real64 const > const & dSettlingFactor_dProppantConc,
                   arrayView2d< real64 const > const & dSettlingFactor_dComponentConc,
                   arrayView1d< real64 const > const & collisionFactor,
                   arrayView1d< real64 const > const & dCollisionFactor_dProppantConc,
                   arrayView1d< integer const > const & isProppantMobile,
                   real64 const (&transmissibility)[MAX_NUM_FLUX_ELEMS],
                   real64 const (&apertureWeight)[MAX_NUM_FLUX_ELEMS],
                   real64 const (&geometricWeight)[MAX_NUM_FLUX_ELEMS],
                   real64 const dt,
                   arraySlice1d< real64 > const & localFlux,
                   arraySlice2d< real64 > const & localFluxJacobian );


  template< localIndex MAX_NUM_FLUX_ELEMS >
  GEOS_HOST_DEVICE
  static void
  computeCellBasedFlux( localIndex const numElems,
                        arraySlice1d< localIndex const > const & stencilElementIndices,
                        arrayView1d< real64 const > const & pres_n,
                        arrayView1d< real64 const > const & gravDepth,
                        arrayView2d< real64 const > const & dens,
                        arrayView2d< real64 const > const & visc,
                        arraySlice1d< R1Tensor const > const & cellCenterToEdgeCenters,
                        real64 const (&transmissibility)[MAX_NUM_FLUX_ELEMS],
                        real64 const (&apertureWeight)[MAX_NUM_FLUX_ELEMS],
                        real64 const (&geometricWeight)[MAX_NUM_FLUX_ELEMS],
                        arrayView2d< real64 > const & cellBasedFlux );
};

struct ProppantPackVolumeKernel
{

  using FlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::elementAperture,
                      fields::elementVolume,
                      fields::proppant::cellBasedFlux,
                      fields::proppant::isProppantMobile,
                      fields::proppant::isProppantBoundary >;

  using ParticleFluidAccessors =
    StencilMaterialAccessors< constitutive::ParticleFluidBase, fields::particlefluid::settlingFactor >;

  using SlurryFluidAccessors =
    StencilMaterialAccessors< constitutive::SlurryFluidBase,
                              fields::singlefluid::density,
                              fields::slurryfluid::fluidDensity,
                              fields::slurryfluid::fluidViscosity >;


  template< typename VIEWTYPE >
  using ElementView = ElementRegionManager::ElementView< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  static void
  launchProppantPackVolumeCalculation( SurfaceElementStencil const & stencil,
                                       real64 const dt,
                                       real64 const proppantDensity,
                                       real64 const proppantDiameter,
                                       real64 const maxProppantConcentration,
                                       R1Tensor const & unitGravityVector,
                                       real64 const criticalShieldsNumber,
                                       real64 const fricitonCoefficient,
                                       ElementView< arrayView1d< real64 const > > const & settlingFactor,
                                       ElementView< arrayView2d< real64 const > > const & density,
                                       ElementView< arrayView2d< real64 const > > const & fluidDensity,
                                       ElementView< arrayView2d< real64 const > > const & fluidViscosity,
                                       ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                       ElementView< arrayView1d< integer const > > const & isProppantBoundaryElement,
                                       ElementView< arrayView1d< real64 const > > const & aperture,
                                       ElementView< arrayView1d< real64 const > > const & volume,
                                       ElementView< arrayView1d< integer const > > const & elemGhostRank,
                                       ElementView< arrayView2d< real64 const > > const & cellBasedFlux,
                                       ElementView< arrayView1d< real64 > > const & conc,
                                       ElementView< arrayView1d< real64 > > const & proppantPackVolFrac,
                                       ElementView< arrayView1d< real64 > > const & proppantExcessPackVolume,
                                       ElementView< arrayView1d< real64 > > const & proppantLiftFlux );

  static void
  launchProppantPackVolumeUpdate( SurfaceElementStencil const & stencil,
                                  R1Tensor const & unitGravityVector,
                                  real64 const maxProppantConcentration,
                                  ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                  ElementView< arrayView1d< real64 const > > const & proppantExcessPackVolume,
                                  ElementView< arrayView1d< real64 > > const & conc,
                                  ElementView< arrayView1d< real64 > > const & proppantPackVolFrac );

  GEOS_HOST_DEVICE
  inline
  static
  void
  computeProppantPackVolume( localIndex const numElems,
                             real64 const dt,
                             real64 const proppantDensity,
                             real64 const proppantDiameter,
                             real64 const maxProppantConcentration,
                             R1Tensor const & unitGravityVector,
                             real64 const criticalShieldsNumber,
                             real64 const frictionCoefficient,
                             arraySlice1d< localIndex const > const & stencilElementIndices,
                             arraySlice1d< real64 const > const & stencilWeights,
                             arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                             arrayView1d< real64 const > const & settlingFactor,
                             arrayView2d< real64 const > const & density,
                             arrayView2d< real64 const > const & fluidDensity,
                             arrayView2d< real64 const > const &,
                             arrayView1d< real64 const > const & volume,
                             arrayView1d< real64 const > const & aperture,
                             arrayView1d< integer const > const & elemGhostRank,
                             arrayView1d< integer const > const & isProppantBoundaryElement,
                             arrayView1d< integer const > const & isProppantMobile,
                             arrayView2d< real64 const > const & cellBasedFlux,
                             arrayView1d< real64 > const & conc,
                             arrayView1d< real64 > const & proppantPackVolFrac,
                             arrayView1d< real64 > const & proppantExcessPackVolume,
                             arrayView1d< real64 > const & proppantLiftFlux );

  GEOS_HOST_DEVICE
  inline
  static
  void
  updateProppantPackVolume( localIndex const numElems,
                            arraySlice1d< localIndex const > const & stencilElementIndices,
                            arraySlice1d< real64 const > const & stencilWeights,
                            arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                            R1Tensor const & unitGravityVector,
                            real64 const maxProppantConcentration,
                            arrayView1d< integer const > const & isProppantMobile,
                            arrayView1d< real64 const > const & proppantExcessPackVolume,
                            arrayView1d< real64 > const & conc,
                            arrayView1d< real64 > const & proppantPackVolFrac );
};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numComp,
                      ElementSubRegionBase const & subRegion,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numComp( numComp ),
    m_volume( subRegion.getElementVolume() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 const normalizer = LvArray::math::max( m_minNormalizer, m_volume[ei] );
    for( integer idof = 0; idof < m_numComp; ++idof )
    {
      real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / normalizer;
      if( valMass > stack.localValue[0] )
      {
        stack.localValue[0] = valMass;
      }
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    real64 const normalizer = LvArray::math::max( m_minNormalizer, m_volume[ei] );
    for( integer idof = 0; idof < m_numComp; ++idof )
    {
      stack.localValue[0] += m_localResidual[stack.localRow + idof] * m_localResidual[stack.localRow + idof];
      stack.localNormalizer[0] += normalizer;
    }
  }


protected:

  /// Number of fluid components
  integer const m_numComp;

  /// View on the element volume
  arrayView1d< real64 const > const m_volume;
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
   * @param[in] numComp the number of components
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   integer const numComp,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               numComp, subRegion, minNormalizer );
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



} // namespace proppantTransportKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORTKERNELS_HPP_
