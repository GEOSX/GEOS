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
 * @file AquiferBCKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_AQUIFERBCKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_AQUIFERBCKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"

namespace geos
{

namespace singlePhaseFVMKernels
{

/******************************** AquiferBCKernel ********************************/

/**
 * @brief Functions to assemble aquifer boundary condition contributions to residual and Jacobian
 */
struct AquiferBCKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  GEOS_HOST_DEVICE
  static void
  compute( real64 const & aquiferVolFlux,
           real64 const & dAquiferVolFlux_dPres,
           real64 const & aquiferDens,
           real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & dt,
           real64 & localFlux,
           real64 & localFluxJacobian )
  {
    if( aquiferVolFlux > 0 ) // aquifer is upstream
    {
      localFlux -= dt * aquiferVolFlux * aquiferDens;
      localFluxJacobian -= dt * dAquiferVolFlux_dPres * aquiferDens;
    }
    else // reservoir is upstream
    {
      localFlux -= dt * aquiferVolFlux * dens;
      localFluxJacobian -= dt * (dAquiferVolFlux_dPres * dens + aquiferVolFlux * dDens_dPres);
    }
  }

  static void
  launch( BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const & aquiferDens,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & pres_n,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          real64 const & timeAtBeginningOfStep,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    using Order = BoundaryStencil::Order;

    BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
    BoundaryStencil::WeightContainerViewConstType const & weight = stencil.getWeights();

    forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      // working variables
      real64 localFlux = 0.0;
      real64 localFluxJacobian = 0.0;

      localIndex const er  = seri( iconn, Order::ELEM );
      localIndex const esr = sesri( iconn, Order::ELEM );
      localIndex const ei  = sefi( iconn, Order::ELEM );
      real64 const areaFraction = weight( iconn, Order::ELEM );

      // compute the aquifer influx rate using the pressure influence function and the aquifer props
      real64 dAquiferVolFlux_dPres = 0.0;
      real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                              dt,
                                                              pres[er][esr][ei],
                                                              pres_n[er][esr][ei],
                                                              gravCoef[er][esr][ei],
                                                              areaFraction,
                                                              dAquiferVolFlux_dPres );

      // compute the phase/component aquifer flux
      AquiferBCKernel::compute( aquiferVolFlux,
                                dAquiferVolFlux_dPres,
                                aquiferDens,
                                dens[er][esr][ei][0],
                                dDens_dPres[er][esr][ei][0],
                                dt,
                                localFlux,
                                localFluxJacobian );

      // Add to residual/jacobian
      if( ghostRank[er][esr][ei] < 0 )
      {
        globalIndex const globalRow = dofNumber[er][esr][ei];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux );
        localMatrix.addToRow< parallelDeviceAtomic >( localRow,
                                                      &dofNumber[er][esr][ei],
                                                      &localFluxJacobian,
                                                      1 );
      }
    } );
  }

};

} // namespace singlePhaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_AQUIFERBCKERNEL_HPP
