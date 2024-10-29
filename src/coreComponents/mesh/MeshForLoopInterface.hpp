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

/**
 * @file MeshForLoopInterface.hpp
 */

#ifndef GEOS_MESH_MESHFORLOOPINTERFACE_HPP
#define GEOS_MESH_MESHFORLOOPINTERFACE_HPP

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geos
{

/**
 * @brief Loop over all elements in a geos::MeshLevel.
 * @tparam POLICY The execution policy for the loop over elements in a geos::ElementSubRegionBase.
 * @tparam LAMBDA The type of lambda function to execute for each element.
 * @param mesh The geos::MeshLevel that will have all of its elements processed by this function.
 * @param lambda The type of lambda function to execute for each element.
 */
template< class POLICY=serialPolicy, typename LAMBDA=void >
void forAllElemsInMesh( MeshLevel const & mesh, LAMBDA && lambda )
{

  ElementRegionManager const & elemManager = mesh.getElemManager();
  elemManager.forElementSubRegionsComplete< ElementSubRegionBase >( [&] ( localIndex const er,
                                                                          localIndex const esr,
                                                                          ElementRegionBase const &,
                                                                          ElementSubRegionBase const & subRegion )
  {
    forAll< POLICY >( subRegion.size(), [&]( localIndex const k ) { lambda( er, esr, k ); } );
  } );
}

/**
 * @brief @return Return the minimum location/indices for a value condition specified by @p lambda.
 * @tparam LAMBDA The type of the lambda function to be used to specify the minimum condition.
 * @param mesh The geos::MeshLevel that will have all of its elements processed by this function.
 * @param lambda  A lambda function that returns as value that will be used in the minimum comparison.
 */
template< typename LAMBDA >
auto
minLocOverElemsInMesh( MeshLevel const & mesh, LAMBDA && lambda )
{
  using NUMBER = decltype( lambda( 0, 0, 0 ) );

  NUMBER minVal = std::numeric_limits< NUMBER >::max();
  localIndex minReg = -1, minSubreg = -1, minIndex = -1;

  ElementRegionManager const & elemManager = mesh.getElemManager();

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase const & elemRegion = elemManager.getRegion( er );

    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr, CellElementSubRegion const & subRegion )
    {
      localIndex const size = subRegion.size();
      for( localIndex k = 0; k < size; ++k )
      {
        NUMBER const val = lambda( er, esr, k );
        if( val < minVal )
        {
          minVal = val;
          minReg = er;
          minSubreg = esr;
          minIndex = k;
        }
      }
    } );
  }

  return std::make_pair( minVal, std::make_tuple( minReg, minSubreg, minIndex ));
}

} // namespace geos

#endif // GEOS_MESH_MESHFORLOOPINTERFACE_HPP
