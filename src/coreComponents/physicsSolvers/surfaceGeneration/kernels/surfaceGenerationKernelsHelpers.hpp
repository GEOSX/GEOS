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
 * @file surfaceGenerationKernelsHelpers.hpp
 */


#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

namespace geos
{

namespace surfaceGenerationKernelsHelpers
{

GEOS_HOST_DEVICE
inline void computeNodalForce( real64 const ( &stress) [ 6 ],
                               real64 const ( &dNdX) [ 3 ],
                               real64 const detJ,
                               real64 ( & force ) [ 3 ] )
{

  force[ 0 ] -= ( stress[ 0 ] * dNdX[ 0 ] +
                  stress[ 5 ] * dNdX[ 1 ] +
                  stress[ 4 ] * dNdX[ 2 ] ) * detJ;
  force[ 1 ] -= ( stress[ 5 ] * dNdX[ 0 ] +
                  stress[ 1 ] * dNdX[ 1 ] +
                  stress[ 3 ] * dNdX[ 2 ] ) * detJ;
  force[ 2 ] -= ( stress[ 4 ] * dNdX[ 0 ] +
                  stress[ 3 ] * dNdX[ 1 ] +
                  stress[ 2 ] * dNdX[ 2 ] ) * detJ;
}

GEOS_HOST_DEVICE
inline void scaleNodalForce( real64 const bulkModulus,
                             real64 const shearModulus,
                             real64 ( & force ) [ 3 ] )
{
  real64 const YoungModulus = 9 * bulkModulus * shearModulus / ( 3 * bulkModulus + shearModulus );
  real64 const poissonRatio = ( 3 * bulkModulus - 2 * shearModulus ) / ( 2 * ( 3 * bulkModulus + shearModulus ) );

  LvArray::tensorOps::scale< 3 >( force, YoungModulus );
  LvArray::tensorOps::scale< 3 >( force, 1.0 / (1 - poissonRatio * poissonRatio) );
}


}

}
