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
 * @file SolidMechanicsConformingContactKernelsHelper.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSCONFORMINGCONTACTKERNELSHELPER_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSCONFORMINGCONTACTKERNELSHELPER_HPP_

#include "common/GeosxMacros.hpp"

namespace geos
{

namespace solidMechanicsConformingContactKernelsHelper
{

template< int I_SIZE,
          int J_SIZE,
          int NUM_NODES >
GEOS_HOST_DEVICE
inline
void accumulateAtuLocalOperator( real64 ( & matrix )[I_SIZE][J_SIZE],
                                 real64 ( & N )[NUM_NODES],
                                 int const ( &perm )[NUM_NODES],
                                 real64 const detJ )
{
  for( int a=0; a < NUM_NODES; ++a )
  {
    for( int i=0; i < I_SIZE; ++i )
    {
      matrix[i][ a*3 + i ] += -N[ perm[ a ] ] * detJ;
      matrix[i][ 3*NUM_NODES + a*3 + i ] += N[ perm[ a ] ] * detJ;
    }
  }
}

template< int I_SIZE,
          int J_SIZE,
          int NUM_SUPPORTS >
GEOS_HOST_DEVICE
inline
void assembleStrainOperator( real64 ( & strainMatrix )[I_SIZE][J_SIZE],
                             real64 ( & dNdX )[NUM_SUPPORTS][3] )
{
  LvArray::tensorOps::fill< I_SIZE, J_SIZE >( strainMatrix, 0 );  //make 0
  for( int a=0; a < NUM_SUPPORTS; ++a )
  {

    strainMatrix[0][a*3 + 0] = dNdX[a][0];
    strainMatrix[1][a*3 + 1] = dNdX[a][1];
    strainMatrix[2][a*3 + 2] = dNdX[a][2];

    strainMatrix[3][a*3 + 1] = dNdX[a][2];
    strainMatrix[3][a*3 + 2] = dNdX[a][1];

    strainMatrix[4][a*3 + 0] = dNdX[a][2];
    strainMatrix[4][a*3 + 2] = dNdX[a][0];

    strainMatrix[5][a*3 + 0] = dNdX[a][1];
    strainMatrix[5][a*3 + 1] = dNdX[a][0];
  }
}

} // solidMechanicsConformingContactKernelsHelper

} // geos

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSCONFORMINGCONTACTKERNELSHELPER_HPP_ */
