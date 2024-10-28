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
 * @file EFEMKernelsHelper.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMKERNELSHELPER_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMKERNELSHELPER_HPP_


#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "common/GeosxMacros.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/output.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace solidMechanicsEFEMKernelsHelper
{

template< int NUM_NODES >
GEOS_HOST_DEVICE
inline
void computeHeavisideFunction( integer (& heaviside)[NUM_NODES],
                               real64 (& X)[NUM_NODES][3],
                               arraySlice1d< real64 const > const normalVector,
                               arraySlice1d< real64 const > const elementCenter )
{
  for( int a=0; a < NUM_NODES; a++ )
  {
    real64 distanceVector[3];
    LvArray::tensorOps::copy< 3 >( distanceVector, X[a] );
    LvArray::tensorOps::subtract< 3 >( distanceVector, elementCenter );

    heaviside[a] = LvArray::tensorOps::AiBi< 3 >( distanceVector, normalVector ) > 0 ? 1 : 0;
  }

}



template< int I_SIZE,
          int J_SIZE,
          int NUM_NODES >
GEOS_HOST_DEVICE
inline
void assembleStrainOperator( real64 ( & strainMatrix )[I_SIZE][J_SIZE],
                             real64 ( & dNdX )[NUM_NODES][3] )
{
  LvArray::tensorOps::fill< I_SIZE, J_SIZE >( strainMatrix, 0 );  //make 0
  for( int a=0; a < NUM_NODES; ++a )
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

template< int NUM_NODES >
GEOS_HOST_DEVICE
inline
void assembleCompatibilityOperator( real64 ( & compMatrix )[3][6],
                                    arraySlice1d< real64 const > const & nVec,
                                    arraySlice1d< real64 const > const & tVec1,
                                    arraySlice1d< real64 const > const & tVec2,
                                    integer ( & heavisideFun )[NUM_NODES],
                                    real64 ( & dNdX )[NUM_NODES][3] )
{
  //GEOS_MARK_FUNCTION;

  // Fill in compatibility operator

  // 1. construct mvector sum(dNdX(a) * H(a)) value for each Gauss point
  real64 mVec[3];
  LvArray::tensorOps::fill< 3 >( mVec, 0 );
  for( integer a=0; a<NUM_NODES; ++a )
  {
    mVec[0] -= dNdX[a][0] * heavisideFun[a];
    mVec[1] -= dNdX[a][1] * heavisideFun[a];
    mVec[2] -= dNdX[a][2] * heavisideFun[a];
  }

  // 2. fill in the operator itself
  // sym(n dyadic m), sym(m dyadic t1) and sym (m dyadic t2)
  LvArray::tensorOps::symRij_eq_AiBj_plus_AjBi< 3 >( compMatrix[0], mVec, nVec );
  LvArray::tensorOps::symRij_eq_AiBj_plus_AjBi< 3 >( compMatrix[1], mVec, tVec1 );
  LvArray::tensorOps::symRij_eq_AiBj_plus_AjBi< 3 >( compMatrix[2], mVec, tVec2 );

  // scale by 0.5 the diagonal entries (it's like a strain)
  for( int ii = 0; ii<3; ii++ )
  {
    for( int jj=0; jj<3; jj++ )
    {
      compMatrix[ii][jj] *= 0.5;
    }
  }

}

GEOS_HOST_DEVICE
inline
void assembleEquilibriumOperator( real64 ( & eqMatrix )[3][6],
                                  arraySlice1d< real64 const > const nVec,
                                  arraySlice1d< real64 const > const tVec1,
                                  arraySlice1d< real64 const > const tVec2,
                                  real64 const hInv )
{
  // (n dyadic n), sym(n dyadic t1) and sym (n dyadic t2)
  LvArray::tensorOps::symRij_eq_AiBj_plus_AjBi< 3 >( eqMatrix[0], nVec, nVec );
  LvArray::tensorOps::symRij_eq_AiBj_plus_AjBi< 3 >( eqMatrix[1], nVec, tVec1 );
  LvArray::tensorOps::symRij_eq_AiBj_plus_AjBi< 3 >( eqMatrix[2], nVec, tVec2 );

  for( int ii = 0; ii<3; ii++ )
  {
    for( int jj=0; jj<3; jj++ )
    {
      eqMatrix[ii][jj] *= 0.5;
    }
  }

  LvArray::tensorOps::scale< 3, 6 >( eqMatrix, -hInv );
}

}
} // geos


#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_CONTACT_EFEMKERNELSHELPER_HPP_ */
