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
 * @file SolidModelDiscretizationOps.hpp
 */


#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPS_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPS_HPP_

#include "common/DataTypes.hpp"

namespace geos
{
namespace constitutive
{

/**
 * Base class for objects that assist in performing operations for a
 * discretization method using constitutive data.
 * Derived objects will define specific implementations of methods
 * appropriately for the corresponding constitutive relation.
 */
struct SolidModelDiscretizationOps
{
  /**
   * @brief Compute inner product matrix for solid mechanics
   * @tparam NUM_SUPPORT POINTS Number of support points (nodes) for this element
   * @tparam BASIS_GRADIENT Finite element shape function gradients type
   * @tparam CBF Callback function
   * @param gradN Finite Element shape function gradients
   * @param elementStiffness Local stiffness matrix
   * @param callbackFunction The callback function
   */
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT,
            typename CBF >
  GEOS_HOST_DEVICE
  void BTDB( BASIS_GRADIENT const & gradN,
             real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3],
             CBF && callbackFunction );

  /**
   * @brief Compute upper portion of inner product matrix for solid mechanics, assuming D is symmetric
   * @tparam NUM_SUPPORT POINTS Number of support points (nodes) for this element
   * @tparam BASIS_GRADIENT Finite element shape function gradients type
   * @tparam CBF Callback function
   * @param gradN Finite Element shape function gradients
   * @param elementStiffness Local stiffness matrix
   * @param callbackFunction The callback function
   */
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT,
            typename CBF >
  GEOS_HOST_DEVICE
  void upperBTDB( BASIS_GRADIENT const & gradN,
                  real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3],
                  CBF && callbackFunction );

  /**
   * @brief Copy upper entries to lower half of inner product matrix
   * @tparam NUM_SUPPORT_POINTS Number of support points (nodes) for this element
   * @param elementStiffness Local stiffness matrix
   */
  template< int NUM_SUPPORT_POINTS >
  GEOS_HOST_DEVICE
  static
  void fillLowerBTDB( real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  /**
   * @brief Compute diagonal of inner product matrix for solid mechanics
   * @tparam NUM_SUPPORT POINTS Number of support points (nodes) for this element
   * @tparam BASIS_GRADIENT Finite element shape function gradients type
   * @tparam CBF Callback function
   * @param gradN Finite Element shape function gradients
   * @param diagElementStiffness Local stiffness matrix diagonal
   * @param callbackFunction The callback function
   */
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT,
            typename CBF >
  GEOS_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 ( &diagElementStiffness )[NUM_SUPPORT_POINTS*3],
                 CBF && callbackFunction );

  /**
   * @brief Compute row-summed diagonal of inner product matrix for solid mechanics
   * @tparam NUM_SUPPORT POINTS Number of support points (nodes) for this element
   * @tparam BASIS_GRADIENT Finite element shape function gradients type
   * @tparam CBF Callback function
   * @param gradN Finite Element shape function gradients
   * @param diagSumElementStiffness Local stiffness matrix diagonal
   * @param callbackFunction The callback function
   */
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT,
            typename CBF >
  GEOS_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 ( &diagSumElementStiffness )[NUM_SUPPORT_POINTS*3],
                       CBF && callbackFunction );

};


/// @copydoc SolidModelDiscretizationOps::BTDB
template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT,
          typename CBF >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOps::BTDB( BASIS_GRADIENT const & gradN,
                                        real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3],
                                        CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=0; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };
      callbackFunction( a, b, gradNa_gradNb, elementStiffness );
    }
  }
}


/// @copydoc SolidModelDiscretizationOps::upperBTDB
template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT,
          typename CBF >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOps::upperBTDB( BASIS_GRADIENT const & gradN,
                                             real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3],
                                             CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };
      callbackFunction( a, b, gradNa_gradNb, elementStiffness );
    }
  }
}


/// @copydoc SolidModelDiscretizationOps::fillLowerBTDB
template< int NUM_SUPPORT_POINTS >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOps::fillLowerBTDB( real64 ( & elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
{
  for( int row=1; row<NUM_SUPPORT_POINTS*3; ++row )
  {
    for( int col=0; col<row; ++col )
    {
      elementStiffness[row][col] = elementStiffness[col][row];
    }
  }
}


/// @copydoc SolidModelDiscretizationOps::diagBTDB
template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT,
          typename CBF >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOps::diagBTDB( BASIS_GRADIENT const & gradN,
                                            real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS*3],
                                            CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    real64 const gradN_gradN[3] = { gradN[a][0] * gradN[a][0], gradN[a][1] * gradN[a][1], gradN[a][2] * gradN[a][2] };

    callbackFunction( a, gradN_gradN, diagElementStiffness );
  }
}


/// @copydoc SolidModelDiscretizationOps::diagRowSumBTDB
template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT,
          typename CBF >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOps::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                  real64 ( & diagSumElementStiffness )[NUM_SUPPORT_POINTS*3],
                                                  CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };
      callbackFunction( a, b, gradNa_gradNb, diagSumElementStiffness );
    }
  }
}

}
}


#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPS_HPP_ */
