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
 * @file SolidModelDiscretizationOpsFullyAnisotroipic.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIAZTIONOPSFULLYANISOTROPIC_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIAZTIONOPSFULLYANISOTROPIC_HPP_

#include "SolidModelDiscretizationOps.hpp"

namespace geos
{
namespace constitutive
{


struct SolidModelDiscretizationOpsFullyAnisotroipic : public SolidModelDiscretizationOps // TODO: spelling error,
                                                                                         // convert to "General" anyway
{
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void BTDB( BASIS_GRADIENT const & gradN,
             real64 const & detJxW,
             real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void upperBTDB( BASIS_GRADIENT const & gradN,
                  real64 const & detJxW,
                  real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 const & detJxW,
                 real64 ( &diagElementStiffness )[NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 const & detJxW,
                       real64 ( &diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] );

  GEOS_HOST_DEVICE
  inline
  void scaleParams( real64 const scale )
  {
    LvArray::tensorOps::scale< 6, 6 >( m_c, scale );
  }


  real64 m_c[6][6];
};


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsFullyAnisotroipic::BTDB( BASIS_GRADIENT const & gradN,
                                                         real64 const & detJxW,
                                                         real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=0; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };


      elementStiffness[a*3+0][b*3+0] = elementStiffness[a*3+0][b*3+0] +
                                       ( m_c[0][0] * gradNa_gradNb[0][0] + m_c[0][5] * gradNa_gradNb[0][1] + m_c[0][4] * gradNa_gradNb[0][2] +
                                         m_c[5][0] * gradNa_gradNb[1][0] + m_c[5][5] * gradNa_gradNb[1][1] + m_c[5][4] * gradNa_gradNb[1][2] +
                                         m_c[4][0] * gradNa_gradNb[2][0] + m_c[4][5] * gradNa_gradNb[2][1] + m_c[4][4] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+0][b*3+1] = elementStiffness[a*3+0][b*3+1] +
                                       (m_c[0][5] * gradNa_gradNb[0][0] + m_c[0][1] * gradNa_gradNb[0][1] + m_c[0][3] * gradNa_gradNb[0][2] +
                                        m_c[5][5] * gradNa_gradNb[1][0] + m_c[5][1] * gradNa_gradNb[1][1] + m_c[5][3] * gradNa_gradNb[1][2] +
                                        m_c[4][5] * gradNa_gradNb[2][0] + m_c[4][1] * gradNa_gradNb[2][1] + m_c[4][3] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+0][b*3+2] = elementStiffness[a*3+0][b*3+2] +
                                       (m_c[0][4] * gradNa_gradNb[0][0] + m_c[0][3] * gradNa_gradNb[0][1] + m_c[0][2] * gradNa_gradNb[0][2] +
                                        m_c[5][4] * gradNa_gradNb[1][0] + m_c[5][3] * gradNa_gradNb[1][1] + m_c[5][2] * gradNa_gradNb[1][2] +
                                        m_c[4][4] * gradNa_gradNb[2][0] + m_c[4][3] * gradNa_gradNb[2][1] + m_c[4][2] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+1][b*3+0] = elementStiffness[a*3+1][b*3+0] +
                                       (m_c[5][0] * gradNa_gradNb[0][0] + m_c[5][5] * gradNa_gradNb[0][1] + m_c[5][4] * gradNa_gradNb[0][2] +
                                        m_c[1][0] * gradNa_gradNb[1][0] + m_c[1][5] * gradNa_gradNb[1][1] + m_c[1][4] * gradNa_gradNb[1][2] +
                                        m_c[3][0] * gradNa_gradNb[2][0] + m_c[3][5] * gradNa_gradNb[2][1] + m_c[3][4] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+1][b*3+1] = elementStiffness[a*3+1][b*3+1] +
                                       (m_c[5][5] * gradNa_gradNb[0][0] + m_c[5][1] * gradNa_gradNb[0][1] + m_c[5][3] * gradNa_gradNb[0][2] +
                                        m_c[1][5] * gradNa_gradNb[1][0] + m_c[1][1] * gradNa_gradNb[1][1] + m_c[1][3] * gradNa_gradNb[1][2] +
                                        m_c[3][5] * gradNa_gradNb[2][0] + m_c[3][1] * gradNa_gradNb[2][1] + m_c[3][3] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+1][b*3+2] = elementStiffness[a*3+1][b*3+2] +
                                       (m_c[5][4] * gradNa_gradNb[0][0] + m_c[5][3] * gradNa_gradNb[0][1] + m_c[5][2] * gradNa_gradNb[0][2] +
                                        m_c[1][4] * gradNa_gradNb[1][0] + m_c[1][3] * gradNa_gradNb[1][1] + m_c[1][2] * gradNa_gradNb[1][2] +
                                        m_c[3][4] * gradNa_gradNb[2][0] + m_c[3][3] * gradNa_gradNb[2][1] + m_c[3][2] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+2][b*3+0] = elementStiffness[a*3+2][b*3+0] +
                                       (m_c[4][0] * gradNa_gradNb[0][0] + m_c[4][5] * gradNa_gradNb[0][1] + m_c[4][4] * gradNa_gradNb[0][2] +
                                        m_c[3][0] * gradNa_gradNb[1][0] + m_c[3][5] * gradNa_gradNb[1][1] + m_c[3][4] * gradNa_gradNb[1][2] +
                                        m_c[2][0] * gradNa_gradNb[2][0] + m_c[2][5] * gradNa_gradNb[2][1] + m_c[2][4] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+2][b*3+1] = elementStiffness[a*3+2][b*3+1] +
                                       (m_c[4][5] * gradNa_gradNb[0][0] + m_c[4][1] * gradNa_gradNb[0][1] + m_c[4][3] * gradNa_gradNb[0][2] +
                                        m_c[3][5] * gradNa_gradNb[1][0] + m_c[3][1] * gradNa_gradNb[1][1] + m_c[3][3] * gradNa_gradNb[1][2] +
                                        m_c[2][5] * gradNa_gradNb[2][0] + m_c[2][1] * gradNa_gradNb[2][1] + m_c[2][3] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+2][b*3+2] = elementStiffness[a*3+2][b*3+2] +
                                       (m_c[4][4] * gradNa_gradNb[0][0] + m_c[4][3] * gradNa_gradNb[0][1] + m_c[4][2] * gradNa_gradNb[0][2] +
                                        m_c[3][4] * gradNa_gradNb[1][0] + m_c[3][3] * gradNa_gradNb[1][1] + m_c[3][2] * gradNa_gradNb[1][2] +
                                        m_c[2][4] * gradNa_gradNb[2][0] + m_c[2][3] * gradNa_gradNb[2][1] + m_c[2][2] * gradNa_gradNb[2][2] ) * detJxW;
    }
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsFullyAnisotroipic::upperBTDB( BASIS_GRADIENT const & gradN,
                                                              real64 const & detJxW,
                                                              real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };


      elementStiffness[a*3+0][b*3+0] = elementStiffness[a*3+0][b*3+0] +
                                       ( m_c[0][0] * gradNa_gradNb[0][0] + m_c[0][5] * gradNa_gradNb[0][1] + m_c[0][4] * gradNa_gradNb[0][2] +
                                         m_c[5][0] * gradNa_gradNb[1][0] + m_c[5][5] * gradNa_gradNb[1][1] + m_c[5][4] * gradNa_gradNb[1][2] +
                                         m_c[4][0] * gradNa_gradNb[2][0] + m_c[4][5] * gradNa_gradNb[2][1] + m_c[4][4] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+0][b*3+1] = elementStiffness[a*3+0][b*3+1] +
                                       (m_c[0][5] * gradNa_gradNb[0][0] + m_c[0][1] * gradNa_gradNb[0][1] + m_c[0][3] * gradNa_gradNb[0][2] +
                                        m_c[5][5] * gradNa_gradNb[1][0] + m_c[5][1] * gradNa_gradNb[1][1] + m_c[5][3] * gradNa_gradNb[1][2] +
                                        m_c[4][5] * gradNa_gradNb[2][0] + m_c[4][1] * gradNa_gradNb[2][1] + m_c[4][3] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+0][b*3+2] = elementStiffness[a*3+0][b*3+2] +
                                       (m_c[0][4] * gradNa_gradNb[0][0] + m_c[0][3] * gradNa_gradNb[0][1] + m_c[0][2] * gradNa_gradNb[0][2] +
                                        m_c[5][4] * gradNa_gradNb[1][0] + m_c[5][3] * gradNa_gradNb[1][1] + m_c[5][2] * gradNa_gradNb[1][2] +
                                        m_c[4][4] * gradNa_gradNb[2][0] + m_c[4][3] * gradNa_gradNb[2][1] + m_c[4][2] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+1][b*3+0] = elementStiffness[a*3+1][b*3+0] +
                                       (m_c[5][0] * gradNa_gradNb[0][0] + m_c[5][5] * gradNa_gradNb[0][1] + m_c[5][4] * gradNa_gradNb[0][2] +
                                        m_c[1][0] * gradNa_gradNb[1][0] + m_c[1][5] * gradNa_gradNb[1][1] + m_c[1][4] * gradNa_gradNb[1][2] +
                                        m_c[3][0] * gradNa_gradNb[2][0] + m_c[3][5] * gradNa_gradNb[2][1] + m_c[3][4] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+1][b*3+1] = elementStiffness[a*3+1][b*3+1] +
                                       (m_c[5][5] * gradNa_gradNb[0][0] + m_c[5][1] * gradNa_gradNb[0][1] + m_c[5][3] * gradNa_gradNb[0][2] +
                                        m_c[1][5] * gradNa_gradNb[1][0] + m_c[1][1] * gradNa_gradNb[1][1] + m_c[1][3] * gradNa_gradNb[1][2] +
                                        m_c[3][5] * gradNa_gradNb[2][0] + m_c[3][1] * gradNa_gradNb[2][1] + m_c[3][3] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+1][b*3+2] = elementStiffness[a*3+1][b*3+2] +
                                       (m_c[5][4] * gradNa_gradNb[0][0] + m_c[5][3] * gradNa_gradNb[0][1] + m_c[5][2] * gradNa_gradNb[0][2] +
                                        m_c[1][4] * gradNa_gradNb[1][0] + m_c[1][3] * gradNa_gradNb[1][1] + m_c[1][2] * gradNa_gradNb[1][2] +
                                        m_c[3][4] * gradNa_gradNb[2][0] + m_c[3][3] * gradNa_gradNb[2][1] + m_c[3][2] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+2][b*3+0] = elementStiffness[a*3+2][b*3+0] +
                                       (m_c[4][0] * gradNa_gradNb[0][0] + m_c[4][5] * gradNa_gradNb[0][1] + m_c[4][4] * gradNa_gradNb[0][2] +
                                        m_c[3][0] * gradNa_gradNb[1][0] + m_c[3][5] * gradNa_gradNb[1][1] + m_c[3][4] * gradNa_gradNb[1][2] +
                                        m_c[2][0] * gradNa_gradNb[2][0] + m_c[2][5] * gradNa_gradNb[2][1] + m_c[2][4] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+2][b*3+1] = elementStiffness[a*3+2][b*3+1] +
                                       (m_c[4][5] * gradNa_gradNb[0][0] + m_c[4][1] * gradNa_gradNb[0][1] + m_c[4][3] * gradNa_gradNb[0][2] +
                                        m_c[3][5] * gradNa_gradNb[1][0] + m_c[3][1] * gradNa_gradNb[1][1] + m_c[3][3] * gradNa_gradNb[1][2] +
                                        m_c[2][5] * gradNa_gradNb[2][0] + m_c[2][1] * gradNa_gradNb[2][1] + m_c[2][3] * gradNa_gradNb[2][2] ) * detJxW;
      elementStiffness[a*3+2][b*3+2] = elementStiffness[a*3+2][b*3+2] +
                                       (m_c[4][4] * gradNa_gradNb[0][0] + m_c[4][3] * gradNa_gradNb[0][1] + m_c[4][2] * gradNa_gradNb[0][2] +
                                        m_c[3][4] * gradNa_gradNb[1][0] + m_c[3][3] * gradNa_gradNb[1][1] + m_c[3][2] * gradNa_gradNb[1][2] +
                                        m_c[2][4] * gradNa_gradNb[2][0] + m_c[2][3] * gradNa_gradNb[2][1] + m_c[2][2] * gradNa_gradNb[2][2] ) * detJxW;
    }
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsFullyAnisotroipic::diagBTDB( BASIS_GRADIENT const & gradN,
                                                             real64 const & detJxW,
                                                             real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS *3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    real64 const gradN_gradN[3][3] =
    { { gradN[a][0] * gradN[a][0], gradN[a][0] * gradN[a][1], gradN[a][0] * gradN[a][2] },
      { gradN[a][1] * gradN[a][0], gradN[a][1] * gradN[a][1], gradN[a][1] * gradN[a][2] },
      { gradN[a][2] * gradN[a][0], gradN[a][2] * gradN[a][1], gradN[a][2] * gradN[a][2] } };
    diagElementStiffness[a*3+0] = diagElementStiffness[a*3+0] +
                                  (m_c[0][0] * gradN_gradN[0][0] + m_c[0][5] * gradN_gradN[0][1] + m_c[0][4] * gradN_gradN[0][2] +
                                   m_c[5][0] * gradN_gradN[1][0] + m_c[5][5] * gradN_gradN[1][1] + m_c[5][4] * gradN_gradN[1][2] +
                                   m_c[4][0] * gradN_gradN[2][0] + m_c[4][5] * gradN_gradN[2][1] + m_c[4][4] * gradN_gradN[2][2] ) * detJxW;
    diagElementStiffness[a*3+1] = diagElementStiffness[a*3+1] +
                                  (m_c[5][5] * gradN_gradN[0][0] + m_c[5][1] * gradN_gradN[0][1] + m_c[5][3] * gradN_gradN[0][2] +
                                   m_c[1][5] * gradN_gradN[1][0] + m_c[1][1] * gradN_gradN[1][1] + m_c[1][3] * gradN_gradN[1][2] +
                                   m_c[3][5] * gradN_gradN[2][0] + m_c[3][1] * gradN_gradN[2][1] + m_c[3][3] * gradN_gradN[2][2] ) * detJxW;
    diagElementStiffness[a*3+2] = diagElementStiffness[a*3+2] +
                                  (m_c[4][4] * gradN_gradN[0][0] + m_c[4][3] * gradN_gradN[0][1] + m_c[4][2] * gradN_gradN[0][2] +
                                   m_c[3][4] * gradN_gradN[1][0] + m_c[3][3] * gradN_gradN[1][1] + m_c[3][2] * gradN_gradN[1][2] +
                                   m_c[2][4] * gradN_gradN[2][0] + m_c[2][3] * gradN_gradN[2][1] + m_c[2][2] * gradN_gradN[2][2] ) * detJxW;
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsFullyAnisotroipic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                                   real64 const & detJxW,
                                                                   real64 ( & diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };


      diagSumElementStiffness[a*3+0] = diagSumElementStiffness[a*3+0] +
                                       ( m_c[0][0] * gradNa_gradNb[0][0] + m_c[0][5] * gradNa_gradNb[0][1] + m_c[0][4] * gradNa_gradNb[0][2] +
                                         m_c[5][0] * gradNa_gradNb[1][0] + m_c[5][5] * gradNa_gradNb[1][1] + m_c[5][4] * gradNa_gradNb[1][2] +
                                         m_c[4][0] * gradNa_gradNb[2][0] + m_c[4][5] * gradNa_gradNb[2][1] + m_c[4][4] * gradNa_gradNb[2][2] +
                                         m_c[0][5] * gradNa_gradNb[0][0] + m_c[0][1] * gradNa_gradNb[0][1] + m_c[0][3] * gradNa_gradNb[0][2] +
                                         m_c[5][5] * gradNa_gradNb[1][0] + m_c[5][1] * gradNa_gradNb[1][1] + m_c[5][3] * gradNa_gradNb[1][2] +
                                         m_c[4][5] * gradNa_gradNb[2][0] + m_c[4][1] * gradNa_gradNb[2][1] + m_c[4][3] * gradNa_gradNb[2][2] +
                                         m_c[0][4] * gradNa_gradNb[0][0] + m_c[0][3] * gradNa_gradNb[0][1] + m_c[0][2] * gradNa_gradNb[0][2] +
                                         m_c[5][4] * gradNa_gradNb[1][0] + m_c[5][3] * gradNa_gradNb[1][1] + m_c[5][2] * gradNa_gradNb[1][2] +
                                         m_c[4][4] * gradNa_gradNb[2][0] + m_c[4][3] * gradNa_gradNb[2][1] + m_c[4][2] * gradNa_gradNb[2][2] ) * detJxW;
      diagSumElementStiffness[a*3+1] = diagSumElementStiffness[a*3+1] +
                                       (m_c[5][0] * gradNa_gradNb[0][0] + m_c[5][5] * gradNa_gradNb[0][1] + m_c[5][4] * gradNa_gradNb[0][2] +
                                        m_c[1][0] * gradNa_gradNb[1][0] + m_c[1][5] * gradNa_gradNb[1][1] + m_c[1][4] * gradNa_gradNb[1][2] +
                                        m_c[3][0] * gradNa_gradNb[2][0] + m_c[3][5] * gradNa_gradNb[2][1] + m_c[3][4] * gradNa_gradNb[2][2] +
                                        m_c[5][5] * gradNa_gradNb[0][0] + m_c[5][1] * gradNa_gradNb[0][1] + m_c[5][3] * gradNa_gradNb[0][2] +
                                        m_c[1][5] * gradNa_gradNb[1][0] + m_c[1][1] * gradNa_gradNb[1][1] + m_c[1][3] * gradNa_gradNb[1][2] +
                                        m_c[3][5] * gradNa_gradNb[2][0] + m_c[3][1] * gradNa_gradNb[2][1] + m_c[3][3] * gradNa_gradNb[2][2] +
                                        m_c[5][4] * gradNa_gradNb[0][0] + m_c[5][3] * gradNa_gradNb[0][1] + m_c[5][2] * gradNa_gradNb[0][2] +
                                        m_c[1][4] * gradNa_gradNb[1][0] + m_c[1][3] * gradNa_gradNb[1][1] + m_c[1][2] * gradNa_gradNb[1][2] +
                                        m_c[3][4] * gradNa_gradNb[2][0] + m_c[3][3] * gradNa_gradNb[2][1] + m_c[3][2] * gradNa_gradNb[2][2] ) * detJxW;
      diagSumElementStiffness[a*3+2] = diagSumElementStiffness[a*3+2] +
                                       (m_c[4][0] * gradNa_gradNb[0][0] + m_c[4][5] * gradNa_gradNb[0][1] + m_c[4][4] * gradNa_gradNb[0][2] +
                                        m_c[3][0] * gradNa_gradNb[1][0] + m_c[3][5] * gradNa_gradNb[1][1] + m_c[3][4] * gradNa_gradNb[1][2] +
                                        m_c[2][0] * gradNa_gradNb[2][0] + m_c[2][5] * gradNa_gradNb[2][1] + m_c[2][4] * gradNa_gradNb[2][2] +
                                        m_c[4][5] * gradNa_gradNb[0][0] + m_c[4][1] * gradNa_gradNb[0][1] + m_c[4][3] * gradNa_gradNb[0][2] +
                                        m_c[3][5] * gradNa_gradNb[1][0] + m_c[3][1] * gradNa_gradNb[1][1] + m_c[3][3] * gradNa_gradNb[1][2] +
                                        m_c[2][5] * gradNa_gradNb[2][0] + m_c[2][1] * gradNa_gradNb[2][1] + m_c[2][3] * gradNa_gradNb[2][2] +
                                        m_c[4][4] * gradNa_gradNb[0][0] + m_c[4][3] * gradNa_gradNb[0][1] + m_c[4][2] * gradNa_gradNb[0][2] +
                                        m_c[3][4] * gradNa_gradNb[1][0] + m_c[3][3] * gradNa_gradNb[1][1] + m_c[3][2] * gradNa_gradNb[1][2] +
                                        m_c[2][4] * gradNa_gradNb[2][0] + m_c[2][3] * gradNa_gradNb[2][1] + m_c[2][2] * gradNa_gradNb[2][2] ) * detJxW;
    }
  }
}

}
}


#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIAZTIONOPSANISOTROPIC_HPP_ */
