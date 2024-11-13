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
 * @file SolidModelDiscretizationOpsTransverseIsotropic.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSTRANSVERSEISOTROPIC_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSTRANSVERSEISOTROPIC_HPP_

#include "SolidModelDiscretizationOps.hpp"

namespace geos
{
namespace constitutive
{

/// Transverse isotropic implementation of discOps concept
struct SolidModelDiscretizationOpsTransverseIsotropic : public SolidModelDiscretizationOps
{
  /// @copydoc SolidModelDiscretizationOpsIsotropic::BTDB
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void BTDB( BASIS_GRADIENT const & gradN,
             real64 const & detJxW,
             real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  /// @copydoc SolidModelDiscretizationOpsIsotropic::upperBTDB
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void upperBTDB( BASIS_GRADIENT const & gradN,
                  real64 const & detJxW,
                  real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  /// @copydoc SolidModelDiscretizationOpsIsotropic::diagBTDB
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 const & detJxW,
                 real64 ( &diagElementStiffness )[NUM_SUPPORT_POINTS*3] );

  /// @copydoc SolidModelDiscretizationOpsIsotropic::diagRowSumBTDB
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 const & detJxW,
                       real64 ( &diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] );

  /// @copydoc SolidModelDiscretizationOpsIsotropic::scaleParams
  GEOS_HOST_DEVICE
  inline
  void scaleParams( real64 const scale )
  {
    m_c11 *= scale;
    m_c13 *= scale;
    m_c33 *= scale;
    m_c44 *= scale;
    m_c66 *= scale;
  }

  real64 m_c11; ///< (1,1) element of TI stiffness matrix
  real64 m_c13; ///< (1,3) element of TI stiffness matrix
  real64 m_c33; ///< (3,3) element of TI stiffness matrix
  real64 m_c44; ///< (4,4) element of TI stiffness matrix
  real64 m_c66; ///< (6,6) element of TI stiffness matrix
};

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsTransverseIsotropic::BTDB( BASIS_GRADIENT const & gradN,
                                                           real64 const & detJxW,
                                                           real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const c12 = (m_c11 - 2 * m_c66) * detJxW;
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c13 = this->m_c13 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;

  SolidModelDiscretizationOps::BTDB< NUM_SUPPORT_POINTS >( gradN,
                                                           elementStiffness,
                                                           [ c11,
                                                             c13,
                                                             c33,
                                                             c44,
                                                             c66,
                                                             c12 ] GEOS_HOST_DEVICE
                                                             ( int const a,
                                                             int const b,
                                                             real64 const (&gradNa_gradNb)[3][3],
                                                             real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    elementStiffness[a*3+0][b*3+0] = elementStiffness[a*3+0][b*3+0] + c11 * gradNa_gradNb[0][0] + c66 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+0][b*3+1] = elementStiffness[a*3+0][b*3+1] + c12 * gradNa_gradNb[0][1] + c66 * gradNa_gradNb[1][0];
    elementStiffness[a*3+0][b*3+2] = elementStiffness[a*3+0][b*3+2] + c13 * gradNa_gradNb[0][2] + c44 * gradNa_gradNb[2][0];
    elementStiffness[a*3+1][b*3+0] = elementStiffness[a*3+1][b*3+0] + c66 * gradNa_gradNb[0][1] + c12 * gradNa_gradNb[1][0];
    elementStiffness[a*3+1][b*3+1] = elementStiffness[a*3+1][b*3+1] + c66 * gradNa_gradNb[0][0] + c11 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+1][b*3+2] = elementStiffness[a*3+1][b*3+2] + c13 * gradNa_gradNb[1][2] + c44 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+0] = elementStiffness[a*3+2][b*3+0] + c44 * gradNa_gradNb[0][2] + c13 * gradNa_gradNb[2][0];
    elementStiffness[a*3+2][b*3+1] = elementStiffness[a*3+2][b*3+1] + c44 * gradNa_gradNb[1][2] + c13 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+2] = elementStiffness[a*3+2][b*3+2] + c44 * gradNa_gradNb[0][0] + c44 * gradNa_gradNb[1][1] + c33 * gradNa_gradNb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsTransverseIsotropic::upperBTDB( BASIS_GRADIENT const & gradN,
                                                                real64 const & detJxW,
                                                                real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const c12 = (m_c11 - 2 * m_c66) * detJxW;
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c13 = this->m_c13 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;

  SolidModelDiscretizationOps::upperBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                elementStiffness,
                                                                [ c11,
                                                                  c13,
                                                                  c33,
                                                                  c44,
                                                                  c66,
                                                                  c12 ] GEOS_HOST_DEVICE
                                                                  ( int const a,
                                                                  int const b,
                                                                  real64 const (&gradNa_gradNb)[3][3],
                                                                  real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    elementStiffness[a*3+0][b*3+0] = elementStiffness[a*3+0][b*3+0] + c11 * gradNa_gradNb[0][0] + c66 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+0][b*3+1] = elementStiffness[a*3+0][b*3+1] + c12 * gradNa_gradNb[0][1] + c66 * gradNa_gradNb[1][0];
    elementStiffness[a*3+0][b*3+2] = elementStiffness[a*3+0][b*3+2] + c13 * gradNa_gradNb[0][2] + c44 * gradNa_gradNb[2][0];
    elementStiffness[a*3+1][b*3+0] = elementStiffness[a*3+1][b*3+0] + c66 * gradNa_gradNb[0][1] + c12 * gradNa_gradNb[1][0];
    elementStiffness[a*3+1][b*3+1] = elementStiffness[a*3+1][b*3+1] + c66 * gradNa_gradNb[0][0] + c11 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+1][b*3+2] = elementStiffness[a*3+1][b*3+2] + c13 * gradNa_gradNb[1][2] + c44 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+0] = elementStiffness[a*3+2][b*3+0] + c44 * gradNa_gradNb[0][2] + c13 * gradNa_gradNb[2][0];
    elementStiffness[a*3+2][b*3+1] = elementStiffness[a*3+2][b*3+1] + c44 * gradNa_gradNb[1][2] + c13 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+2] = elementStiffness[a*3+2][b*3+2] + c44 * gradNa_gradNb[0][0] + c44 * gradNa_gradNb[1][1] + c33 * gradNa_gradNb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsTransverseIsotropic::diagBTDB( BASIS_GRADIENT const & gradN,
                                                               real64 const & detJxW,
                                                               real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS *3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;
  SolidModelDiscretizationOps::diagBTDB< NUM_SUPPORT_POINTS,
                                         3 >( gradN,
                                              diagElementStiffness,
                                              [ c11,
                                                c33,
                                                c44,
                                                c66 ] GEOS_HOST_DEVICE
                                                ( const int a,
                                                const real64 (& gradN_gradN)[3],
                                                real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagElementStiffness[a*3+0] = diagElementStiffness[a*3+0] + c11 * gradN_gradN[0] + c66 * gradN_gradN[1] + c44 * gradN_gradN[2];
    diagElementStiffness[a*3+1] = diagElementStiffness[a*3+1] + c66 * gradN_gradN[0] + c11 * gradN_gradN[1] + c44 * gradN_gradN[2];
    diagElementStiffness[a*3+2] = diagElementStiffness[a*3+2] + c44 * gradN_gradN[0] + c44 * gradN_gradN[1] + c33 * gradN_gradN[2];
  } );
}

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsTransverseIsotropic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                                     real64 const & detJxW,
                                                                     real64 ( & diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] )
{
  real64 const c12 = (m_c11 - 2 * m_c66) * detJxW;
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c13 = this->m_c13 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;

  SolidModelDiscretizationOps::diagRowSumBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                     diagSumElementStiffness,
                                                                     [ c11,
                                                                       c13,
                                                                       c33,
                                                                       c44,
                                                                       c66,
                                                                       c12 ] GEOS_HOST_DEVICE
                                                                       ( int const a,
                                                                       real64 const (&gradNa_gradNb)[3][3],
                                                                       real64 (& diagSumElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagSumElementStiffness[a*3+0] = diagSumElementStiffness[a*3+0] +
                                     c11 * gradNa_gradNb[0][0] +
                                     c66 * gradNa_gradNb[1][1] +
                                     c44 * gradNa_gradNb[2][2] +
                                     c12 * gradNa_gradNb[0][1] +
                                     c66 * gradNa_gradNb[1][0] +
                                     c13 * gradNa_gradNb[0][2] +
                                     c44 * gradNa_gradNb[2][0];
    diagSumElementStiffness[a*3+1] = diagSumElementStiffness[a*3+1] +
                                     c66 * gradNa_gradNb[0][1] +
                                     c12 * gradNa_gradNb[1][0] +
                                     c66 * gradNa_gradNb[0][0] +
                                     c11 * gradNa_gradNb[1][1] +
                                     c44 * gradNa_gradNb[2][2] +
                                     c13 * gradNa_gradNb[1][2] +
                                     c44 * gradNa_gradNb[2][1];
    diagSumElementStiffness[a*3+2] = diagSumElementStiffness[a*3+2] +
                                     c44 * gradNa_gradNb[0][2] +
                                     c13 * gradNa_gradNb[2][0] +
                                     c44 * gradNa_gradNb[1][2] +
                                     c13 * gradNa_gradNb[2][1] +
                                     c44 * gradNa_gradNb[0][0] +
                                     c44 * gradNa_gradNb[1][1] +
                                     c33 * gradNa_gradNb[2][2];
  } );
}
#if __GNUC__
#pragma GCC diagnostic pop
#endif

}
}


#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSTRANSVERSEISOTROPIC_HPP_ */
