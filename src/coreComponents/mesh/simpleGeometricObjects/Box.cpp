/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file Box.cpp
 * @brief Generate a box geometry.
 * @param Maximum (x,y,z) coordinates of the box
 * @param Minimum (x,y,z) coordinates of the box
 * @param The strike angle of the box
 */

#include "Box.hpp"
#include "LvArray/src/genericTensorOps.hpp"

namespace geos
{
using namespace dataRepository;

Box::Box( const string & name, Group * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_min{ 0.0, 0.0, 0.0 },
  m_max{ 0.0, 0.0, 0.0 },
  m_strikeAngle{ 0.0 },
  m_boxCenter{ 0.0, 0.0, 0.0 },
  m_cosStrike{ 0.0 },
  m_sinStrike{ 0.0 }
{
  registerWrapper( viewKeyStruct::xMinString(), &m_min ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Minimum (x,y,z) coordinates of the box" );

  registerWrapper( viewKeyStruct::xMaxString(), &m_max ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum (x,y,z) coordinates of the box" );

  registerWrapper( viewKeyStruct::strikeAngleString(), &m_strikeAngle ).
    setApplyDefaultValue( -90.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The strike angle of the box" );

  registerWrapper( viewKeyStruct::boxCenterString(), &m_boxCenter );
  registerWrapper( viewKeyStruct::cosStrikeString(), &m_cosStrike );
  registerWrapper( viewKeyStruct::sinStrikeString(), &m_sinStrike );
}

Box::~Box()
{}



void Box::postInputInitialization()
{
  LvArray::tensorOps::copy< 3 >( m_boxCenter, m_min );
  LvArray::tensorOps::add< 3 >( m_boxCenter, m_max );
  LvArray::tensorOps::scale< 3 >( m_boxCenter, 0.5 );

  // reordering min and max to fit former interface but so that user can input any of the four diagonals
  for( int i = 0; i < m_max.SIZE; ++i )
  {
    if( m_max[i]<m_min[i] )
    {
      std::swap( m_max[i], m_min[i] );
      GEOS_LOG_RANK_0( GEOS_FMT( "Reordering box definition for {} component as {} < {} ", i, m_max[i], m_min[i] ) );
    }
  }

  m_strikeAngle += 90; // Counterclockwise from x-axis
  if( std::fabs( m_strikeAngle ) > 1e-20 )
  {
    GEOS_ERROR_IF( (m_max[0]-m_min[0]) < (m_max[1]-m_min[1]),
                   getDataContext() << ": When a strike angle is specified, the box is supposed to" <<
                   " represent a plane normal to the y direction. This box seems to be too thick." );

    m_cosStrike = std::cos( m_strikeAngle / 180 *M_PI );
    m_sinStrike = std::sin( m_strikeAngle / 180 *M_PI );
  }
}

bool Box::isCoordInObject( real64 const ( &coord ) [3] ) const
{
  real64 coord0[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( coord );
  if( std::fabs( m_strikeAngle ) >= 1e-20 )
  {
    real64 coordR[3];
    LvArray::tensorOps::subtract< 3 >( coord0, m_boxCenter );
    coordR[0] =  coord0[0] * m_cosStrike + coord0[1] * m_sinStrike;
    coordR[1] = -coord0[0] * m_sinStrike + coord0[1] * m_cosStrike;
    coordR[2] = coord0[2];
    LvArray::tensorOps::add< 3 >( coordR, m_boxCenter );
    LvArray::tensorOps::copy< 3 >( coord0, coordR );
  }
  for( int i = 0; i < 3; ++i )
  {
    if( coord0[i] < m_min[i] || coord0[i] > m_max[i] )
    {
      return false;
    }
  }
  return true;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Box, string const &, Group * const )

} /* namespace geos */
