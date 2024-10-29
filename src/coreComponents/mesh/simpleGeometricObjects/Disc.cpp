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
 * @file Disc.cpp
 */

#include "Disc.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
using namespace dataRepository;

Disc::Disc( const string & name, Group * const parent ):
  PlanarGeometricObject( name, parent ),
  m_center{ 0.0, 0.0, 0.0 },
  m_radius( 1.0 ),
  m_tolerance()
{
  registerWrapper( viewKeyStruct::centerString(), &m_center ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "(x,y,z) coordinates of the center of the disc" );

  registerWrapper( viewKeyStruct::radiusString(), &m_radius ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Radius of the disc." );

  registerWrapper( viewKeyStruct::toleranceString(), &m_tolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 1e-5 ).
    setDescription( "Tolerance to determine if a point sits on the disc or not. "
                    "It is relative to the maximum dimension of the disc." );

}

Disc::~Disc()
{}

void Disc::postInputInitialization()
{
  // Make sure that you have an orthonormal basis.
  LvArray::tensorOps::normalize< 3 >( m_normal );
  m_tolerance = m_tolerance * m_radius;

}

bool Disc::isCoordInObject( real64 const ( &coord ) [3] ) const
{
  bool isInside = true;

  //gets the vector from coord to the center
  real64 dummy[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( coord );
  LvArray::tensorOps::subtract< 3 >( dummy, m_center );

  // 1. Check if point is on the plane of the disc
  if( std::abs( LvArray::tensorOps::AiBi< 3 >( dummy, m_normal ) ) > m_tolerance )
  {
    isInside = false;
  }

  // 2. Check if it is inside the disc
  if( LvArray::tensorOps::l2NormSquared< 3 >( dummy ) > m_radius*m_radius )
  {
    isInside = false;
  }

  return isInside;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Disc, string const &, Group * const )

} /* namespace geos */
