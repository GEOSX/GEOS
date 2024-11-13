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
 * @file Rectangle.cpp
 */

#include "Rectangle.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
using namespace dataRepository;

Rectangle::Rectangle( const string & name, Group * const parent ):
  PlanarGeometricObject( name, parent ),
  m_origin{ 0.0, 0.0, 0.0 },
  m_tolerance()
{
  registerWrapper( viewKeyStruct::originString(), &m_origin ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Origin point (x,y,z) of the plane (basically, any point on the plane)" );

  registerWrapper( viewKeyStruct::dimensionsString(), &m_dimensions ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Length and width of the bounded plane" );

  registerWrapper( viewKeyStruct::toleranceString(), &m_tolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 1e-5 ).
    setDescription( "Tolerance to determine if a point sits on the plane or not. "
                    "It is relative to the maximum dimension of the plane." );


  m_points.resize( 4, 3 );
}

//constructor given two points, used for 2.5D problems
Rectangle::Rectangle( const real64 oldX, const real64 oldY,
                      const real64 newX, const real64 newY,
                      const string & name, Group * const parent ):
  PlanarGeometricObject( name, parent ),
  m_origin{ 0.0, 0.0, 0.0 },
  m_tolerance( 1e-5 )
{
  m_origin = { (oldX + newX)/2.0, (oldY + newY)/2.0, 0.0};
  m_normal = { -(newY-oldY), newX-oldX, 0.0 };
  m_lengthVector = { newX-oldX, newY-oldY, 0.0 };
  m_widthVector = { 0.0, 0.0, 1.0 };
  real64 norm = std::sqrt( pow( newX-oldX, 2 )+pow( newY-oldY, 2 ));
  m_dimensions.resize( 2 );
  m_dimensions[0] = norm+1e-4; //small tolerance to ensure that both ends are contained in the plane - TODO: try to use m_tolerance
  m_dimensions[1] = 5; //TODO: this is arbitrary, it only needs to be larger than the z thickness in the 2.5D model

  m_points.resize( 4, 3 );
  //parent->registerGroup< Rectangle >( name, std::move( this ) );
  this->postInputInitialization();
}

Rectangle::~Rectangle()
{}

void Rectangle::postInputInitialization()
{
  // Make sure that you have an orthonormal basis.
  LvArray::tensorOps::normalize< 3 >( m_normal );
  LvArray::tensorOps::normalize< 3 >( m_lengthVector );
  LvArray::tensorOps::normalize< 3 >( m_widthVector );

  m_tolerance = m_tolerance * std::min( m_dimensions[0], m_dimensions[1] );

  //Check if they are all orthogonal
  real64 vector[ 3 ];
  LvArray::tensorOps::crossProduct( vector, m_lengthVector, m_widthVector );

  GEOS_ERROR_IF( std::fabs( std::fabs( LvArray::tensorOps::AiBi< 3 >( m_normal, vector )) - 1 ) > orthoNormalBaseTolerance
                 || std::fabs( LvArray::tensorOps::AiBi< 3 >( m_widthVector, m_lengthVector )) > orthoNormalBaseTolerance,
                 "Error: the 3 vectors provided in the Rectangle do not form an orthonormal basis!" );

  GEOS_ERROR_IF( m_dimensions.size() != 2, "Error: Need to provide both length and width!" );

  findRectangleLimits();
}

void Rectangle::findRectangleLimits()
{
  real64 lengthVec[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_lengthVector );
  real64 widthVec[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_widthVector );

  LvArray::tensorOps::scale< 3 >( lengthVec, 0.5 * m_dimensions[0] );
  LvArray::tensorOps::scale< 3 >( widthVec, 0.5 * m_dimensions[1] );

  for( int i = 0; i < 4; i++ )
  {
    LvArray::tensorOps::copy< 3 >( m_points[i], m_origin );
  }

  LvArray::tensorOps::subtract< 3 >( m_points[0], lengthVec );
  LvArray::tensorOps::subtract< 3 >( m_points[0], widthVec );

  LvArray::tensorOps::add< 3 >( m_points[1], lengthVec );
  LvArray::tensorOps::subtract< 3 >( m_points[1], widthVec );

  LvArray::tensorOps::add< 3 >( m_points[2], lengthVec );
  LvArray::tensorOps::add< 3 >( m_points[2], widthVec );

  LvArray::tensorOps::subtract< 3 >( m_points[3], lengthVec );
  LvArray::tensorOps::add< 3 >( m_points[3], widthVec );

  if( getLogLevel() > 1 )
  {
    GEOS_LOG_RANK_0( "Point A: " << m_points[0] );
    GEOS_LOG_RANK_0( "Point B: " << m_points[1] );
    GEOS_LOG_RANK_0( "Point C: " << m_points[2] );
    GEOS_LOG_RANK_0( "Point D: " << m_points[3] );
  }
}

bool Rectangle::isCoordInObject( real64 const ( &coord ) [3] ) const
{
  bool isInside = true;

  real64 dummy[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( coord );
  LvArray::tensorOps::subtract< 3 >( dummy, m_origin );

  // 1. Check if point is on the plane
  if( std::abs( LvArray::tensorOps::AiBi< 3 >( dummy, m_normal ) ) < m_tolerance )
  {
    real64 vec[ 3 ]   = LVARRAY_TENSOROPS_INIT_LOCAL_3( coord );
    real64 abVec[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_points[1] );
    real64 adVec[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_points[3] );

    LvArray::tensorOps::subtract< 3 >( vec, m_points[0] );
    LvArray::tensorOps::subtract< 3 >( abVec, m_points[0] );
    LvArray::tensorOps::subtract< 3 >( adVec, m_points[0] );

    real64 const abDotProd = LvArray::tensorOps::AiBi< 3 >( vec, abVec );
    real64 const adDotProd = LvArray::tensorOps::AiBi< 3 >( vec, adVec );

    // 2. Check if it is inside the rectangle
    if( abDotProd < 0 || abDotProd > LvArray::tensorOps::l2NormSquared< 3 >( abVec ) )
    {
      isInside = false;
    }

    if( adDotProd < 0 || adDotProd > LvArray::tensorOps::l2NormSquared< 3 >( adVec ) )
    {
      isInside = false;
    }

  }
  else
  {
    isInside = false;
  }

  return isInside;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Rectangle, string const &, Group * const )

} /* namespace geos */
