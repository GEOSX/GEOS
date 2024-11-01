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
 * @file Cylinder.cpp
 * @brief Generate a cylindrical geometry.
 * @param Center point of one (upper or lower) face of the cylinder
 * @param Center point of the other face of the cylinder
 * @param Radius of the cylinder
 */

#include "Cylinder.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
using namespace dataRepository;

Cylinder::Cylinder( const string & name, Group * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_point1{ 0.0, 0.0, 0.0 },
  m_point2{ 0.0, 0.0, 0.0 },
  m_radius{ 0.0 }
{
  registerWrapper( viewKeyStruct::point1String(), &m_point1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Center point of the first face of the cylinder" );

  registerWrapper( viewKeyStruct::point2String(), &m_point2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Center point of the second face of the cylinder" );

  registerWrapper( viewKeyStruct::radiusString(), &m_radius ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Outer radius of the cylinder" );

  registerWrapper( viewKeyStruct::innerRadiusString(), &m_innerRadius ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Inner radius of the annulus" );

}

Cylinder::~Cylinder()
{}


bool Cylinder::isCoordInObject( real64 const ( &targetPt ) [3] ) const
{
  bool rval = false;

  // Assuming the cylinder is defined by (pt1,pt2,innerRadius,outerRadius),
  // we check that the target point is inside using the following formulas
  //
  // 1) Check that the target point lies between the planes of the two circular facets of the cylinder
  // ( \vec{targetPt} - \vec{pt1} ) \cdot ( \vec{pt2} - \vec{pt1} ) > 0
  // ( \vec{targetPt} - \vec{pt2} ) \cdot ( \vec{pt2} - \vec{pt1} ) < 0
  //
  // 2) Check that the target point lies inside the curved surface of the cylinder
  // \frac{| ( \vec{targetPt} - \vec{pt1} ) x ( \vec{pt2} - \vec{pt1} ) |}{| \vec{pt2} - \vec{pt1} |} \geq innerRadius
  // \frac{| ( \vec{targetPt} - \vec{pt1} ) x ( \vec{pt2} - \vec{pt1} ) |}{| \vec{pt2} - \vec{pt1} |} < outerRadius

  real64 pt2Pt1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_point2 );
  LvArray::tensorOps::subtract< 3 >( pt2Pt1, m_point1 );

  real64 targetPtPt1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( targetPt );
  LvArray::tensorOps::subtract< 3 >( targetPtPt1, m_point1 );

  real64 targetPtPt2[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( targetPt );
  LvArray::tensorOps::subtract< 3 >( targetPtPt2, m_point2 );

  real64 const dotProd_pt2Pt1_targetPtPt1 = LvArray::tensorOps::AiBi< 3 >( pt2Pt1, targetPtPt1 );
  real64 const dotProd_pt2Pt1_targetPtPt2 = LvArray::tensorOps::AiBi< 3 >( pt2Pt1, targetPtPt2 );

  real64 crossProd[3]{};
  LvArray::tensorOps::crossProduct( crossProd, targetPtPt1, pt2Pt1 );
  real64 const radius = LvArray::tensorOps::l2Norm< 3 >( crossProd ) / LvArray::tensorOps::l2Norm< 3 >( pt2Pt1 );

  if( radius < m_radius &&
      radius >= m_innerRadius &&
      dotProd_pt2Pt1_targetPtPt1 > 0 &&
      dotProd_pt2Pt1_targetPtPt2 < 0 )
  {
    rval = true;
  }

  return rval;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Cylinder, string const &, Group * const )

} /* namespace geos */
