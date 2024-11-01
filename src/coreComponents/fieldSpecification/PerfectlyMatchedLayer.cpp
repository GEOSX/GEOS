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

/*
 * PerfectlyMatchedLayer.cpp
 *
 */

#include "PerfectlyMatchedLayer.hpp"

namespace geos
{
using namespace dataRepository;

PerfectlyMatchedLayer::PerfectlyMatchedLayer( string const & name, Group * const parent ):
  FieldSpecificationBase( name, parent )
{
  registerWrapper( viewKeyStruct::xMinString(), &m_xMin ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { -LvArray::NumericLimits< real32 >::max,
                            -LvArray::NumericLimits< real32 >::max,
                            -LvArray::NumericLimits< real32 >::max } ).
    setDescription( "Minimum (x,y,z) coordinates of the inner PML boundaries" );

  registerWrapper( viewKeyStruct::xMaxString(), &m_xMax ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { LvArray::NumericLimits< real32 >::max,
                            LvArray::NumericLimits< real32 >::max,
                            LvArray::NumericLimits< real32 >::max } ).
    setDescription( "Maximum (x,y,z) coordinates of the inner PML boundaries" );

  registerWrapper( viewKeyStruct::reflectivityString(), &m_reflectivity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.001 ).
    setDescription( "Desired reflectivity of the PML region, used to compute the damping profile" );

  registerWrapper( viewKeyStruct::thicknessMinXYZString(), &m_thicknessMinXYZ ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {-1.0, -1.0, -1.0} ).
    setDescription( "Thickness of the PML region, at left, front, and top sides, used to compute the damping profile" );

  registerWrapper( viewKeyStruct::thicknessMaxXYZString(), &m_thicknessMaxXYZ ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {-1.0, -1.0, -1.0} ).
    setDescription( "Thickness of the PML region, at right, back, and bottom sides, used to compute the damping profile" );

  registerWrapper( viewKeyStruct::waveSpeedMinXYZString(), &m_waveSpeedMinXYZ ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {-1.0, -1.0, -1.0} ).
    setDescription( "Wave speed in the PML, at left, front, and top sides, used to compute the damping profile" );

  registerWrapper( viewKeyStruct::waveSpeedMaxXYZString(), &m_waveSpeedMaxXYZ ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {-1.0, -1.0, -1.0} ).
    setDescription( "Wave speed in the PML, at right, back, and bottom sides, used to compute the damping profile" );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::initialConditionString() ).
    setInputFlag( InputFlags::FALSE );
  initialCondition( false ); // to make sure this is not called by applyInitialConditions

}

void PerfectlyMatchedLayer::postInputInitialization()
{
  GEOS_THROW_IF( (m_xMax[0]<m_xMin[0] || m_xMax[1]<m_xMin[1] || m_xMax[2]<m_xMin[2]),
                 getCatalogName() << " " << getDataContext() << " "
                                  << viewKeyStruct::xMinString()
                                  << " must be smaller than "
                                  << viewKeyStruct::xMaxString(),
                 InputError );

  GEOS_THROW_IF( (m_reflectivity<=0 || m_reflectivity>1),
                 getCatalogName() << " " << getDataContext() << " "
                                  << viewKeyStruct::reflectivityString()
                                  << " must satisfy 0 < reflectivity <= 1",
                 InputError );

  GEOS_LOG_RANK_0_IF( (m_xMin[0]<smallestXMin || m_xMin[1]<smallestXMin || m_xMin[2]<smallestXMin),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::xMinString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_xMax[0]>largestXMax || m_xMax[1]>largestXMax || m_xMax[2]>largestXMax),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::xMaxString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_thicknessMinXYZ[0]<0 || m_thicknessMinXYZ[1]<0 || m_thicknessMinXYZ[2]<0),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::thicknessMinXYZString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_thicknessMaxXYZ[0]<0 || m_thicknessMaxXYZ[1]<0 || m_thicknessMaxXYZ[2]<0),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::thicknessMaxXYZString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_waveSpeedMinXYZ[0]<0 || m_waveSpeedMinXYZ[1]<0 || m_waveSpeedMinXYZ[2]<0),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::waveSpeedMinXYZString()
                                       << " will be computed internally" );

  GEOS_LOG_RANK_0_IF( (m_waveSpeedMaxXYZ[0]<0 || m_waveSpeedMaxXYZ[1]<0 || m_waveSpeedMaxXYZ[2]<0),
                      getCatalogName() << " " << getDataContext() << " "
                                       << viewKeyStruct::waveSpeedMaxXYZString()
                                       << " will be computed internally" );
}


REGISTER_CATALOG_ENTRY( FieldSpecificationBase, PerfectlyMatchedLayer, string const &, Group * const )

} /* namespace geos */
