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

/**
 * @file AquiferBoundaryCondition.cpp
 */

#include "AquiferBoundaryCondition.hpp"

#include "mesh/DomainPartition.hpp"

namespace geos
{

using namespace dataRepository;

AquiferBoundaryCondition::AquiferBoundaryCondition( string const & name, Group * parent )
  : FieldSpecificationBase( name, parent ),
  m_waterPhaseIndex( -1 ),
  m_cumulativeFlux( 0.0 )
{
  registerWrapper( viewKeyStruct::aquiferPorosityString(), &m_porosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer porosity" );

  registerWrapper( viewKeyStruct::aquiferPermeabilityString(), &m_permeability ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer permeability [m^2]" );

  registerWrapper( viewKeyStruct::aquiferInitialPressureString(), &m_initialPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer initial pressure [Pa]" );

  registerWrapper( viewKeyStruct::aquiferWaterViscosityString(), &m_viscosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer water viscosity [Pa.s]" );

  registerWrapper( viewKeyStruct::aquiferWaterDensityString(), &m_density ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer water density [kg.m^-3]" );

  registerWrapper( viewKeyStruct::aquiferTotalCompressibilityString(), &m_totalCompressibility ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer total compressibility (rock and fluid) [Pa^-1]" );

  registerWrapper( viewKeyStruct::allowAllPhasesIntoAquiferString(), &m_allowAllPhasesIntoAquifer ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to allow all phases to flow into the aquifer. \n"
                    "This flag only matters for the configuration in which flow is from reservoir to aquifer. \n"
                    "    - If the flag is equal to 1, then all phases, including non-aqueous phases, are allowed to flow into the aquifer. \n "
                    "    - If the flag is equal to 0, then only the water phase is allowed to flow into the aquifer. \n"
                    "If you are in a configuration in which flow is from reservoir to aquifer "
                    "and you expect non-aqueous phases to saturate the reservoir cells next to the aquifer, set this flag to 1. \n"
                    "This keyword is ignored for single-phase flow simulations" );

  registerWrapper( viewKeyStruct::aquiferWaterPhaseComponentFractionString(), &m_phaseComponentFraction ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Aquifer water phase component fraction. This keyword is ignored for single-phase flow simulations." );

  registerWrapper( viewKeyStruct::aquiferWaterPhaseComponentNamesString(), &m_phaseComponentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Aquifer water phase component names. This keyword is ignored for single-phase flow simulations." );

  registerWrapper( viewKeyStruct::aquiferElevationString(), &m_elevation ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer elevation (positive going upward) [m]" );

  registerWrapper( viewKeyStruct::aquiferThicknessString(), &m_thickness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer thickness [m]" );

  registerWrapper( viewKeyStruct::aquiferInnerRadiusString(), &m_innerRadius ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Aquifer inner radius [m]" );

  registerWrapper( viewKeyStruct::aquiferAngleString(), &m_angle ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Angle subtended by the aquifer boundary from the center of the reservoir [degress]" );

  registerWrapper( viewKeyStruct::pressureInfluenceFunctionNameString(), &m_pressureInfluenceFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the table describing the pressure influence function\n. "
                    "If not provided, we use a default pressure influence function" );

  registerWrapper( viewKeyStruct::cumulativeFluxString(), &m_cumulativeFlux ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::objectPathString() ).
    setInputFlag( InputFlags::FALSE );
  setObjectPath( "faceManager" );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::FALSE );

}

void AquiferBoundaryCondition::postInputInitialization()
{
  GEOS_THROW_IF_LE_MSG( m_permeability, 0.0,
                        getCatalogName() << " " << getDataContext() <<
                        ": the aquifer permeability cannot be equal to zero or negative",
                        InputError );

  if( m_pressureInfluenceFunctionName.empty() )
  {
    setupDefaultPressureInfluenceFunction();
  }
  else
  {
    FunctionManager const & functionManager = FunctionManager::getInstance();
    GEOS_THROW_IF( !functionManager.hasGroup( m_pressureInfluenceFunctionName ),
                   getCatalogName() << " " << getDataContext() <<
                   ": the pressure influence table " << m_pressureInfluenceFunctionName << " could not be found",
                   InputError );

    TableFunction const & pressureInfluenceFunction = functionManager.getGroup< TableFunction >( m_pressureInfluenceFunctionName );
    GEOS_THROW_IF( pressureInfluenceFunction.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                   getCatalogName() << " " << getDataContext() <<
                   ": The interpolation method for the pressure influence function table " <<
                   pressureInfluenceFunction.getDataContext() <<
                   " should be TableFunction::InterpolationType::Linear",
                   InputError );
  }

  computeTimeConstant();
  computeInfluxConstant();

  GEOS_THROW_IF_LE_MSG( m_timeConstant, 0.0,
                        getCatalogName() << " " << getDataContext() <<
                        ": the aquifer time constant is equal to zero or negative, the simulation cannot procede",
                        InputError );

  GEOS_THROW_IF_LE_MSG( m_influxConstant, 0.0,
                        getCatalogName() << " " << getDataContext() <<
                        ": the aquifer influx constant is equal to zero or negative, the simulation cannot procede",
                        InputError );

  GEOS_THROW_IF_NE_MSG( m_phaseComponentFraction.size(), m_phaseComponentNames.size(),
                        getCatalogName() << " " << getDataContext() <<
                        ": the sizes of " << viewKeyStruct::aquiferWaterPhaseComponentFractionString() <<
                        " and " << viewKeyStruct::aquiferWaterPhaseComponentNamesString() << " are inconsistent",
                        InputError );

}

void AquiferBoundaryCondition::setupDefaultPressureInfluenceFunction()
{
  // default table; see Eclipse or Intersect documentation

  array1d< array1d< real64 > > dimensionlessTime;
  dimensionlessTime.resize( 1 );
  dimensionlessTime[0].resize( 42 );
  dimensionlessTime[0][0] = 0.01;
  dimensionlessTime[0][1] = 0.05;
  dimensionlessTime[0][2] = 0.1;
  dimensionlessTime[0][3] = 0.15;
  dimensionlessTime[0][4] = 0.2;
  dimensionlessTime[0][5] = 0.25;
  dimensionlessTime[0][6] = 0.3;
  dimensionlessTime[0][7] = 0.4;
  dimensionlessTime[0][8] = 0.5;
  dimensionlessTime[0][9] = 0.6;
  dimensionlessTime[0][10] = 0.7;
  dimensionlessTime[0][11] = 0.8;
  dimensionlessTime[0][12] = 0.9;
  dimensionlessTime[0][13] = 1.0;
  dimensionlessTime[0][14] = 1.5;
  dimensionlessTime[0][15] = 2.0;
  dimensionlessTime[0][16] = 2.5;
  dimensionlessTime[0][17] = 3.0;
  dimensionlessTime[0][18] = 4.0;
  dimensionlessTime[0][19] = 5.0;
  dimensionlessTime[0][20] = 6.0;
  dimensionlessTime[0][21] = 7.0;
  dimensionlessTime[0][22] = 8.0;
  dimensionlessTime[0][23] = 9.0;
  dimensionlessTime[0][24] = 10.0;
  dimensionlessTime[0][25] = 15.0;
  dimensionlessTime[0][26] = 20.0;
  dimensionlessTime[0][27] = 25.0;
  dimensionlessTime[0][28] = 30.0;
  dimensionlessTime[0][29] = 40.0;
  dimensionlessTime[0][30] = 50.0;
  dimensionlessTime[0][31] = 60.0;
  dimensionlessTime[0][32] = 70.0;
  dimensionlessTime[0][33] = 80.0;
  dimensionlessTime[0][34] = 90.0;
  dimensionlessTime[0][35] = 100.0;
  dimensionlessTime[0][36] = 200.0;
  dimensionlessTime[0][37] = 800.0;
  dimensionlessTime[0][38] = 1600.0;
  dimensionlessTime[0][39] = 3200.0;
  dimensionlessTime[0][40] = 6400.0;
  dimensionlessTime[0][41] = 12800.0;

  array1d< real64 > pressureInfluence;
  pressureInfluence.resize( 42 );
  pressureInfluence[0] = 0.112;
  pressureInfluence[1] = 0.229;
  pressureInfluence[2] = 0.315;
  pressureInfluence[3] = 0.376;
  pressureInfluence[4] = 0.424;
  pressureInfluence[5] = 0.469;
  pressureInfluence[6] = 0.503;
  pressureInfluence[7] = 0.564;
  pressureInfluence[8] = 0.616;
  pressureInfluence[9] = 0.659;
  pressureInfluence[10] = 0.702;
  pressureInfluence[11] = 0.735;
  pressureInfluence[12] = 0.772;
  pressureInfluence[13] = 0.802;
  pressureInfluence[14] = 0.927;
  pressureInfluence[15] = 1.02;
  pressureInfluence[16] = 1.101;
  pressureInfluence[17] = 1.169;
  pressureInfluence[18] = 1.275;
  pressureInfluence[19] = 1.362;
  pressureInfluence[20] = 1.436;
  pressureInfluence[21] = 1.5;
  pressureInfluence[22] = 1.556;
  pressureInfluence[23] = 1.604;
  pressureInfluence[24] = 1.651;
  pressureInfluence[25] = 1.829;
  pressureInfluence[26] = 1.96;
  pressureInfluence[27] = 2.067;
  pressureInfluence[28] = 2.147;
  pressureInfluence[29] = 2.282;
  pressureInfluence[30] = 2.388;
  pressureInfluence[31] = 2.476;
  pressureInfluence[32] = 2.55;
  pressureInfluence[33] = 2.615;
  pressureInfluence[34] = 2.672;
  pressureInfluence[35] = 2.723;
  pressureInfluence[36] = 3.0537;
  pressureInfluence[37] = 3.7468;
  pressureInfluence[38] = 4.0934;
  pressureInfluence[39] = 4.44;
  pressureInfluence[40] = 4.7866;
  pressureInfluence[41] = 5.1331;

  FunctionManager & functionManager = FunctionManager::getInstance();
  m_pressureInfluenceFunctionName = getName() + "_pressureInfluence_table";
  TableFunction * const pressureInfluenceTable =
    dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), m_pressureInfluenceFunctionName ) );
  pressureInfluenceTable->setTableCoordinates( dimensionlessTime, { units::Dimensionless } );
  pressureInfluenceTable->setTableValues( pressureInfluence, units::Dimensionless );
  pressureInfluenceTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );

}

void AquiferBoundaryCondition::setGravityVector( R1Tensor const & gravityVector )
{
  GEOS_LOG_RANK_0_IF( ( !isZero( gravityVector[0] ) || !isZero( gravityVector[1] ) ),
                      catalogName() << " " << getDataContext() <<
                      "The gravity vector specified in this simulation (" << gravityVector[0] << " " << gravityVector[1] << " " << gravityVector[2] <<
                      ") is not aligned with the z-axis. \n" <<
                      "But, the pressure difference between reservoir and aquifer uses " << viewKeyStruct::aquiferElevationString() <<
                      " and assumes that the gravity vector is aligned with the z-axis. \n" <<
                      "As a result, aquifer calculations may be inacurrate." );

  m_gravityVector = gravityVector;
}

void AquiferBoundaryCondition::computeTimeConstant()
{
  // equation 5.3 of the Eclipse TD
  m_timeConstant = m_viscosity * m_porosity * m_totalCompressibility * m_innerRadius * m_innerRadius;
  m_timeConstant /= m_permeability;
}

void AquiferBoundaryCondition::computeInfluxConstant()
{
  // equation 5.4 of the Eclipse TD, including the constant 6.283 of the Carter-Tracy model
  m_influxConstant = 6.283 * m_thickness * ( m_angle / 360.0 ) * m_porosity * m_totalCompressibility * m_innerRadius * m_innerRadius;
}

AquiferBoundaryCondition::KernelWrapper AquiferBoundaryCondition::createKernelWrapper() const
{
  FunctionManager const & functionManager = FunctionManager::getInstance();
  TableFunction const & pressureInfluenceFunction = functionManager.getGroup< TableFunction >( m_pressureInfluenceFunctionName );

  return AquiferBoundaryCondition::KernelWrapper( m_initialPressure,
                                                  m_density,
                                                  m_elevation * m_gravityVector[2],
                                                  m_timeConstant,
                                                  m_influxConstant,
                                                  m_cumulativeFlux,
                                                  pressureInfluenceFunction.createKernelWrapper() );
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, AquiferBoundaryCondition, string const &, Group * const )


} /* namespace geos */
