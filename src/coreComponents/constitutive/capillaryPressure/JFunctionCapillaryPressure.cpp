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
 * @file JFunctionCapillaryPressure.cpp
 */

#include "JFunctionCapillaryPressure.hpp"

#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressureHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

JFunctionCapillaryPressure::JFunctionCapillaryPressure( std::string const & name,
                                                        Group * const parent )
  : CapillaryPressureBase( name, parent ),
  m_wettingNonWettingSurfaceTension( -1 ),
  m_wettingIntermediateSurfaceTension( -1 ),
  m_nonWettingIntermediateSurfaceTension( -1 )
{
  registerWrapper( viewKeyStruct::wettingNonWettingJFuncTableNameString(), &m_wettingNonWettingJFuncTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "J-function table (dimensionless) for the pair (wetting phase, non-wetting phase)\n"
                    "Note that this input is only used for two-phase flow.\n"
                    "If you want to do a three-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingIntermediateJFuncTableNameString() ) +
                    " and " +
                    string( viewKeyStruct::nonWettingIntermediateJFuncTableNameString() ) +
                    " to specify the table names." );

  registerWrapper( viewKeyStruct::wettingIntermediateJFuncTableNameString(), &m_wettingIntermediateJFuncTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "J-function table (dimensionless) for the pair (wetting phase, intermediate phase)\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingJFuncTableNameString() ) +
                    " to specify the table names." );

  registerWrapper( viewKeyStruct::nonWettingIntermediateJFuncTableNameString(), &m_nonWettingIntermediateJFuncTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "J-function table (dimensionless) for the pair (non-wetting phase, intermediate phase)\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingJFuncTableNameString() ) +
                    " to specify the table names." );

  registerWrapper( viewKeyStruct::wettingNonWettingSurfaceTensionString(), &m_wettingNonWettingSurfaceTension ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface tension [N/m] for the pair (wetting phase, non-wetting phase)\n"
                    "If you have a value in [dyne/cm], divide it by 1000 to obtain the value in [N/m]\n"
                    "Note that this input is only used for two-phase flow.\n"
                    "If you want to do a three-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingIntermediateSurfaceTensionString() ) +
                    " and " +
                    string( viewKeyStruct::nonWettingIntermediateSurfaceTensionString() ) +
                    " to specify the surface tensions." );

  registerWrapper( viewKeyStruct::wettingIntermediateSurfaceTensionString(), &m_wettingIntermediateSurfaceTension ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface tension [N/m] for the pair (wetting phase, intermediate phase)\n"
                    "If you have a value in [dyne/cm], divide it by 1000 to obtain the value in [N/m]\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingSurfaceTensionString() ) +
                    " to specify the surface tensions." );

  registerWrapper( viewKeyStruct::nonWettingIntermediateSurfaceTensionString(), &m_nonWettingIntermediateSurfaceTension ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface tension [N/m] for the pair (non-wetting phase, intermediate phase)\n"
                    "If you have a value in [dyne/cm], divide it by 1000 to obtain the value in [N/m]\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingSurfaceTensionString() ) +
                    " to specify the surface tensions." );

  registerWrapper( viewKeyStruct::porosityExponentString(), &m_porosityExponent ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Porosity exponent" );

  registerWrapper( viewKeyStruct::permeabilityExponentString(), &m_permeabilityExponent ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Permeability exponent" );

  registerWrapper( viewKeyStruct::permeabilityDirectionString(), &m_permeabilityDirection ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Permeability direction. Options are:\n" +
                    toString( PermeabilityDirection::XY ) + " - use the average of the permeabilities in the x and y directions,\n" +
                    toString( PermeabilityDirection::X ) + " - only use the permeability in the x direction,\n" +
                    toString( PermeabilityDirection::Y ) + " - only use the permeability in the y direction,\n" +
                    toString( PermeabilityDirection::Z ) + " - only use the permeability in the z direction." );

  registerField( fields::cappres::jFuncMultiplier{}, &m_jFuncMultiplier );

  registerWrapper( viewKeyStruct::jFunctionWrappersString(), &m_jFuncKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void JFunctionCapillaryPressure::postInputInitialization()
{
  CapillaryPressureBase::postInputInitialization();

  integer const numPhases = m_phaseNames.size();
  GEOS_THROW_IF( numPhases != 2 && numPhases != 3,
                 GEOS_FMT( "{}: the expected number of fluid phases is either two, or three",
                           getFullName() ),
                 InputError );

  if( numPhases == 2 )
  {
    GEOS_THROW_IF( m_wettingNonWettingJFuncTableName.empty(),
                   GEOS_FMT( "{}: for a two-phase flow simulation, we must use {} to specify the J-function table for the pair (wetting phase, non-wetting phase)",
                             getFullName(),
                             viewKeyStruct::wettingNonWettingJFuncTableNameString() ),
                   InputError );
    GEOS_THROW_IF( m_wettingNonWettingSurfaceTension <= 0,
                   GEOS_FMT( "{}: for a two-phase flow simulation, we must use {} to specify the surface tension for the pair (wetting phase, non-wetting phase)",
                             getFullName(),
                             viewKeyStruct::wettingNonWettingSurfaceTensionString() ),
                   InputError );
  }
  else if( numPhases == 3 )
  {
    GEOS_THROW_IF( m_wettingIntermediateJFuncTableName.empty() || m_nonWettingIntermediateJFuncTableName.empty(),
                   GEOS_FMT( "{}: for a three-phase flow simulation, we must use {} to specify the J-function table"
                             "for the pair (wetting phase, intermediate phase), "
                             "and {} to specify the J-function table for the pair (non-wetting phase, intermediate phase)",
                             getFullName(),
                             viewKeyStruct::wettingIntermediateJFuncTableNameString(),
                             viewKeyStruct::nonWettingIntermediateJFuncTableNameString()  ),
                   InputError );
    GEOS_THROW_IF( m_wettingIntermediateSurfaceTension <= 0 || m_nonWettingIntermediateSurfaceTension <= 0,
                   GEOS_FMT( "{}: for a three-phase flow simulation, we must use {} to specify the surface tension"
                             "for the pair (wetting phase, intermediate phase), "
                             "and {} to specify the J-function table for the pair (non-wetting phase, intermediate phase)",
                             getFullName(),
                             viewKeyStruct::wettingIntermediateSurfaceTensionString(),
                             viewKeyStruct::nonWettingIntermediateSurfaceTensionString()  ),
                   InputError );
  }
}

void JFunctionCapillaryPressure::initializePreSubGroups()
{
  CapillaryPressureBase::initializePreSubGroups();

  integer const numPhases = m_phaseNames.size();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  if( numPhases == 2 )
  {
    GEOS_THROW_IF( !functionManager.hasGroup( m_wettingNonWettingJFuncTableName ),
                   GEOS_FMT( "{}: the table function named {} could not be found",
                             getFullName(),
                             m_wettingNonWettingJFuncTableName ),
                   InputError );
    TableFunction const & jFuncTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingJFuncTableName );
    bool const jFuncMustBeIncreasing = ( m_phaseOrder[PhaseType::WATER] < 0 )
      ? true   // pc on the gas phase, function must be increasing
      : false; // pc on the water phase, function must be decreasing
    TableCapillaryPressureHelpers::validateCapillaryPressureTable( jFuncTable, getFullName(), jFuncMustBeIncreasing );
  }
  else if( numPhases == 3 )
  {
    GEOS_THROW_IF( !functionManager.hasGroup( m_wettingIntermediateJFuncTableName ),
                   GEOS_FMT( "{}: the table function named {} could not be found",
                             getFullName(),
                             m_wettingIntermediateJFuncTableName ),
                   InputError );
    TableFunction const & jFuncTableWI = functionManager.getGroup< TableFunction >( m_wettingIntermediateJFuncTableName );
    TableCapillaryPressureHelpers::validateCapillaryPressureTable( jFuncTableWI, getFullName(), false );

    GEOS_THROW_IF( !functionManager.hasGroup( m_nonWettingIntermediateJFuncTableName ),
                   GEOS_FMT( "{}: the table function named {} could not be found",
                             getFullName(),
                             m_nonWettingIntermediateJFuncTableName ),
                   InputError );
    TableFunction const & jFuncTableNWI = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateJFuncTableName );
    TableCapillaryPressureHelpers::validateCapillaryPressureTable( jFuncTableNWI, getFullName(), true );
  }
}

void JFunctionCapillaryPressure::initializeRockState( arrayView2d< real64 const > const & initialPorosity,
                                                      arrayView3d< real64 const > const & initialPermeability ) const
{
  // in this function, we compute and save the cell-wise J-function multipliers as a function of porosity and permeability
  saveConvergedRockState( initialPorosity, initialPermeability );
}

void JFunctionCapillaryPressure::saveConvergedRockState( arrayView2d< real64 const > const & convergedPorosity,
                                                         arrayView3d< real64 const > const & convergedPermeability ) const
{
  // in this function, we compute and save the cell-wise J-function multipliers as a function of porosity and permeability
  // we assume an explicit dependence of capillary pressure on porosity and permeability

  localIndex const numElems = m_phaseCapPressure.size( 0 );
  integer const numPhases = numFluidPhases();

  real64 const porosityExponent = m_porosityExponent;
  real64 const permeabilityExponent = m_permeabilityExponent;
  PermeabilityDirection const permeabilityDirection = m_permeabilityDirection;

  real64 const wettingNonWettingSurfaceTension = m_wettingNonWettingSurfaceTension;
  real64 const wettingIntermediateSurfaceTension = m_wettingIntermediateSurfaceTension;
  real64 const nonWettingIntermediateSurfaceTension = m_nonWettingIntermediateSurfaceTension;

  arrayView2d< real64 > jFuncMultiplier = m_jFuncMultiplier;

  forAll< parallelDevicePolicy<> >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 permeability = 0;
    if( permeabilityDirection == PermeabilityDirection::XY )
    {
      permeability = 0.5 * ( convergedPermeability[ei][0][0] + convergedPermeability[ei][0][1] );
    }
    else if( permeabilityDirection == PermeabilityDirection::X )
    {
      permeability = convergedPermeability[ei][0][0];
    }
    else if( permeabilityDirection == PermeabilityDirection::Y )
    {
      permeability = convergedPermeability[ei][0][1];
    }
    else if( permeabilityDirection == PermeabilityDirection::Z )
    {
      permeability = convergedPermeability[ei][0][2];
    }

    // here we compute an average of the porosity over quadrature points
    // this average is exact for tets, regular pyramids/wedges/hexes, or for VEM
    real64 porosityAveragedOverQuadraturePoints = 0;
    for( integer i = 0; i < convergedPorosity.size( 1 ); ++i )
    {
      porosityAveragedOverQuadraturePoints += convergedPorosity[ei][i];
    }
    porosityAveragedOverQuadraturePoints /= convergedPorosity.size( 1 );
    porosityAveragedOverQuadraturePoints =
      LvArray::math::max( porosityAveragedOverQuadraturePoints, LvArray::NumericLimits< real64 >::epsilon );

    real64 const porosityOverPermeability = pow( porosityAveragedOverQuadraturePoints, porosityExponent )
                                            / pow( permeability, permeabilityExponent );

    // units:
    // ------
    // porosity []
    // permeability [m^2]
    // surfaceTension [N/m]
    // multiplier [N/m]/([m^2])^1/2 = [N/m^2] = [Pa]

    if( numPhases == 2 )
    {
      jFuncMultiplier[ei][0] = wettingNonWettingSurfaceTension * porosityOverPermeability;
    }
    else if( numPhases == 3 )
    {
      using TPT = JFunctionCapillaryPressure::ThreePhasePairPhaseType;
      jFuncMultiplier[ei][TPT::INTERMEDIATE_WETTING] = wettingIntermediateSurfaceTension * porosityOverPermeability;
      jFuncMultiplier[ei][TPT::INTERMEDIATE_NONWETTING] = nonWettingIntermediateSurfaceTension * porosityOverPermeability;
    }
  } );
}


void JFunctionCapillaryPressure::createAllTableKernelWrappers()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  integer const numPhases = m_phaseNames.size();

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_jFuncKernelWrappers.clear();
  if( numPhases == 2 )
  {
    TableFunction const & jFuncTable = functionManager.getGroup< TableFunction >( m_wettingNonWettingJFuncTableName );
    m_jFuncKernelWrappers.emplace_back( jFuncTable.createKernelWrapper() );
  }
  else if( numPhases == 3 )
  {
    // the assumption used everywhere in this class is that the WI information comes before the NWI information
    TableFunction const & jFuncTableWI = functionManager.getGroup< TableFunction >( m_wettingIntermediateJFuncTableName );
    m_jFuncKernelWrappers.emplace_back( jFuncTableWI.createKernelWrapper() );
    TableFunction const & jFuncTableNWI = functionManager.getGroup< TableFunction >( m_nonWettingIntermediateJFuncTableName );
    m_jFuncKernelWrappers.emplace_back( jFuncTableNWI.createKernelWrapper() );
  }
}

JFunctionCapillaryPressure::KernelWrapper::
  KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & jFuncKernelWrappers,
                 arrayView2d< real64 const > const & jFuncMultiplier,
                 arrayView1d< integer const > const & phaseTypes,
                 arrayView1d< integer const > const & phaseOrder,
                 arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPres,
                 arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPres_dPhaseVolFrac )
  : CapillaryPressureBaseUpdate( phaseTypes,
                                 phaseOrder,
                                 phaseCapPres,
                                 dPhaseCapPres_dPhaseVolFrac ),
  m_jFuncKernelWrappers( jFuncKernelWrappers ),
  m_jFuncMultiplier( jFuncMultiplier )
{}

JFunctionCapillaryPressure::KernelWrapper
JFunctionCapillaryPressure::createKernelWrapper()
{
  createAllTableKernelWrappers();
  return KernelWrapper( m_jFuncKernelWrappers,
                        m_jFuncMultiplier,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

void JFunctionCapillaryPressure::allocateConstitutiveData( dataRepository::Group & parent,
                                                           localIndex const numConstitutivePointsPerParentIndex )
{
  m_jFuncMultiplier.resize( parent.size(), numFluidPhases()-1 );
  CapillaryPressureBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, JFunctionCapillaryPressure, std::string const &, Group * const )

} // namespace constitutive

} // namespace geos
