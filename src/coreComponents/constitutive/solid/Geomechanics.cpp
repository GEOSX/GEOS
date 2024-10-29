/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file Geomechanics.cpp
 */

#include "Geomechanics.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

Geomechanics::Geomechanics( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_b0( 0.0 ),
  m_b1( 0.0 ),
  m_b2( 0.0 ),
  m_b3( 0.0 ),
  m_b4( 0.0 ),
  m_g0( 0.0 ),
  m_g1( 0.0 ),
  m_g2( 0.0 ),
  m_g3( 0.0 ),
  m_g4( 0.0 ),
  m_p0( 0.0 ),
  m_p1( 0.0 ),
  m_p2( 0.0 ),
  m_p3( 0.0 ),
  m_p4( 0.0 ),
  m_peakT1( 0.0 ),
  m_fSlope( 0.0 ),
  m_stren( 0.0 ),
  m_ySlope( 0.0 ),
  m_beta( 1.0 ),
  m_t1RateDependence( 0.0 ),
  m_t2RateDependence( 0.0 ),
  m_fractureEnergyReleaseRate( 0.0 ),
  m_cr( 0.0 ),
  m_fluidBulkModulus(0.0 ),
  m_fluidInitialPressure( 0.0 ),
  m_creep( 0 ),
  m_creepC0( 0.0),
  m_creepC1( 0.0 ),
  m_creepA( 0.0 ),
  m_creepB( 0.0 ),
  m_creepC( 0.0 ),
  m_strainHardeningN( 0.0 ),
  m_strainHardeningK( 0.0 ),
  m_bulkModulus(),
  m_shearModulus(),
  m_velocityGradient(),
  m_plasticStrain(),
  m_damage(),
  m_lengthScale()
{
  // register default values
  registerWrapper( viewKeyStruct::b0String(), &m_b0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 0" );

  registerWrapper( viewKeyStruct::b1String(), &m_b1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 1" );

  registerWrapper( viewKeyStruct::b2String(), &m_b2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 2" );

  registerWrapper( viewKeyStruct::b3String(), &m_b3 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 3" );

  registerWrapper( viewKeyStruct::b4String(), &m_b4 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic bulk modulus parameter 4" );

  registerWrapper( viewKeyStruct::g0String(), &m_g0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent shear shear modulus parameter 0" );

  registerWrapper( viewKeyStruct::g1String(), &m_g1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 1" );

  registerWrapper( viewKeyStruct::g2String(), &m_g2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 2" );

  registerWrapper( viewKeyStruct::g3String(), &m_g3 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 3" );

  registerWrapper( viewKeyStruct::g4String(), &m_g4 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default tangent elastic shear modulus parameter 4" );

  registerWrapper( viewKeyStruct::p0String(), &m_p0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 0" );

  registerWrapper( viewKeyStruct::p1String(), &m_p1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 1" );

  registerWrapper( viewKeyStruct::p2String(), &m_p2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 2" );

  registerWrapper( viewKeyStruct::p3String(), &m_p3 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 3" );

  registerWrapper( viewKeyStruct::p4String(), &m_p4 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crush curve parameter 4" );

  registerWrapper( viewKeyStruct::crString(), &m_cr ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Cap shape parameter" );

  registerWrapper( viewKeyStruct::fluidBulkModulusString(), &m_fluidBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fluid bulk modulus" );

  registerWrapper( viewKeyStruct::fluidInitialPressureString(), &m_fluidInitialPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fluid initial pressure" );

  registerWrapper( viewKeyStruct::t1RateDependenceString(), &m_t1RateDependence ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Rate dependence parameter 1" );

  registerWrapper( viewKeyStruct::t2RateDependenceString(), &m_t2RateDependence ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Rate dependence parameter 2" );

  registerWrapper( viewKeyStruct::peakT1String(), &m_peakT1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Peak T1 shear limit parameter" );

  registerWrapper( viewKeyStruct::fSlopeString(), &m_fSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "F slope shear limit parameter" );

  registerWrapper( viewKeyStruct::strenString(), &m_stren ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Stren shear limit parameter" );

  registerWrapper( viewKeyStruct::ySlopeString(), &m_ySlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Y slope shear limit parameter" );

  registerWrapper( viewKeyStruct::betaString(), &m_beta ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Nonassociativity parameter" );

  registerWrapper( viewKeyStruct::creepString(), &m_creep ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep flag" );

  registerWrapper( viewKeyStruct::creepC0String(), &m_creepC0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep C0 parameter" );

  registerWrapper( viewKeyStruct::creepC1String(), &m_creepC1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep C1 parameter" );

  registerWrapper( viewKeyStruct::creepAString(), &m_creepA ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep A parameter" );
  
  registerWrapper( viewKeyStruct::creepBString(), &m_creepB ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep B parameter" );

  registerWrapper( viewKeyStruct::creepCString(), &m_creepC ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Creep C parameter" );

  registerWrapper( viewKeyStruct::strainHardeningNString(), &m_strainHardeningN ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strain Hardening n parameter" );

  registerWrapper( viewKeyStruct::strainHardeningKString(), &m_strainHardeningK ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strain Hardening K parameter" );

  // register fields
  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Bulk modulus" );
  
  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setInputFlag( InputFlags::FALSE).
    setDescription( "Shear modulus");

  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Velocity gradient" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Plastic strain" );

  registerWrapper( viewKeyStruct::porosityString(), &m_porosity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Porosity" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setApplyDefaultValue( DBL_MIN ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point length scale values" );
}


Geomechanics::~Geomechanics()
{}


void Geomechanics::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_bulkModulus.resize( 0 );
  m_shearModulus.resize( 0 );
  m_velocityGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_porosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
}


void Geomechanics::postInputInitialization()
{
    SolidBase::postInputInitialization();
    GEOS_THROW_IF( m_b0 <= 0.0, "b0 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_b1 <= 0.0, "b1 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_b2 <= 0.0, "b2 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_b3 <= 0.0, "b3 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_b4 <= 0.0, "b4 must be greater than 0", InputError );

    GEOS_THROW_IF( m_g0 < 0.0, "g0 must be greater than or equalt to 0", InputError );
    // GEOS_THROW_IF( m_g1 <= 0.0, "g1 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_g2 <= 0.0, "g2 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_g3 <= 0.0, "g3 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_g4 <= 0.0, "g4 must be greater than 0", InputError );

    GEOS_THROW_IF( m_p0 > 0.0, "p0 must be less than 0", InputError );
    // GEOS_THROW_IF( m_p1 <= 0.0, "p1 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_p2 <= 0.0, "p2 must be greater than 0", InputError );
    GEOS_THROW_IF( m_p3 <= 0.0, "p3 must be greater than 0", InputError );
    // GEOS_THROW_IF( m_p4 <= 0.0, "p4 must be greater than 0", InputError );

    // GEOS_THROW_IF( m_peakT1 <= 0.0, "peakT1 must be greater than 0", InputError );
    GEOS_THROW_IF( m_fSlope < 0.0, "fSlope must be greater than 0", InputError );
    // GEOS_THROW_IF( m_ySlope <= 0.0, "ySlope must be greater than 0", InputError );
    // GEOS_THROW_IF( m_stren <= 0.0, "stren must be greater than 0", InputError );
    GEOS_THROW_IF( m_beta <= 0.0, "beta must be greater than 0", InputError );
    // GEOS_THROW_IF( m_t1RateDependence <= 0.0, "t1RateDependence must be greater than 0", InputError );
    // GEOS_THROW_IF( m_t2RateDependence <= 0.0, "t2RateDependence must be greater than 0", InputError );
    // GEOS_THROW_IF( m_fractureEnergyReleaseRate <= 0.0, "fractureEnergyReleaseRate must be greater than 0", InputError );
    GEOS_THROW_IF( m_cr <= 0.0, "cr must be 0 < CR < 1", InputError );
    // GEOS_THROW_IF( m_fluidBulkModulus <= 0.0, "fluidBulkModulus must be greater than 0", InputError );
    // GEOS_THROW_IF( m_initialFluidPressure <= 0.0, "initialFluidPressure must be greater than 0", InputError );
}


void Geomechanics::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Geomechanics, std::string const &, Group * const )
}
} /* namespace geos */
