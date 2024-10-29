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
 *  @file Graphite.cpp
 */

#include "Graphite.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

Graphite::Graphite( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultYoungModulusTransverse(),
  m_defaultYoungModulusAxial(),
  m_defaultPoissonRatioTransverse(),
  m_defaultPoissonRatioAxialTransverse(),
  m_defaultShearModulusAxialTransverse(),
  m_defaultYoungModulusTransversePressureDerivative(),
  m_defaultYoungModulusAxialPressureDerivative(),
  m_defaultShearModulusAxialTransversePressureDerivative(),
  m_velocityGradient(),
  m_plasticStrain(),
  m_relaxation(),
  m_damage(),
  m_jacobian(),
  m_lengthScale(),
  m_strengthScale(),
  m_failureStrength(),
  m_crackSpeed(),
  m_damagedMaterialFrictionalSlope(),
  m_distortionShearResponseX2(),
  m_distortionShearResponseY1(),
  m_distortionShearResponseY2(),
  m_distortionShearResponseM1(),
  m_inPlaneShearResponseX2(),
  m_inPlaneShearResponseY1(),
  m_inPlaneShearResponseY2(),
  m_inPlaneShearResponseM1(),
  m_coupledShearResponseX2(),
  m_coupledShearResponseY1(),
  m_coupledShearResponseY2(),
  m_coupledShearResponseM1(),
  m_maximumPlasticStrain()
{
  // register default values
  registerWrapper( viewKeyStruct::defaultYoungModulusTransverseString(), &m_defaultYoungModulusTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Transverse Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultYoungModulusAxialString(), &m_defaultYoungModulusAxial ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioTransverseString(), &m_defaultPoissonRatioTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioAxialTransverseString(), &m_defaultPoissonRatioAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial-Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverseString(), &m_defaultShearModulusAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial-Transverse Shear Modulus" );

  registerWrapper( viewKeyStruct::defaultYoungModulusTransversePressureDerivativeString(), &m_defaultYoungModulusTransversePressureDerivative ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Transverse Young's modulus pressure derivative" );

  registerWrapper( viewKeyStruct::defaultYoungModulusAxialPressureDerivativeString(), &m_defaultYoungModulusAxialPressureDerivative ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Axial Young's modulus pressure derivative" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransversePressureDerivativeString(), &m_defaultShearModulusAxialTransversePressureDerivative ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Axial transverse shear modulus pressure derivative" );

  registerWrapper( viewKeyStruct::failureStrengthString(), &m_failureStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum theoretical strength" );

  registerWrapper( viewKeyStruct::crackSpeedString(), &m_crackSpeed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crack speed" );

  registerWrapper( viewKeyStruct::damagedMaterialFrictionalSlopeString(), &m_damagedMaterialFrictionalSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Damaged material frictional slope" );

  registerWrapper( viewKeyStruct::distortionShearResponseX2String(), &m_distortionShearResponseX2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Distortion Shear Response X2" );

  registerWrapper( viewKeyStruct::distortionShearResponseY1String(), &m_distortionShearResponseY1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Distortion Shear Response Y1" );

  registerWrapper( viewKeyStruct::distortionShearResponseY2String(), &m_distortionShearResponseY2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Distortion Shear Response Y2" );

  registerWrapper( viewKeyStruct::distortionShearResponseM1String(), &m_distortionShearResponseM1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Distortion Shear Response M1" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseX2String(), &m_inPlaneShearResponseX2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response X2" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseY1String(), &m_inPlaneShearResponseY1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response Y1" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseY2String(), &m_inPlaneShearResponseY2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response Y2" );

  registerWrapper( viewKeyStruct::inPlaneShearResponseM1String(), &m_inPlaneShearResponseM1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "In Plane Shear Response M1" );

  registerWrapper( viewKeyStruct::coupledShearResponseX2String(), &m_coupledShearResponseX2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response X2" );

  registerWrapper( viewKeyStruct::coupledShearResponseY1String(), &m_coupledShearResponseY1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response Y1" );

  registerWrapper( viewKeyStruct::coupledShearResponseY2String(), &m_coupledShearResponseY2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response Y2" );

  registerWrapper( viewKeyStruct::coupledShearResponseM1String(), &m_coupledShearResponseM1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupled Shear Response M1" );

  registerWrapper( viewKeyStruct::maximumPlasticStrainString(), &m_maximumPlasticStrain ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum plastic strain" );

  // register fields
  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Velocity gradient" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Plastic strain" );

  registerWrapper( viewKeyStruct::relaxationString(), &m_relaxation).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Relaxation" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point jacobian values" );

  registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setApplyDefaultValue( DBL_MIN ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point length scale values" );

  registerWrapper( viewKeyStruct::strengthScaleString(), &m_strengthScale ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point strength scale values" );

  registerWrapper( viewKeyStruct::effectiveBulkModulusString(), &m_effectiveBulkModulus ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Effective bulk modulus for stress control and wavespeed calculations" );
  
  registerWrapper( viewKeyStruct::effectiveShearModulusString(), &m_effectiveShearModulus ).
    setInputFlag( InputFlags::FALSE).
    setDescription( "Effective shear modulus for stress control and wavespeed calculations");

  registerWrapper( viewKeyStruct::materialDirectionString(), &m_materialDirection ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Material direction" );
}


Graphite::~Graphite()
{}


void Graphite::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_effectiveBulkModulus.resize( 0 );
  m_effectiveShearModulus.resize( 0 );
  m_materialDirection.resize( 0, 3 );
  m_velocityGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_relaxation.resize( 0, numConstitutivePointsPerParentIndex );
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_lengthScale.resize( 0 );
  m_strengthScale.resize( 0 );
}


void Graphite::postInputInitialization()
{
  SolidBase::postInputInitialization();

  // GEOS_LOG_RANK_0( "Ez: " << m_defaultYoungModulusAxial << "\n" << 
  //                  "Ep: " << m_defaultYoungModulusTransverse << "\n" << 
  //                  "Nup: " << m_defaultPoissonRatioTransverse << "\n" << 
  //                  "Nuzp: " << m_defaultPoissonRatioAxialTransverse << "\n" << 
  //                  "Gzp: " << m_defaultShearModulusAxialTransverse << "\n" << 
  //                  "dEzdp: " << m_defaultYoungModulusAxialPressureDerivative << "\n" << 
  //                  "dEpdp: " << m_defaultYoungModulusTransversePressureDerivative << "\n" << 
  //                  "dGzpdp: " << m_defaultShearModulusAxialTransversePressureDerivative << "\n" << 
  //                  "sigmaFail: " << m_failureStrength << "\n" << 
  //                  "crackSpeed: " << m_crackSpeed << "\n" << 
  //                  "ds X2: " << m_distortionShearResponseX2 << "\n" << 
  //                  "ds Y1: " << m_distortionShearResponseY1 << "\n" << 
  //                  "ds Y2: " << m_distortionShearResponseY2 << "\n" << 
  //                  "ds M1: " << m_distortionShearResponseM1 << "\n" <<
  //                  "ips X2: " << m_inPlaneShearResponseX2 << "\n" << 
  //                  "ips Y1: " << m_inPlaneShearResponseY1 << "\n" << 
  //                  "ips Y2: " << m_inPlaneShearResponseY2 << "\n" << 
  //                  "ips M1: " << m_inPlaneShearResponseM1 << "\n" << 
  //                  "cs X2: " << m_coupledShearResponseX2 << "\n" << 
  //                  "cs Y1: " << m_coupledShearResponseY1 << "\n" << 
  //                  "cs Y2: " << m_coupledShearResponseY2 << "\n" << 
  //                  "cs M1: " << m_coupledShearResponseM1 << "\n" <<
  //                  "max ep: " << m_maximumPlasticStrain );

  // Add elastic constants check

  GEOS_THROW_IF( m_failureStrength <= 0.0, "Maximum theoretical strength must be greater than 0", InputError );
  GEOS_THROW_IF( m_crackSpeed <= 0.0, "Crack speed must be a positive number.", InputError );
  
  GEOS_THROW_IF( m_damagedMaterialFrictionalSlope < 0.0, "Damaged material frictional slope must be greater than 0", InputError );

  GEOS_THROW_IF( m_distortionShearResponseX2 < 0.0, "Distortion shear response x2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_distortionShearResponseY1 < 0.0, "Distortion shear response y1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_distortionShearResponseY2 < 0.0, "Distortion shear response y2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_distortionShearResponseM1 < 0.0, "Distortion shear response m1 must be a positive number.", InputError );

  GEOS_THROW_IF( m_inPlaneShearResponseX2 < 0.0, "In plane shear response x2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_inPlaneShearResponseY1 < 0.0, "In plane shear response y1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_inPlaneShearResponseY2 < 0.0, "In plane shear response y2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_inPlaneShearResponseM1 < 0.0, "In plane shear response m1 must be a positive number.", InputError );

  GEOS_THROW_IF( m_coupledShearResponseX2 < 0.0, "Coupled shear response x2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_coupledShearResponseY1 < 0.0, "Coupled shear response y1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_coupledShearResponseY2 < 0.0, "Coupled shear response y2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_coupledShearResponseM1 < 0.0, "Coupled shear response m1 must be a positive number.", InputError );

  GEOS_THROW_IF( m_maximumPlasticStrain < 0.0, "Maximum plastic strain must be a positive number.", InputError);
}


void Graphite::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Graphite, std::string const &, Group * const )
}
} /* namespace geos */
