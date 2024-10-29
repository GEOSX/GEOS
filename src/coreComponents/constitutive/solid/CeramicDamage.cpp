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
 *  @file CeramicDamage.cpp
 */

#include "CeramicDamage.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

CeramicDamage::CeramicDamage( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_damage(),
  m_jacobian(),
  m_lengthScale(),
  m_strengthScale(),
  m_porosity(),
  m_referencePorosity(),
  m_tensileStrength(),
  m_compressiveStrength(),
  m_maximumStrength(),
  m_crackSpeed(),
  m_damagedMaterialFrictionSlope(),
  m_thirdInvariantDependence(),
  m_velocityGradient(),
  m_plasticStrain()
{
  // register default values
  registerWrapper( viewKeyStruct::tensileStrengthString(), &m_tensileStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tensile strength" );

  registerWrapper( viewKeyStruct::compressiveStrengthString(), &m_compressiveStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Compressive strength" );

  registerWrapper( viewKeyStruct::maximumStrengthString(), &m_maximumStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum theoretical strength" );

  registerWrapper( viewKeyStruct::crackSpeedString(), &m_crackSpeed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crack speed" );

  // register fields
  registerWrapper( viewKeyStruct::strengthScaleString(), &m_strengthScale ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::LEVEL_0).
    setDescription( "Strength scale" );

  registerWrapper( viewKeyStruct::porosityString(), &m_porosity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Porosity" );

  registerWrapper( viewKeyStruct::referencePorosityString(), &m_referencePorosity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Reference porosity" );

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
    setDescription( "Array of quadrature point damage values" );
  
  registerWrapper( viewKeyStruct::damagedMaterialFrictionSlopeString(), &m_damagedMaterialFrictionSlope ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5773502691896258 ).
    setDescription( "Value of the damaged material friction slope");

  registerWrapper( viewKeyStruct::thirdInvariantDependenceString(), &m_thirdInvariantDependence ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to enable third invariant dependence" );

  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Velocity gradient" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Plastic strain" );
}


CeramicDamage::~CeramicDamage()
{}


void CeramicDamage::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_strengthScale.resize( 0 );
  m_porosity.resize( 0 );
  m_referencePorosity.resize( 0 );
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_velocityGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
}


void CeramicDamage::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  // CC: TODO double check model inputs
  GEOS_THROW_IF( m_tensileStrength < 0.0, "Tensile strength must be a positive number.", InputError );
  GEOS_THROW_IF( m_compressiveStrength < m_tensileStrength, "Compressive strength must be greater than tensile strength.", InputError );
  GEOS_THROW_IF( m_maximumStrength < m_compressiveStrength, "Maximum theoretical strength must be greater than compressive strength.", InputError );
  GEOS_THROW_IF( m_crackSpeed < 0.0, "Crack speed must be a positive number.", InputError );
}


void CeramicDamage::saveConvergedState() const
{
  ElasticIsotropic::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, CeramicDamage, std::string const &, Group * const )
}
} /* namespace geos */
