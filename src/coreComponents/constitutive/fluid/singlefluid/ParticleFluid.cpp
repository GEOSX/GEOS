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
 * @file ParticleFluid.cpp
 */

#include "ParticleFluid.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ParticleFluid::ParticleFluid( string const & name, Group * const parent ):
  ParticleFluidBase( name, parent )
{

  registerWrapper( viewKeyStruct::particleSettlingModelString(), &m_particleSettlingModel ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Particle settling velocity model. Valid options:\n* " + EnumStrings< ParticleSettlingModel >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::proppantDensityString(), &m_proppantDensity ).
    setApplyDefaultValue( 1400.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant density" );

  registerWrapper( viewKeyStruct::fluidViscosityString(), &m_fluidViscosity ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid viscosity" );

  registerWrapper( viewKeyStruct::proppantDiameterString(), &m_proppantDiameter ).
    setApplyDefaultValue( 200e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant diameter" );

  registerWrapper( viewKeyStruct::hinderedSettlingCoefficientString(), &m_hinderedSettlingCoefficient ).
    setApplyDefaultValue( 5.9 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Hindered settling coefficient" );

  registerWrapper( viewKeyStruct::collisionAlphaString(), &m_collisionAlpha ).
    setApplyDefaultValue( 1.27 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Collision alpha coefficient" );

  registerWrapper( viewKeyStruct::slipConcentrationString(), &m_slipConcentration ).
    setApplyDefaultValue( 0.1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Slip concentration" );

  registerWrapper( viewKeyStruct::collisionBetaString(), &m_collisionBeta ).
    setApplyDefaultValue( 1.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Collision beta coefficient" );


  registerWrapper( viewKeyStruct::sphericityString(), &m_sphericity ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Sphericity" );

}

ParticleFluid::~ParticleFluid() = default;

void ParticleFluid::postInputInitialization()
{
  ParticleFluidBase::postInputInitialization();

  GEOS_ERROR_IF( m_proppantDensity < 500.0,
                 "Invalid proppantDensity in ParticleFluid "
                 << getDataContext() << ", which must >= 500.0 " );

  GEOS_ERROR_IF( m_proppantDiameter < 10e-6,
                 "Invalid proppantDiameter in ParticleFluid "
                 << getDataContext() << ", which must >= 10e-6 " );

  GEOS_ERROR_IF( m_hinderedSettlingCoefficient< 0.0 || m_hinderedSettlingCoefficient > 10.0,
                 "Invalid hinderedSettlingCoefficient in ParticleFluid "
                 << getDataContext() << ", which must between 0 and 10 " );

  GEOS_ERROR_IF( m_collisionAlpha < 1.0,
                 "Invalid collisionAlpha in ParticleFluid "
                 << getDataContext() << ", which must >= 1 " );

  GEOS_ERROR_IF( m_collisionBeta < 0.0,
                 "Invalid collisionBeta in ParticleFluid "
                 << getDataContext() << ", which must >= 0" );

  GEOS_ERROR_IF( m_slipConcentration > 0.3,
                 "Invalid slipConcentration in ParticleFluid "
                 << getDataContext() << ", which must <= 0.3" );

  m_packPermeabilityCoef = pow( m_sphericity * m_proppantDiameter, 2.0 ) / 180.0;
}

ParticleFluid::KernelWrapper
ParticleFluid::createKernelWrapper() const
{
  return KernelWrapper( m_particleSettlingModel,
                        m_proppantDensity,
                        m_proppantDiameter,
                        m_hinderedSettlingCoefficient,
                        m_collisionAlpha,
                        m_slipConcentration,
                        m_collisionBeta,
                        m_isCollisionalSlip,
                        m_maxProppantConcentration,
                        m_settlingFactor,
                        m_dSettlingFactor_dPressure,
                        m_dSettlingFactor_dProppantConcentration,
                        m_dSettlingFactor_dComponentConcentration,
                        m_collisionFactor,
                        m_dCollisionFactor_dProppantConcentration,
                        m_proppantPackPermeability );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ParticleFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
