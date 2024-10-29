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
 * @file ProppantSlurryFluid.cpp
 */

#include "ProppantSlurryFluid.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ProppantSlurryFluid::ProppantSlurryFluid( string const & name, Group * const parent ):
  SlurryFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid compressibility" );

  registerWrapper( viewKeyStruct::referenceProppantDensityString(), &m_referenceProppantDensity ).
    setApplyDefaultValue( 1400.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference proppant density" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 1e5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::referenceDensityString(), &m_referenceDensity ).
    setApplyDefaultValue( 1000.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid density" );

  registerWrapper( viewKeyStruct::referenceViscosityString(), &m_referenceViscosity ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid viscosity" );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).
    setApplyDefaultValue( 0.6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum proppant concentration" );

}

ProppantSlurryFluid::~ProppantSlurryFluid() = default;

void ProppantSlurryFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  SlurryFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_density.setValues< serialPolicy >( m_referenceDensity );
  m_viscosity.setValues< serialPolicy >( m_referenceViscosity );
}


void ProppantSlurryFluid::postInputInitialization()
{
  SlurryFluidBase::postInputInitialization();

  GEOS_ERROR_IF_LT_MSG( m_compressibility, 0.0,
                        getFullName() << ": invalid value of " << viewKeyStruct::compressibilityString() );

  GEOS_ERROR_IF_LE_MSG( m_referenceDensity, 0.0,
                        getFullName() << ": invalid value of " << viewKeyStruct::referenceDensityString() );

  GEOS_ERROR_IF_LT_MSG( m_referenceViscosity, 0.0,
                        getFullName() << ": invalid value of " << viewKeyStruct::referenceViscosityString() );

  GEOS_ERROR_IF_LE_MSG( m_maxProppantConcentration, 0.0,
                        getFullName() << ": invalid value of " << viewKeyStruct::maxProppantConcentrationString() );

  GEOS_ERROR_IF_GT_MSG( m_maxProppantConcentration, 1.0,
                        getFullName() << ": invalid value of " << viewKeyStruct::maxProppantConcentrationString() );
}

ProppantSlurryFluid::KernelWrapper
ProppantSlurryFluid::createKernelWrapper()
{
  return KernelWrapper( m_compressibility,
                        m_referenceProppantDensity,
                        m_referencePressure,
                        m_referenceDensity,
                        m_referenceViscosity,
                        m_maxProppantConcentration,
                        m_defaultComponentDensity,
                        m_defaultComponentCompressibility,
                        m_defaultComponentViscosity,
                        m_nIndices,
                        m_Ks,
                        m_isNewtonianFluid,
                        m_density,
                        m_dDensity_dPressure,
                        m_dDensity_dProppantConc,
                        m_dDensity_dCompConc,
                        m_componentDensity,
                        m_dCompDens_dPres,
                        m_dCompDens_dCompConc,
                        m_fluidDensity,
                        m_dFluidDens_dPres,
                        m_dFluidDens_dCompConc,
                        m_fluidViscosity,
                        m_dFluidVisc_dPres,
                        m_dFluidVisc_dCompConc,
                        m_viscosity,
                        m_dViscosity_dPressure,
                        m_dViscosity_dProppantConc,
                        m_dViscosity_dCompConc );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantSlurryFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
