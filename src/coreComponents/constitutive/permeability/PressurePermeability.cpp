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
 * @file PressurePermeability.cpp
 */

#include "PressurePermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


PressurePermeability::PressurePermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::referencePermeabilityComponentsString(), &m_referencePermeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Reference xx, yy and zz components of a diagonal permeability tensor." );

  registerWrapper( viewKeyStruct::pressureDependenceConstantsString(), &m_pressureDependenceConstants ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Pressure dependence coefficients for each permeability component." );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference pressure for the pressure permeability model" );

  registerWrapper( viewKeyStruct::referencePermeabilityString(), &m_referencePermeability ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Reference permeability field" );

  registerWrapper( viewKeyStruct::maxPermeabilityString(), &m_maxPermeability ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max. permeability can be reached." );

  registerWrapper( viewKeyStruct::pressureModelTypeString(), &m_presModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( PressureModelType::Hyperbolic ).
    setDescription( "Type of the pressure dependence model. " );
}

std::unique_ptr< ConstitutiveBase >
PressurePermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void PressurePermeability::postInputInitialization()
{
  for( localIndex i=0; i < 3; i++ )
  {
    GEOS_ERROR_IF( fabs( m_pressureDependenceConstants[i] ) < 1e-15 && m_presModelType == PressureModelType::Hyperbolic,
                   getDataContext() << ": the pressure dependent constant at component " << i << " is too close to zero, which is not allowed for the hyperbolic model." );
  }
}

void PressurePermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  m_referencePermeability.resize( 0, 1, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  integer const numQuad = 1; // NOTE: enforcing 1 quadrature point

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      m_referencePermeability[ei][q][0] =  m_referencePermeabilityComponents[0];
      m_referencePermeability[ei][q][1] =  m_referencePermeabilityComponents[1];
      m_referencePermeability[ei][q][2] =  m_referencePermeabilityComponents[2];
    }
  }
}

void PressurePermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();
  real64 const permComponents[3] = { m_referencePermeabilityComponents[0],
                                     m_referencePermeabilityComponents[1],
                                     m_referencePermeabilityComponents[2] };

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim=0; dim < 3; ++dim )
      {
        // The default value is -1 so if it still -1 it needs to be set to something physical
        if( permView[ei][q][dim] < 0 )
        {
          permView[ei][q][dim] =  permComponents[dim];
        }
      }
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PressurePermeability, string const &, Group * const )

}
} /* namespace geos */
