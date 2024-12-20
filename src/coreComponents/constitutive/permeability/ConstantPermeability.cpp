/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ConstantPermeability.cpp
 */

#include "ConstantPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


ConstantPermeability::ConstantPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent ),
  m_diagonalPermeabilityTensor{-1.0, -1.0, -1.0},
  m_symmetricFullPermeabilityTensor{-1.0, -1.0, -1.0, 0.0, 0.0, 0.0}
{
  registerWrapper( viewKeyStruct::diagonalPermeabilityTensorString(), &m_diagonalPermeabilityTensor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy and zz components of a diagonal permeability tensor." );

  registerWrapper( viewKeyStruct::symmetricFullPermeabilityTensorString(), &m_symmetricFullPermeabilityTensor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy and zz, xz, yz, xy components of a symmetric permeability tensor." );
}

void ConstantPermeability::postProcessInput()
{
  GEOS_ERROR_IF( m_diagonalPermeabilityTensor[0] < 0.0 && m_symmetricFullPermeabilityTensor[0] < 0.0, 
  "Either a diagonal permeability tensor or a full tensor must be provided.");

  GEOS_ERROR_IF( m_diagonalPermeabilityTensor[0] > 0.0 && m_symmetricFullPermeabilityTensor[0] > 0.0, 
  "Only one between a diagonal permeability tensor and a full tensor permeability can be provided.");

  if ( m_diagonalPermeabilityTensor[0] > 0.0 )
  {
    for( int i = 0; i < 3; i++ )
    { 
      m_symmetricFullPermeabilityTensor[i] = m_diagonalPermeabilityTensor[i]; 
    }
    for( int j=4; j < 6; j++ )
    { 
      m_symmetricFullPermeabilityTensor[j] = 0.0; 
    }
  }
}

std::unique_ptr< ConstitutiveBase >
ConstantPermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void ConstantPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  integer const numQuad = 1; // NOTE: enforcing 1 quadrature point

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for ( int c=0; c < 6; c++)
      m_permeability[ei][q][c] =  m_symmetricFullPermeabilityTensor[c];
    }
  }
}

void ConstantPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();
  real64 const permComponents[6] = { m_symmetricFullPermeabilityTensor[0],
                                     m_symmetricFullPermeabilityTensor[1],
                                     m_symmetricFullPermeabilityTensor[2],
                                     m_symmetricFullPermeabilityTensor[3],
                                     m_symmetricFullPermeabilityTensor[4],
                                     m_symmetricFullPermeabilityTensor[5] };

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim=0; dim < 6; ++dim )
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

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ConstantPermeability, string const &, Group * const )

}
} /* namespace geos */
