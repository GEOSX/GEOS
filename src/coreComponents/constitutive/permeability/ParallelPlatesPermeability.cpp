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

/**
 * @file ParallelPlatesPermeability.cpp
 */

#include "ParallelPlatesPermeability.hpp"

#include "constitutive/permeability/PermeabilityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


ParallelPlatesPermeability::ParallelPlatesPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent ),
  m_updateTransversalComponent( true )
{
  registerWrapper( viewKeyStruct::transversalPermeabilityString(), &m_transversalPermeability ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Default value of the permeability normal to the surface. If not specified the permeability is updated using the cubic law. " );

  registerField( fields::permeability::dPerm_dDispJump{}, &m_dPerm_dDispJump );
}

std::unique_ptr< ConstitutiveBase >
ParallelPlatesPermeability::deliverClone( string const & name,
                                          Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void ParallelPlatesPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                           localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dDispJump.resize( 0, 1, 3, 3 );

  if( m_transversalPermeability > -1 )
  {
    m_updateTransversalComponent = false;
  }

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void ParallelPlatesPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();

  real64 const transversalPerm = m_transversalPermeability;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      permView[ei][q][2] = transversalPerm;
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ParallelPlatesPermeability, string const &, Group * const )

}
} /* namespace geos */
