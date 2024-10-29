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
 * @file PermeabilityBase.cpp
 */

#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


PermeabilityBase::PermeabilityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_permeability(),
  m_dPerm_dPressure()
{
  registerField( fields::permeability::permeability{}, &m_permeability );
  registerField( fields::permeability::dPerm_dPressure{}, &m_dPerm_dPressure );
}

std::unique_ptr< ConstitutiveBase >
PermeabilityBase::deliverClone( string const & name,
                                Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void PermeabilityBase::scaleHorizontalPermeability( arrayView1d< real64 const > scalingFactors ) const
{
  localIndex const numElems = m_permeability.size( 0 );
  integer const numQuad = 1; // NOTE: enforcing 1 quadrature point
  for( localIndex ei = 0; ei < numElems; ++ei )
  {
    for( integer q = 0; q < numQuad; ++q )
    {
      m_permeability[ei][q][0] *= scalingFactors[ei];
      m_permeability[ei][q][1] *= scalingFactors[ei];
    }
  }
}

void PermeabilityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_permeability.resize( 0, 1, 3 );
  m_dPerm_dPressure.resize( 0, 1, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

}
} /* namespace geos */
