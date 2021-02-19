/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FracturePermeabilityBase.cpp
 */

#include "../permeability/FracturePermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


FracturePermeabilityBase::FracturePermeabilityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{}

FracturePermeabilityBase::~FracturePermeabilityBase() = default;

std::unique_ptr< ConstitutiveBase >
FracturePermeabilityBase::deliverClone( string const & name,
                                        Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void FracturePermeabilityBase::allocateConstitutiveData( dataRepository::Group * const parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_permeability.resize( parent->size(), numConstitutivePointsPerParentIndex, 1 );
  m_dPerm_dAperture.resize( parent->size(), numConstitutivePointsPerParentIndex, 1 );
}

void FracturePermeabilityBase::postProcessInput()
{}

}
} /* namespace geosx */
