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
 * @file PressurePorosityModel.cpp
 */

#include "PressurePorosityModel.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

PressurePorosityModel::PressurePorosityModel( string const & name, Group * const parent ):
  PorosityBase( name, parent ),
  m_referencePressure(),
  m_compressibility()
{
  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference pressure for solid compressibility" );

  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Solid compressibility" );
}

PressurePorosityModel::~PressurePorosityModel() = default;

std::unique_ptr< ConstitutiveBase >
PressurePorosityModel::deliverClone( string const & name,
                                Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void PressurePorosityModel::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PressurePorosityModel::postProcessInput()
{
  this->getWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() ).
    setApplyDefaultValue( m_defaultReferencePorosity );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PressurePorosityModel, string const &, Group * const )
}
} /* namespace geosx */
