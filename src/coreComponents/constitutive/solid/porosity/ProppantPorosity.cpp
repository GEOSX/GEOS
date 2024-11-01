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
 * @file ProppantPorosity.cpp
 */

#include "constitutive/solid/porosity/ProppantPorosity.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ProppantPorosity::ProppantPorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent ),
  m_maxProppantConcentration()
{
  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum proppant concentration " );
}

ProppantPorosity::~ProppantPorosity() = default;

std::unique_ptr< ConstitutiveBase >
ProppantPorosity::deliverClone( string const & name,
                                Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void ProppantPorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void ProppantPorosity::postInputInitialization()
{
  getField< fields::porosity::referencePorosity >().
    setApplyDefaultValue( m_defaultReferencePorosity );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantPorosity, string const &, Group * const )
}
} /* namespace geos */
