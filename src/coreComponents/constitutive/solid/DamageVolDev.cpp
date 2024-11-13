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
 * @file Damage.cpp
 */

#include "Damage.hpp"
#include "DamageVolDev.hpp"

#include "ElasticIsotropic.hpp"

namespace geos
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
DamageVolDev< BASE >::DamageVolDev( string const & name, Group * const parent ):
  Damage< BASE >( name, parent )
{}

template< typename BASE >
DamageVolDev< BASE >::~DamageVolDev()
{}

typedef DamageVolDev< ElasticIsotropic > DamageVolDevElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageVolDevElasticIsotropic, string const &, Group * const )

}
} /* namespace geos */
