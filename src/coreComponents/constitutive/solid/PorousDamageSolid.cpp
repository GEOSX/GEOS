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
 * @file PorousDamageSolid.cpp
 */

#include "PorousDamageSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "Damage.hpp"
#include "DamageSpectral.hpp"
#include "DamageVolDev.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE >
PorousDamageSolid< SOLID_TYPE >::PorousDamageSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, DamagePermeability >( name, parent )
{}

template< typename SOLID_TYPE >
PorousDamageSolid< SOLID_TYPE >::~PorousDamageSolid() = default;

// Register all PorousDamageSolid model types.
typedef PorousDamageSolid< Damage< ElasticIsotropic > > PorousDamageElasticIsotropic;
typedef PorousDamageSolid< DamageSpectral< ElasticIsotropic > > PorousDamageSpectralElasticIsotropic;
typedef PorousDamageSolid< DamageVolDev< ElasticIsotropic > > PorousDamageVolDevElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDamageElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDamageSpectralElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDamageVolDevElasticIsotropic, string const &, Group * const )

}
} /* namespace geos */
