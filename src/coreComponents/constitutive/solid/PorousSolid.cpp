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
 * @file PorousSolid.cpp
 */

#include "PorousSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "ElasticOrthotropic.hpp"
#include "DelftEgg.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "ModifiedCamClay.hpp"
#include "DuvautLionsSolid.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::PorousSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >( name, parent )
{}

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::~PorousSolid() = default;

template< typename SOLID_TYPE >
void PorousSolid< SOLID_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::initializeState();
}

// Register all PorousSolid model types.
typedef PorousSolid< ElasticIsotropic > PorousElasticIsotropic;
typedef PorousSolid< ElasticTransverseIsotropic > PorousElasticTransverseIsotropic;
typedef PorousSolid< ElasticOrthotropic > PorousElasticOrthotropic;
typedef PorousSolid< DelftEgg > PorousDelftEgg;
typedef PorousSolid< DruckerPrager > PorousDruckerPrager;
typedef PorousSolid< DruckerPragerExtended > PorousDruckerPragerExtended;
// typedef PorousSolid< Damage< ElasticIsotropic > > PorousDamageElasticIsotropic;
// typedef PorousSolid< DamageSpectral< ElasticIsotropic > > PorousDamageSpectralElasticIsotropic;
// typedef PorousSolid< DamageVolDev< ElasticIsotropic > > PorousDamageVolDevElasticIsotropic;
typedef PorousSolid< DuvautLionsSolid< DruckerPrager > > PorousViscoDruckerPrager;
typedef PorousSolid< DuvautLionsSolid< DruckerPragerExtended > > PorousViscoDruckerPragerExtended;
typedef PorousSolid< DuvautLionsSolid< ModifiedCamClay > > PorousViscoModifiedCamClay;
typedef PorousSolid< ModifiedCamClay > PorousModifiedCamClay;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticOrthotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDelftEgg, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerExtended, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousModifiedCamClay, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPragerExtended, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoModifiedCamClay, string const &, Group * const )


}
} /* namespace geos */
