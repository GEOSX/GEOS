
/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DuvautLionsSolid.cpp
 */

#include "DuvautLionsSolid.hpp"
#include "ElasticIsotropic.hpp"
//#include "DelftEgg.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "ModifiedCamClay.hpp"

namespace geos
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
DuvautLionsSolid< BASE >::DuvautLionsSolid( string const & name, Group * const parent ):
  BASE( name, parent ),
  m_relaxationTime()
{

  this->registerWrapper( viewKeyStruct::relaxationTimeString(), &m_relaxationTime ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Relaxation time" );
}


template< typename BASE >
void DuvautLionsSolid< BASE >::postInputInitialization()
{
  BASE::postInputInitialization();
}

template< typename BASE >
void DuvautLionsSolid< BASE >::allocateConstitutiveData( dataRepository::Group & parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  BASE::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

//typedef DuvautLionsSolid< ElasticIsotropic > ViscoElasticIsotropic;
typedef DuvautLionsSolid< DruckerPrager > ViscoDruckerPrager;
typedef DuvautLionsSolid< DruckerPragerExtended > ViscoDruckerPragerExtended;
//typedef DuvautLionsSolid< DelftEgg > ViscoDelftEgg;
typedef DuvautLionsSolid< ModifiedCamClay > ViscoModifiedCamClay;

//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDruckerPragerExtended, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDelftEgg, string const &, Group * const )
//REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ViscoModifiedCamClay, string const &, Group * const )

}
} /* namespace geos */
