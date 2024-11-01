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
 * @file ProppantSolid.cpp
 */

#include "ProppantSolid.hpp"
#include "porosity/ProppantPorosity.hpp"
#include "constitutive/permeability/ProppantPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename PORO_TYPE,
          typename PERM_TYPE >
ProppantSolid< PORO_TYPE, PERM_TYPE >::ProppantSolid( string const & name, Group * const parent ):
  CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >( name, parent )
{}

template< typename PORO_TYPE,
          typename PERM_TYPE >
ProppantSolid< PORO_TYPE, PERM_TYPE >::~ProppantSolid() = default;

typedef ProppantSolid< ProppantPorosity, ProppantPermeability > ProppantSolidModel;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantSolidModel, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
