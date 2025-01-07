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
 * @file FrictionSelector.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_CONTACTSELECTOR_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_CONTACTSELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/contact/HydraulicApertureTable.hpp"
#include "constitutive/contact/BartonBandis.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( HydraulicApertureBase const & contact,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< HydraulicApertureTable,
                               BartonBandis >::execute( contact, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( HydraulicApertureBase & contact,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< HydraulicApertureTable,
                               BartonBandis >::execute( contact, std::forward< LAMBDA >( lambda ) );
}

} /* namespace constitutive */

} /* namespace geos */

#endif // GEOS_CONSTITUTIVE_CONTACT_CONTACTSELECTOR_HPP_
