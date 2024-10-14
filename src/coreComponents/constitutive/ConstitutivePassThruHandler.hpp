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
 * @file ConstitutivePassThruHandler.hpp
 */
#ifndef GEOS_CONSTITUTIVEPASSTHRUHANDLER_HPP
#define GEOS_CONSTITUTIVEPASSTHRUHANDLER_HPP

#include "codingUtilities/traits.hpp"
#include "codingUtilities/RTTypes.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Struct to help in dynamic dispatch based on type of contitutive model.
 * @tparam TYPES pack of derived constitutive types to handle
 */
template< typename ... TYPES >
struct ConstitutivePassThruHandler {};

/**
 * @brief Specialization for an empty type pack (always errors)
 */
template<>
struct ConstitutivePassThruHandler<>
{
  template< typename BASE, typename LAMBDA >
  static void execute( BASE & relation, LAMBDA lambda )
  {
    GEOS_UNUSED_VAR( relation, lambda );
    GEOS_ERROR( "The constitutive model " << relation.getDataContext() << " was not dispatched. " <<
                "The model type does not match the list of supported types." );
  }
};

/**
 * @brief Specialization for a non-empty type pack
 * @tparam TYPE first type to handle
 * @tparam TYPES rest of the type list
 */
template< typename TYPE, typename ... TYPES >
struct ConstitutivePassThruHandler< TYPE, TYPES... >
{
  template< typename BASE, typename LAMBDA >
  static void execute( BASE & relation, LAMBDA && lambda )
  {
    using Derived = add_const_if_t< TYPE, std::is_const< BASE >::value >;

    if( dynamicCast< Derived * >( &relation ) )
    {
      lambda( static_cast< Derived & >( relation ) );
    }
    else
    {
      ConstitutivePassThruHandler< TYPES... >::execute( relation, std::forward< LAMBDA >( lambda ) );
    }
  }
};

}//namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVEPASSTHRUHANDLER_HPP
