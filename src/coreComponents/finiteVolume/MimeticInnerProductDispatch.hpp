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
 * @file MimeticInnerProductDispatch.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTDISPATCH_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTDISPATCH_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiRTInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"
#include "LvArray/src/system.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @struct MimeticInnerProductTypeStrings
 * @brief Struct containing the keys to all the implemented inner products
 */
struct MimeticInnerProductTypeStrings
{
  /// string for the TPFA inner product
  static constexpr auto TPFA      = "TPFA";
  /// string for the quasi-TPFA inner product
  static constexpr auto QuasiTPFA = "quasiTPFA";
  /// string for the quasi-RT inner product
  static constexpr auto QuasiRT   = "quasiRT";
  /// string for the Simple inner product
  static constexpr auto Simple    = "simple";
  /// string for the inner product of Beirao da Veiga, Lipnikov, Manzini
  static constexpr auto BdVLM     = "beiraoDaVeigaLipnikovManzini";
};

/**
 * @brief Dispatch for the selection of the mimetic inner product
 * @tparam LAMBDA the type of the lambda
 * @param input the operator implementing the desired mimetic inner product
 * @param lambda the function that will launch the FluxKernel of the hybrid FVM solver
 */
template< typename LAMBDA >
void
mimeticInnerProductDispatch( MimeticInnerProductBase const & input,
                             LAMBDA && lambda )
{
  if( auto const * const ptr1 = dynamic_cast< TPFAInnerProduct const * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 = dynamic_cast< QuasiTPFAInnerProduct const * >(&input) )
  {
    lambda( *ptr2 );
  }
  else if( auto const * const ptr3 = dynamic_cast< QuasiRTInnerProduct const * >(&input) )
  {
    lambda( *ptr3 );
  }
  else if( auto const * const ptr4 = dynamic_cast< SimpleInnerProduct const * >(&input) )
  {
    lambda( *ptr4 );
  }
  else if( auto const * const ptr5 = dynamic_cast< BdVLMInnerProduct const * >(&input) )
  {
    lambda( *ptr5 );
  }
  else
  {
    GEOS_ERROR( "mimeticInnerProductDispatch() is not implemented for input of " << LvArray::system::demangleType( input ) );
  }
}

/**
 * @brief Dispatch for the selection of the mimetic inner product
 * @tparam LAMBDA the type of the lambda
 * @param input the operator implementing the desired mimetic inner product
 * @param lambda the function that will launch the FluxKernel of the hybrid FVM solver
 */
template< typename LAMBDA >
void
mimeticInnerProductDispatch( MimeticInnerProductBase & input,
                             LAMBDA && lambda )
{
  if( auto * const ptr1 = dynamic_cast< TPFAInnerProduct * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto * const ptr2 = dynamic_cast< QuasiTPFAInnerProduct * >(&input) )
  {
    lambda( *ptr2 );
  }
  else if( auto * const ptr3 = dynamic_cast< QuasiRTInnerProduct * >(&input) )
  {
    lambda( *ptr3 );
  }
  else if( auto * const ptr4 = dynamic_cast< SimpleInnerProduct * >(&input) )
  {
    lambda( *ptr4 );
  }
  else if( auto * const ptr5 = dynamic_cast< BdVLMInnerProduct * >(&input) )
  {
    lambda( *ptr5 );
  }
  else
  {
    GEOS_ERROR( "mimeticInnerProductDispatch() is not implemented for input of " << LvArray::system::demangleType( input ) );
  }
}

/**
 * @brief Dispatch for the selection of the mimetic inner product (limited number of possible templates).
 *        The purpose of this function is to reduce the number of possible templates to speed up the compilation
 *        of CompositionalMultiphaseHybridFVM
 * @tparam LAMBDA the type of the lambda
 * @param input the operator implementing the desired mimetic inner product
 * @param lambda the function that will launch the FluxKernel of the hybrid FVM solver
 */
template< typename LAMBDA >
void
mimeticInnerProductReducedDispatch( MimeticInnerProductBase const & input,
                                    LAMBDA && lambda )
{
  if( auto const * const ptr1 = dynamic_cast< TPFAInnerProduct const * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 = dynamic_cast< BdVLMInnerProduct const * >(&input) )
  {
    lambda( *ptr2 );
  }
  else
  {
    GEOS_ERROR( "mimeticInnerProductReducedDispatch() is not implemented for input of " << LvArray::system::demangleType( input ) );
  }
}

/**
 * @brief Dispatch for the selection of the mimetic inner product (limited number of possible templates).
 *        The purpose of this function is to reduce the number of possible templates to speed up the compilation
 *        of CompositionalMultiphaseHybridFVM.
 * @tparam LAMBDA the type of the lambda
 * @param input the operator implementing the desired mimetic inner product
 * @param lambda the function that will launch the FluxKernel of the hybrid FVM solver
 */
template< typename LAMBDA >
void
mimeticInnerProductReducedDispatch( MimeticInnerProductBase & input,
                                    LAMBDA && lambda )
{
  if( auto * const ptr1 = dynamic_cast< TPFAInnerProduct * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto * const ptr2 = dynamic_cast< BdVLMInnerProduct * >(&input) )
  {
    lambda( *ptr2 );
  }
  else
  {
    GEOS_ERROR( "mimeticInnerProductReducedDispatch() is not supported for input of " << LvArray::system::demangleType( input ) );
  }
}


} // end namespace mimeticInnerProduct

} // end namespace geos

#endif /* GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTDISPATCH_HPP_ */
