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
 * @file FiniteElementDispatch.hpp
 */

#ifndef GEOS_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_
#define GEOS_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_

#include "common/GeosxMacros.hpp"
#include "elementFormulations/ConformingVirtualElementOrder1.hpp"
#include "elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp"
#include "elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp"
#include "elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp"
#if !defined( GEOS_USE_HIP )
#include "elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif
#include "elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"
#include "elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp"
#include "LvArray/src/system.hpp"

#define FE_1_TYPES \
  finiteElement::H1_Hexahedron_Lagrange1_GaussLegendre2, \
  finiteElement::H1_Wedge_Lagrange1_Gauss6, \
  finiteElement::H1_Tetrahedron_Lagrange1_Gauss1, \
  finiteElement::H1_Pyramid_Lagrange1_Gauss5

#define GL_FE_TYPES \
  finiteElement::Q1_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q2_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q3_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q4_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q5_Hexahedron_Lagrange_GaussLobatto

#if defined( GEOS_DISPATCH_VEM )

#define VEM_1_TYPES \
  finiteElement::H1_Tetrahedron_VEM_Gauss1, \
  finiteElement::H1_Prism5_VEM_Gauss1, \
  finiteElement::H1_Prism6_VEM_Gauss1, \
  finiteElement::H1_Prism7_VEM_Gauss1, \
  finiteElement::H1_Prism8_VEM_Gauss1, \
  finiteElement::H1_Prism9_VEM_Gauss1, \
  finiteElement::H1_Prism10_VEM_Gauss1

// can only compile these when not using cce+rocm
#define VEM_2_TYPES \
  finiteElement::H1_Hexahedron_VEM_Gauss1, \
  finiteElement::H1_Wedge_VEM_Gauss1, \
  finiteElement::H1_Prism11_VEM_Gauss1

#if !defined( GEOS_USE_HIP )
#define VEM_TYPES VEM_1_TYPES, VEM_2_TYPES
#else
#define VEM_TYPES VEM_1_TYPES
#endif

#define BASE_FE_TYPES FE_1_TYPES, VEM_TYPES

#else

#define BASE_FE_TYPES FE_1_TYPES

#endif

#if !defined( GEOS_USE_HIP )
// can only compile GL_FE_TYPES when not using cce+rocm
#define ALL_FE_TYPES BASE_FE_TYPES, GL_FE_TYPES
#else
#define ALL_FE_TYPES BASE_FE_TYPES
#endif



#define FE_TYPES_2D \
  finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre2, \
  finiteElement::H1_TriangleFace_Lagrange1_Gauss1

#define BASE_FE_TYPES_2D FE_TYPES_2D

namespace geos
{
namespace finiteElement
{

/**
 * @brief Helper structure for dynamic finite element dispatch
 * @tparam FE_TYPES list of finite element types to handle
 */
template< typename ... FE_TYPES >
struct FiniteElementDispatchHandler {};

/**
 * @brief Finite elemnt dispatch fall-through in case all casts were unsuccessful: always error
 */
template<>
struct FiniteElementDispatchHandler<>
{
  template< typename LAMBDA >
  static void
  dispatch3D( FiniteElementBase const & input,
              LAMBDA && GEOS_UNUSED_PARAM( lambda ) )
  {
    GEOS_ERROR( "finiteElement::dispatch3D() is not implemented for input of "<<typeid(input).name() );
    GEOS_UNUSED_VAR( input );
  }

  template< typename LAMBDA >
  static void
  dispatch3D( FiniteElementBase & input,
              LAMBDA && GEOS_UNUSED_PARAM( lambda ) )
  {
    GEOS_ERROR( "finiteElement::dispatch3D() is not implemented for input of "<<typeid(input).name() );
    GEOS_UNUSED_VAR( input );
  }

  template< typename LAMBDA >
  static void
  dispatch2D( FiniteElementBase const & input,
              LAMBDA && GEOS_UNUSED_PARAM( lambda ) )
  {
    GEOS_ERROR( "finiteElement::dispatch2D() is not implemented for input of: "<<LvArray::system::demangleType( &input ) );
    GEOS_UNUSED_VAR( input );
  }
};

/**
 * @brief Structure for recursive finite element dispatch implementation
 * @tparam FE_TYPE first finite element type to handle
 * @tparam FE_TYPES following finite element types to handle
 */
template< typename FE_TYPE, typename ... FE_TYPES >
struct FiniteElementDispatchHandler< FE_TYPE, FE_TYPES... >
{
  template< typename LAMBDA >
  static void
  dispatch3D( FiniteElementBase const & input,
              LAMBDA && lambda )
  {
    if( auto const * const ptr = dynamic_cast< FE_TYPE const * >(&input) )
    {
      lambda( *ptr );
    }
    else
    {
      FiniteElementDispatchHandler< FE_TYPES... >::dispatch3D( input, lambda );
    }
  }

  template< typename LAMBDA >
  static void
  dispatch3D( FiniteElementBase & input,
              LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< FE_TYPE * >(&input) )
    {
      lambda( *ptr );
    }
    else
    {
      FiniteElementDispatchHandler< FE_TYPES... >::dispatch3D( input, lambda );
    }
  }

  template< typename LAMBDA >
  static void
  dispatch2D( FiniteElementBase const & input,
              LAMBDA && lambda )
  {
    if( auto const * const ptr = dynamic_cast< FE_TYPE const * >(&input) )
    {
      lambda( *ptr );
    }
    else
    {
      FiniteElementDispatchHandler< FE_TYPES... >::dispatch2D( input, lambda );
    }
  }
};



template< typename LAMBDA >
void
dispatchlowOrder3D( FiniteElementBase const & input,
                    LAMBDA && lambda )
{
  if( auto const * const ptr1 = dynamic_cast< H1_Hexahedron_Lagrange1_GaussLegendre2 const * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 = dynamic_cast< H1_Wedge_Lagrange1_Gauss6 const * >(&input) )
  {
    lambda( *ptr2 );
  }
  else if( auto const * const ptr3 = dynamic_cast< H1_Tetrahedron_Lagrange1_Gauss1 const * >(&input) )
  {
    lambda( *ptr3 );
  }
  else if( auto const * const ptr4 = dynamic_cast< H1_Pyramid_Lagrange1_Gauss5 const * >(&input) )
  {
    lambda( *ptr4 );
  }
  else
  {
    GEOS_ERROR( "finiteElement::dispatchlowOrder3D() is not implemented for input of "<<LvArray::system::demangleType( &input ) );
  }
}

} // namespace finiteElement

} // namespace geos



#endif /* GEOS_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_ */
