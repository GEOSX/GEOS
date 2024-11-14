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
 * @file FiniteElementSpace.cpp
 */

// TODO make this not dependent on this header...need better key implementation
#include "FiniteElementDiscretization.hpp"

#include "common/GeosxMacros.hpp"

namespace geos
{
using namespace dataRepository;
using namespace finiteElement;


FiniteElementDiscretization::FiniteElementDiscretization( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::orderString(), &m_order ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The order of the finite element basis." );

  registerWrapper( viewKeyStruct::formulationString(), &m_formulation ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( Formulation::Default ).
    setDescription( "Specifier to indicate any specialized formuations. "
                    "For instance, one of the many enhanced assumed strain "
                    "methods of the Hexahedron parent shape would be indicated "
                    "here" );

  registerWrapper( viewKeyStruct::useVemString(), &m_useVem ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Specifier to indicate whether to force the use of VEM" );
}

FiniteElementDiscretization::~FiniteElementDiscretization()
{}


void FiniteElementDiscretization::postInputInitialization()
{
  GEOS_ERROR_IF( m_useVem < 0 || m_useVem > 1,
                 getDataContext() << ": The flag useVirtualElements can be either 0 or 1" );
}

std::unique_ptr< FiniteElementBase >
FiniteElementDiscretization::factory( ElementType const parentElementShape ) const
{
  if( m_order==1 )
  {
    switch( parentElementShape )
    {
      case ElementType::Triangle:      return std::make_unique< H1_TriangleFace_Lagrange1_Gauss1 >();
      case ElementType::Quadrilateral: return std::make_unique< H1_QuadrilateralFace_Lagrange1_GaussLegendre2 >();
      // On polyhedra where FEM are available, we use VEM only if useVirtualElements is set to 1 in
      // the input file.
      case ElementType::Tetrahedron:
      {
        if( m_useVem == 1 )
        {
          return std::make_unique< H1_Tetrahedron_VEM_Gauss1 >();
        }
        else if( m_formulation == Formulation::DG )
        {
          return std::make_unique< BB1_Tetrahedron >();
        }
        else
        {
          return std::make_unique< H1_Tetrahedron_Lagrange1_Gauss1 >();
        }
      }
      case ElementType::Pyramid:
      {
        if( m_useVem == 1 )
        {
          return std::make_unique< H1_Pyramid_VEM_Gauss1 >();
        }
        else
        {
          return std::make_unique< H1_Pyramid_Lagrange1_Gauss5 >();
        }
      }
      case ElementType::Wedge:
      {
        if( m_useVem == 1 )
        {
#if !defined( GEOS_USE_HIP )
          return std::make_unique< H1_Wedge_VEM_Gauss1 >();
#else
          GEOS_ERROR( "Cannot compile this with HIP active." );
          return nullptr;
#endif
        }
        else
        {
          return std::make_unique< H1_Wedge_Lagrange1_Gauss6 >();
        }
      }
      case ElementType::Hexahedron:
      {
        if( m_useVem == 1 )
        {
#if !defined( GEOS_USE_HIP )
          return std::make_unique< H1_Hexahedron_VEM_Gauss1 >();
#else
          GEOS_ERROR( "Cannot compile this with HIP active." );
          return nullptr;
#endif
        }
        else if( m_formulation == Formulation::SEM )
        {
#if !defined( GEOS_USE_HIP )
          return std::make_unique< Q1_Hexahedron_Lagrange_GaussLobatto >();
#else
          GEOS_ERROR( "Cannot compile this with HIP active." );
          return nullptr;
#endif
        }
        else
        {
          return std::make_unique< H1_Hexahedron_Lagrange1_GaussLegendre2 >();
        }
      }
      // On more general polyhedra, we always use VEM
      case ElementType::Prism5:
      {
        return std::make_unique< H1_Prism5_VEM_Gauss1 >();
      }
      case ElementType::Prism6:
      {
        return std::make_unique< H1_Prism6_VEM_Gauss1 >();
      }
      case ElementType::Prism7:
      {
        return std::make_unique< H1_Prism7_VEM_Gauss1 >();
      }
      case ElementType::Prism8:
      {
        return std::make_unique< H1_Prism8_VEM_Gauss1 >();
      }
      case ElementType::Prism9:
      {
        return std::make_unique< H1_Prism9_VEM_Gauss1 >();
      }
      case ElementType::Prism10:
      {
        return std::make_unique< H1_Prism10_VEM_Gauss1 >();
      }
#if !defined( GEOS_USE_HIP )
      case ElementType::Prism11:
      {
        return std::make_unique< H1_Prism11_VEM_Gauss1 >();
      }
#endif
      default:
      {
        GEOS_ERROR( getDataContext() << ": Element type " << parentElementShape << " does not have an associated element formulation." );
      }
    }
    return {};
  }

  if( m_order==2 )
  {
    switch( parentElementShape )
    {
#if !defined( GEOS_USE_HIP )
      case ElementType::Hexahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::SEM,
                       getDataContext() << ": Element type Hexahedron with order 2 available" <<
                       " only when using the Spectral Element Method" );
        return std::make_unique< Q2_Hexahedron_Lagrange_GaussLobatto >();
#else
      GEOS_ERROR( "Cannot compile this with HIP active." );
#endif
      case ElementType::Tetrahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::DG,
                       getDataContext() << ": Element type Tetrahedron with order 2 available" <<
                       " only when using the Discontinuous Galerkin Method" );
        return std::make_unique< BB2_Tetrahedron >();
      default:
      {
        GEOS_ERROR( getDataContext() << ": Element type " << parentElementShape << " does not" <<
                    " have an associated element formulation." );
      }
    }
    return {};
  }

  if( m_order==3 )
  {
    switch( parentElementShape )
    {
#if !defined( GEOS_USE_HIP )
      case ElementType::Hexahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::SEM,
                       getDataContext() << ": Element type Hexahedron with order 3 available" <<
                       " only when using the Spectral Element Method" );
        return std::make_unique< Q3_Hexahedron_Lagrange_GaussLobatto >();
#else
      GEOS_ERROR( "Cannot compile this with HIP active." );
#endif
      case ElementType::Tetrahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::DG,
                       getDataContext() << ": Element type Tetrahedron with order 3 available" <<
                       " only when using the Discontinuous Galerkin Method" );
        return std::make_unique< BB3_Tetrahedron >();
      default:
      {
        GEOS_ERROR( getDataContext() << ": Element type " << parentElementShape << " does not" <<
                    " have an associated element formulation." );
      }
    }
    return {};
  }

  if( m_order==4 )
  {
    switch( parentElementShape )
    {
#if !defined( GEOS_USE_HIP )
      case ElementType::Hexahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::SEM,
                       getDataContext() << ": Element type Hexahedron with order 4 available only" <<
                       " when using the Spectral Element Method" );
        return std::make_unique< Q4_Hexahedron_Lagrange_GaussLobatto >();
#else
      GEOS_ERROR( "Cannot compile this with HIP active." );
#endif
      case ElementType::Tetrahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::DG,
                       getDataContext() << ": Element type Tetrahedron with order 4 available" <<
                       " only when using the Discontinuous Galerkin Method" );
        return std::make_unique< BB4_Tetrahedron >();
      default:
      {
        GEOS_ERROR( getDataContext() << ": Element type " << parentElementShape << " does not have" <<
                    " an associated element formulation." );
      }
    }
    return {};
  }

  if( m_order==5 )
  {
    switch( parentElementShape )
    {
#if !defined( GEOS_USE_HIP )
      case ElementType::Hexahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::SEM,
                       getDataContext() << ": Element type Hexahedron with order 5 available only" <<
                       " when using the Spectral Element Method" );
        return std::make_unique< Q5_Hexahedron_Lagrange_GaussLobatto >();
#else
      GEOS_ERROR( "Cannot compile this with HIP active." );
#endif
      case ElementType::Tetrahedron:
        GEOS_ERROR_IF( m_formulation != Formulation::DG,
                       getDataContext() << ": Element type Tetrahedron with order 5 available" <<
                       " only when using the Discontinuous Galerkin Method" );
        return std::make_unique< BB5_Tetrahedron >();
      default:
      {
        GEOS_ERROR( getDataContext() << ": Element type " << parentElementShape << " does not have" <<
                    " an associated element formulation." );
      }
    }
    return {};
  }
  GEOS_ERROR( getDataContext() << ": Element type " << parentElementShape << " does not have an" <<
              " associated element formulation." );
  return {};
}

REGISTER_CATALOG_ENTRY( Group, FiniteElementDiscretization, string const &, Group * const )

} /* namespace geos */
