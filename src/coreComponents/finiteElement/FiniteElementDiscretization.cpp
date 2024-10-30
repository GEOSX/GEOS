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
  switch( m_formulation )
  {
    case Formulation::Default:
      return createDefaultlElement( parentElementShape );
    case Formulation::SEM:
      return createSpectralElement( parentElementShape );
    default:
      GEOS_ERROR( getDataContext() << ": Formulation " << m_formulation << " is not supported." );
  }
}

// TOOD fix broken refactoring

std::unique_ptr<FiniteElementBase> createSpectralElement(ElementType const type) const {
  switch( type )
  {
    case ElementType::Voxel 
    {
      // @TODO -> need to check if voxels are cubes
      if (order == 1) return std::make_unique< Q1_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 2) return std::make_unique< Q2_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 3) return std::make_unique< Q3_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 4) return std::make_unique< Q4_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 5) return std::make_unique< Q5_Hexahedron_Lagrange_GaussLobatto >();
      GEOS_ERROR( getDataContext() << ": Element Voxel does not have an associated SEM formulation for order " << order << "." );
    }
    case ElementType::Hexahedron:
    {
      if (order == 1) return std::make_unique< Q1_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 2) return std::make_unique< Q2_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 3) return std::make_unique< Q3_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 4) return std::make_unique< Q4_Hexahedron_Lagrange_GaussLobatto >();
      if (order == 5) return std::make_unique< Q5_Hexahedron_Lagrange_GaussLobatto >();
      GEOS_ERROR( getDataContext() << ": Element Hexahedron does not have an associated SEM formulation for order " << order << "." );
    }
    default:
      GEOS_ERROR( getDataContext() << ": Element type " << type << " does not have an associated SEM formulation." );
  }
  return {};
}

std::unique_ptr<FiniteElementBase> createDefaultlElement(ElementType const type) const {
  // On polyhedra where FEM are available, we use VEM only if useVirtualElements is set to 1 in
  // the input file. On more general polyhedra (Prism), we always use VEM
  if( m_useVem == 1 )
  {
    switch( type )
    {
      case ElementType::Triangle:      return std::make_unique< H1_TriangleFace_Lagrange1_Gauss1 >();
      case ElementType::Quadrilateral: return std::make_unique< H1_QuadrilateralFace_Lagrange1_GaussLegendre2 >();
      case ElementType::Tetrahedron:   return std::make_unique< H1_Tetrahedron_VEM_Gauss1 >();
      case ElementType::Pyramid:       return std::make_unique< H1_Pyramid_VEM_Gauss1 >();
      case ElementType::Wedge:         return std::make_unique< H1_Wedge_VEM_Gauss1 >();
      case ElementType::Voxel:         return std::make_unique< H1_Hexahedron_VEM_Gauss1 >();
      case ElementType::Hexahedron:    return std::make_unique< H1_Hexahedron_VEM_Gauss1 >();
      case ElementType::Prism5:        return std::make_unique< H1_Prism5_VEM_Gauss1 >();
      case ElementType::Prism6:        return std::make_unique< H1_Prism6_VEM_Gauss1 >();
      case ElementType::Prism7:        return std::make_unique< H1_Prism7_VEM_Gauss1 >();
      case ElementType::Prism8:        return std::make_unique< H1_Prism8_VEM_Gauss1 >();
      case ElementType::Prism9:        return std::make_unique< H1_Prism9_VEM_Gauss1 >();
      case ElementType::Prism10:       return std::make_unique< H1_Prism10_VEM_Gauss1 >();
      case ElementType::Prism11:       return std::make_unique< H1_Prism11_VEM_Gauss1 >();
      default:
        GEOS_ERROR( getDataContext() << ": Element type " << type << " does not have an associated default element formulation." );
    }
  }
  else
  {
    switch( type )
    {
      case ElementType::Triangle:      return std::make_unique< H1_TriangleFace_Lagrange1_Gauss1 >();
      case ElementType::Quadrilateral: return std::make_unique< H1_QuadrilateralFace_Lagrange1_GaussLegendre2 >();
      case ElementType::Tetrahedron:   return std::make_unique< H1_Tetrahedron_Lagrange1_Gauss1 >();
      case ElementType::Pyramid:       return std::make_unique< H1_Pyramid_Lagrange1_Gauss5 >();
      case ElementType::Wedge:         return std::make_unique< H1_Wedge_Lagrange1_Gauss6 >();
      case ElementType::Voxel:         return std::make_unique< H1_Hexahedron_Lagrange1_GaussLegendre2 >();
      case ElementType::Hexahedron:    return std::make_unique< H1_Hexahedron_Lagrange1_GaussLegendre2 >();
      case ElementType::Prism5:        return std::make_unique< H1_Prism5_VEM_Gauss1 >();
      case ElementType::Prism6:        return std::make_unique< H1_Prism6_VEM_Gauss1 >();
      case ElementType::Prism7:        return std::make_unique< H1_Prism7_VEM_Gauss1 >();
      case ElementType::Prism8:        return std::make_unique< H1_Prism8_VEM_Gauss1 >();
      case ElementType::Prism9:        return std::make_unique< H1_Prism9_VEM_Gauss1 >();
      case ElementType::Prism10:       return std::make_unique< H1_Prism10_VEM_Gauss1 >();
      case ElementType::Prism11:       return std::make_unique< H1_Prism11_VEM_Gauss1 >();
      default:
        GEOS_ERROR( getDataContext() << ": Element type " << type << " does not have an associated default element formulation." );
    }
  }
  return {};
}

REGISTER_CATALOG_ENTRY( Group, FiniteElementDiscretization, string const &, Group * const )

} /* namespace geos */
