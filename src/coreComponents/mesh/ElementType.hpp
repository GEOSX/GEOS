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
 * @file ElementType.hpp
 */

#ifndef GEOS_MESH_ELEMENTTYPE_HPP
#define GEOS_MESH_ELEMENTTYPE_HPP

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

/**
 * @brief Denotes type of cell/element shape
 */
enum class ElementType : integer
{
  Vertex,        ///< Single-node vertex element
  Line,          ///< Two-node line segment
  Triangle,      ///< Three-node triangle
  Quadrilateral, ///< Four-node quadrilateral
  Polygon,       ///< General polygonal element
  Tetrahedron,   ///< Four-node tetrahedral element
  Pyramid,       ///< Five-node pyramid element
  Wedge,         ///< Six-node wedge element
  Hexahedron,    ///< Eight-node hexahedral element
  Prism5,        ///< Ten-node pentagonal prism
  Prism6,        ///< Twelve-node hexagonal prism
  Prism7,        ///< Heptagonal prism
  Prism8,        ///< Octagonal prism
  Prism9,        ///< Nonagonal prism
  Prism10,       ///< Decagonal prism
  Prism11,       ///< Hendecagonal prism
  Polyhedron,    ///< General polyhedral element
  // NOTE: If you add anything below Polyhedron,
  // don't forget to update numElementTypes() below.
};

/**
 * @brief @return number of supported element types
 * @note this MUST be updated if a new element type is inserted after Polyhedron
 */
inline constexpr integer numElementTypes()
{
  return static_cast< integer >( ElementType::Polyhedron ) + 1;
}

/**
 * @brief Get number of spatial dimensions of element type
 * @param elementType type of element
 * @return number of spatial dimensions (1-3)
 */
inline int getElementDim( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Vertex:        return 0;
    case ElementType::Line:          return 1;
    case ElementType::Triangle:
    case ElementType::Quadrilateral:
    case ElementType::Polygon:       return 2;
    case ElementType::Tetrahedron:
    case ElementType::Pyramid:
    case ElementType::Wedge:
    case ElementType::Hexahedron:
    case ElementType::Prism5:
    case ElementType::Prism6:
    case ElementType::Prism7:
    case ElementType::Prism8:
    case ElementType::Prism9:
    case ElementType::Prism10:
    case ElementType::Prism11:
    case ElementType::Polyhedron:    return 3;
  }
  return 0;
}

/// Strings for ElementType
ENUM_STRINGS( ElementType,
              "Vertex",
              "BEAM",
              "C2D3",
              "C2D4",
              "Polygon",
              "C3D4",
              "C3D5",
              "C3D6",
              "C3D8",
              "PentagonalPrism",
              "HexagonalPrism",
              "HeptagonalPrism",
              "OctagonalPrism",
              "NonagonalPrism",
              "DecagonalPrism",
              "HendecagonalPrism",
              "Polyhedron" );

/**
 * @brief Returns a string describing the element.
 * @param[in] type The element type.
 * @return The name.
 * @warning This information will be visible in the input file... Consider refactoring with great care.
 */
inline string getElementTypeName( ElementType const type )
{
  switch( type )
  {
    case ElementType::Hexahedron:  return "hexahedra";
    case ElementType::Tetrahedron: return "tetrahedra";
    case ElementType::Wedge:       return "wedges";
    case ElementType::Pyramid:     return "pyramids";
    case ElementType::Prism5:      return "pentagonalPrisms";
    case ElementType::Prism6:      return "hexagonalPrisms";
    case ElementType::Prism7:      return "heptagonalPrisms";
    case ElementType::Prism8:      return "octagonalPrisms";
    case ElementType::Prism9:      return "nonagonalPrisms";
    case ElementType::Prism10:     return "decagonalPrisms";
    case ElementType::Prism11:     return "hendecagonalPrisms";
    case ElementType::Polyhedron:  return "polyhedra";
    default:
    {
      GEOS_ERROR( "Element type '" << type << "' is not supported" );
      return {};
    }
  }
}

/// String available for mesh errors
inline auto constexpr generalMeshErrorAdvice = "Consider checking the validity of your mesh with "
                                               "the `mesh_doctor` GEOS python tools (documentation at "
                                               "https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/pythonTools/mesh_doctor.html).";

} // namespace geos

#endif //GEOS_MESH_ELEMENTTYPE_HPP
