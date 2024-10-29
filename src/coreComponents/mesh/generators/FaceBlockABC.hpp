/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_FACEBLOCKABC_HPP
#define GEOS_FACEBLOCKABC_HPP

#include "CellBlockUtilities.hpp"
#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @brief Block of 2d elements (geometrical surfaces in 3d).
 *
 * @details The @p FaceBlockABC represents a zone of 2d (@e i.e. surfacic) elements.
 * Unlike its volumic equivalent @p CellBlockABC, @p FaceBlockABC is @e not @e homogeneous
 * and may hold 2d elements of multiple types (triangles, quadrangles...).
 * @details In this class, we'll use the term @e 2d @e element for the elements of the @p FaceBlockABC,
 * which are geometrical surfaces (in 3d).
 * In the same way, we'll use the wording @e 2d @e face
 * to refer to the 2d boundaries of the @e 2d @e elements.
 * The @e 2d @e face are geometrical segments (in 3d).
 */
class FaceBlockABC : public dataRepository::Group
{
public:
  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  FaceBlockABC( string const & name,
                Group * const parent ):
    Group( name, parent )
  { }

  /**
   * @brief Get the number of 2d elements (geometrical surfaces in 3d).
   * @return Number of 2d elements in the current face block.
   */
  virtual localIndex num2dElements() const = 0;

  /**
   * @brief Get the number of 2d faces (geometrical segments in 3d).
   * @return Number of 2d faces in the current face block.
   */
  virtual localIndex num2dFaces() const = 0;

  /**
   * @brief Get the nodes of each 2d element (geometrical surfaces in 3d).
   * @return The mapping of first dimension @p num2dElements.
   * Second dimension depends on the number of nodes of each 2d element.
   * Element numbering is local to the FaceBlockABC. Node numbering is local to the rank.
   *
   * @details Mesh is supposed to be conformal, so each 2d element touches two 3d faces (or only one on boundaries).
   * The returned mapping provides @e all the nodes of all those 3d faces, in no particular order.
   * For example, a triangular 2d element will provide 6 nodes.
   */
  virtual ArrayOfArrays< localIndex > get2dElemToNodes() const = 0;

  /**
   * @brief Get the 3d edges of each 2d element (geometrical surfaces in 3d).
   * @return The mapping of first dimension @p num2dElements.
   * Second dimension depends on the number of edges of each 2d element.
   * Element numbering is local to the @p FaceBlockABC. Edges numbering is local to the rank.
   *
   * @details On faults and fractures, because nodes are duplicated (and collocated),
   * there are two collocated 3d edges as well.
   * But the returned mapping only refers to one unique 3d edge; the other one is simply not provided.
   *
   * @note One unique 3d edge gets chosen and the other one is simply ignored.
   * And this choice is consistent with the @p get2dFaceToEdge mapping.
   * That is, when one 3d edge gets selected over its collocated twin,
   * it's guaranteed that the same 3d edge will always be selected over the ignored one,
   * in all the mappings provided by this @p FaceBlockABC.
   * There is no particular rule to define which of the two 3d edges is going to be selected,
   * but this rule is consistently used throughout the different mapping computations.
   */
  virtual ArrayOfArrays< localIndex > get2dElemToEdges() const = 0;

  /**
   * @brief Get the 3d faces (of the volumic mesh) that are aside each 2d element (geometrical surfaces in 3d) of the @p FaceBlockABC.
   * @return A mapping of size @p num2dFaces * 2 (each 2d element matches two 3d faces).
   * Element numbering is local to the @p FaceBlockABC. Faces numbering is local to the rank.
   *
   * @details Mesh is supposed to be conformal, so each 2d element touches two 3d faces (or only one on boundaries).
   * In the case where the 2d element lies on a boundary, the second mapped value will be @p -1.
   *
   * @note Both mappings @p get2dElemToFaces and @p get2dElemToElems are consistent.
   * For the same given 2d element, the 3d face at index 0 (or 1) in the @p get2dElemToFaces mapping
   * will be part of the boundary of the 3d element at the same index 0 (or 1) in the @p get2dElemToElems mapping.
   */
  virtual ArrayOfArrays< localIndex > get2dElemToFaces() const = 0;

  /**
   * @brief Get the 3d elements that are aside each 2d element (geometrical surfaces in 3d) of the @p FaceBlockABC.
   * @return The mapping of dimension @p num2dElements * 2.
   * 2d element numbering is local to the @p FaceBlockABC.
   * 3d element numbering is local to the @p CellBlockABC.
   *
   * @details Mesh is supposed to be conformal, so each 2d element touches two 3d elements (or only one on boundaries).
   * In the case where the 2d element lies on a boundary, the second mapped value will be @p -1.
   *
   * @note Both mappings @p get2dElemToFaces and @p get2dElemToElems are consistent.
   * For the same given 2d element, the 3d face at index 0 (or 1) in the @p get2dElemToFaces mapping
   * will be part of the boundary of the 3d element at the same index 0 (or 1) in the @p get2dElemToElems mapping.
   */
  virtual ToCellRelation< ArrayOfArrays< localIndex > > get2dElemToElems() const = 0;

  /**
   * @brief Returns the collocated nodes for each node of each 2d element of the @p FaceBlockABC.
   * @return The bucket of collocated nodes.
   * Indices of the 2d elements (first dimension) local to the @p FaceBlockABC.
   * The size of the first dimension is equal to @p num2dElements.
   * The size of the second dimension is the number of nodes in the 2d element (e.g. 3 for a triangle).
   *
   * @details Each node of the @p FaceBlockABC is pointing to other nodes which are collocated.
   * Those other nodes are meant to be nodes of neighboring 3d cells.
   * All the collocated nodes of each node of each 2d element of the @p FaceBlockABC are gathered in the same bucket.
   * @warning There is no guarantee that the nodes for each 2d element are provided any order.
   * As well, there is no guarantee that buckets of collocated nodes are provided in any order.
   */
  virtual ArrayOfArrays< array1d< globalIndex > > get2dElemsToCollocatedNodesBuckets() const = 0;

  /**
   * @brief Get @e one 3d edge equivalent for each 2d faces (geometrical edges in 3d).
   * @return The mapping of size @p num2dFaces.
   * 2d face numbering is local to the FaceBlockABC.
   * 3d edge numbering is local to the rank.
   *
   * @details The current @p FaceBlockABC has its own 2d faces (which are segments in terms of 3d geometry).
   * Those segments are edges in the volumic mesh.
   * The returned mapping tells which 3d edge is the geometrical equivalent of the 2d face of the @p FaceBlockABC.
   * There can be two collocated 3d edges (because for faults and fractures situations, nodes are duplicated).
   * But the returned mapping only refers to one unique 3d edge; the other is simply not provided.
   *
   * @note Instead of returning only @e one 3d edge, we could actually return 2.
   * Note that both @p get2dElemToEdges and @p get2dFaceToEdge are consistent.
   * When one 3d edge gets selected over its collocated twin,
   * it's guaranteed that the same 3d edge will always be selected and the same always ignored,
   * in all the mappings provided by the @p FaceBlockABC.
   * There is no particular rule to define which of the two edges is going to be selected,
   * but this rule is consistently used throughout the different mapping computations.
   */
  virtual array1d< localIndex > get2dFaceToEdge() const = 0;

  /**
   * @brief Get the 2d element(s) (geometrical surfaces in 3d) connected to each 2d face (geometrical 3d segment).
   * @return The mapping of first dimension @p num2dFaces.
   * 2d face and 2d element numberings are both local to the @p FaceBlockABC.
   */
  virtual ArrayOfArrays< localIndex > get2dFaceTo2dElems() const = 0;

  /**
   * @brief Get local to global map for the 2d elements.
   * @return The mapping relationship as an array.
   */
  virtual array1d< globalIndex > localToGlobalMap() const = 0;
};

}

#endif // include guard
