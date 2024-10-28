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
 * @file CellElementRegionSelector.hpp
 */

#ifndef GEOS_MESH_CELLELEMENTREGIONSELECTOR_HPP_
#define GEOS_MESH_CELLELEMENTREGIONSELECTOR_HPP_

#include "mesh/CellElementRegion.hpp"

namespace geos
{

/**
 * @class CellElementRegionSelector
 *
 * CellElementRegionSelector allows a CellElementRegion to safely select cell-blocks according to the user input.
 * It handles all user input checks and throws an exception in the event of inconsistency.
 */
class CellElementRegionSelector
{
public:

  /**
   * @brief Construct a new CellElementRegionSelector.
   * @param cellBlocks a Group containing all the available cell-blocks.
   * @param cellBlocksRegion A map of the cellblocks name lists for each region attributes value.
   */
  CellElementRegionSelector( dataRepository::Group const & cellBlocks,
                             std::map< integer, std::set< string > > const & cellBlocksRegion );

  /**
   * @brief Select the mesh cell-blocks for the specified region following the user inputs.
   * @throw an InputError if the user setting are inconsistant.
   * @param region the region for which we want to select the cell-blocks.
   * @return the selected cell-blocks names.
   */
  std::set< string > buildCellBlocksSelection( CellElementRegion const & region );

  /**
   * @throw An InputError if region cell-blocks selections is inconsistent:
   *        - cell-block is in more than one region,
   *        - orphan cell-block
   * @todo For now, multiple regions per cell is not supported (ElementRegionManager::getCellBlockToSubRegionMap()).
   *       We could refactor the CellElementRegion & Mesh classes so regions are mapped to cell-blocks IN the mesh (and potentially
   *       to multiple regions per cell). So, for external meshes, the cell-blocks would no longer be exposed to the final user.
   */
  void checkSelectionConsistency() const;

private:

  /// @brief The separator between the regionAttribute value and the type of the shape in a given cellBlock.
  static constexpr string_view cellBlockTypeSeparator = "_";

  /// @brief A map that link every cell-block name to the CellElementRegion(s) that references it.
  std::map< string, std::vector< CellElementRegion const * > > m_cellBlocksOwners;

  /// @brief A map that link every region attribute values to the CellElementRegion(s) that references it.
  std::map< string, std::vector< CellElementRegion const * > > m_regionAttributesOwners;

  /**
   * @brief A map of the cellblocks name lists for each region attributes value. Internal attribute type
   *        is integer to facilitate comparison with cellBlock qualifiers.
   */
  std::map< string, std::set< string > const & > m_regionAttributesCellBlocks;

  /**
   * @return A set of the FNMatch pattern from the provided lists.
   * @param region The region for which we collect the match-patterns.
   * @param requestedAttributeValues The user requested attribute values. They will get converted to
   *                                 match-pattern that select the coresponding cell-blocks.
   * @param requestedMatchPatterns The match patterns that the user requested explicitely.
   * @throw An InputError if the attribute values does not exist in the mesh.
   */
  std::set< string > buildMatchPatterns( CellElementRegion const & region,
                                         std::set< string > const & attributeValues,
                                         std::set< string > const & matchPatterns ) const;

  /**
   * @return A set of the cell-blocks that the provided match-patterns select.
   * @param region The region for which we collect the cell-blocks.
   * @param matchPatterns The FNMatch pattern
   * @throw An InputError if a FNMatch pattern does not select any cell-block.
   */
  std::set< string > getMatchingCellblocks( CellElementRegion const & region,
                                            string_view matchPattern ) const;

  /**
   * @brief A set of the cell-blocks that the provided match-patterns select.
   * @param region The region for which we collect the cell-blocks.
   * @param matchPatterns The FNMatch pattern set.
   * @throw An InputError if a FNMatch pattern does not select any cell-block.
   */
  void verifyRequestedCellBlocks( CellElementRegion const & region,
                                  std::set< string > const & cellBlockNames ) const;

  /**
   * @brief Effectively select the specified cell-blocks & region attribute for the specified region.
   * @param region The region for which when want to select cell-blocks
   * @param attributeValues The attribute values we want to select (can be empty).
   * @param cellBlockNames The cell-block names we want to select (can be empty).
   */
  void registerRegionSelection( CellElementRegion const & region,
                                std::set< string > const & cellBlockNames,
                                std::set< string > const & attributeValues );

};

} /* namespace geos */

#endif /* GEOS_MESH_CELLELEMENTREGIONSELECTOR_HPP_ */
