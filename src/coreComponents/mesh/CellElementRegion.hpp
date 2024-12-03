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
 * @file CellElementRegion.hpp
 */

#ifndef GEOS_MESH_CELLELEMENTREGION_HPP_
#define GEOS_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geos
{

/**
 * @class CellElementRegion
 *
 * The CellElementRegion class contains the functionality to support the concept of a CellElementRegion in the element
 * hierarchy. CellElementRegion derives from ElementRegionBase and has an entry in the ObjectManagerBase catalog.
 */
class CellElementRegion : public ElementRegionBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{


  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  CellElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  CellElementRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~CellElementRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static string catalogName()
  { return "CellElementRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return catalogName(); }

  ///@}

  /**
   * @name Generation of the cell element mesh region
   */
  ///@{


  /**
   * @return List of user-requested mesh cellBlocks names.
   * @note the list may be incomplete / illegible if CellElementRegionSelectorhas not been used on
   *       the instance.
   */
  arrayView1d< string const > getCellBlockNames() const
  { return m_cellBlockNames.toViewConst(); }

  /**
   * @brief Select a cellBlock by its name for generateMesh().
   * @param cellBlockName string containing the cell block region name.
   */
  void addCellBlockName( string const & cellBlockName )
  {
    m_cellBlockNames.emplace_back( cellBlockName );
  }

  /**
   * @brief Select cellBlocks by their names for generateMesh().
   * @param cellBlockNames array of string containing the cell block region names.
   */
  template< typename StringContainerType >
  void setCellBlockNames( StringContainerType const & cellBlockNames )
  {
    m_cellBlockNames.clear();
    for( auto const & name: cellBlockNames )
    {
      m_cellBlockNames.emplace_back( name );
    }
  }

  /**
   * @brief register every cellBlocks that is requested in the cellBlockNames list.
   * @note Assume that the cellBlockNames list is filled valid `cellBlocks` names
   *       (use CellElementRegionSelector to verify & fill this list after the user-requests).
   * @param cellBlocks Cell blocks from where the mesh is extracted.
   */
  virtual void generateMesh( Group const & cellBlocks ) override;

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return String key for the coarsening ratio
    static constexpr char const * coarseningRatioString() {return "coarseningRatio"; }

    /// @return String key for the user-requested mesh cellBlocks qualifiers: cellblock names, cellblock match patterns, attribute values.
    static constexpr char const * sourceCellBlockNamesString() {return "cellBlocks"; }
  };

private:

  /// @brief List of user-requested mesh cellBlocks qualifiers: cellblock names, cellblock match patterns, attribute values.
  string_array m_cellBlockNames;

  /// @brief Coarsening ratio
  real64 m_coarseningRatio;


  /**
   * @return all cell-block names entries from m_cellBlockAttributeValues,
   * m_cellBlockMatchPatterns and m_cellBlockNames.
   * @param cellBlocks the input mesh cell-block list
   */
  std::set< string > computeSelectedCellBlocks( std::set< string > const & cellBlocksNames ) const;

};

} /* namespace geos */

#endif /* GEOS_MESH_CELLELEMENTREGION_HPP_ */
