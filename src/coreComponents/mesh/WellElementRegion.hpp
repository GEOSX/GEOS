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
 * @file WellElementRegion.hpp
 *
 */

#ifndef GEOS_MESH_WELLELEMENTREGION_HPP_
#define GEOS_MESH_WELLELEMENTREGION_HPP_

#include "mesh/ElementRegionBase.hpp"
#include "mesh/generators/LineBlockABC.hpp"

namespace geos
{

class MeshLevel;

/**
 * @class WellElementRegion
 * @brief This class specializes the element region for the case
 *        of a well. This class is also in charge of starting the
 *        construction of the well data structure in GenerateWell.
 */
class WellElementRegion : public ElementRegionBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  WellElementRegion( string const & name, Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~WellElementRegion() override;

  /**
   * @brief Deleted default constructor.
   */
  WellElementRegion() = delete;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName()
  { return "WellElementRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return catalogName(); }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Set the name of the InternalWellGenerator object of this well.
   * @param[in] name the name of the InternalWellGenerator object
   */
  void setWellGeneratorName( string const & name ) { m_wellGeneratorName = name; }

  /**
   * @brief Get the name of the InternalWellGenerator object of this well.
   * @return the name of the InternalWellGenerator object
   */
  string const & getWellGeneratorName() const { return m_wellGeneratorName; }

  /**
   * @brief Set the name of the WellControls object of this well.
   * @param name the name of the WellControls object
   */
  void setWellControlsName( string const & name ) { m_wellControlsName = name; }

  /**
   * @brief Get the name of the subRegion.
   * @return the name of the subRegion object
   */
  string const & getSubRegionName() const { return m_subRegionName; }

  ///@}

  /**
   * @name Construction of the well connectivity
   */
  ///@{

  /**
   * @brief Not implemented, this task is performed in GenerateWell.
   */
  void generateMesh( Group const & ) override {}

  /**
   * @brief Build the local well elements and perforations from global well geometry.
   * @param[in] mesh the mesh object (single level only)
   * @param[in] lineBlock the LineBlockABC containing the global well topology
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */
  void generateWell( MeshLevel & mesh,
                     LineBlockABC const & lineBlock,
                     globalIndex nodeOffsetGlobal,
                     globalIndex elemOffsetGlobal );

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return String key for the well control name
    static constexpr char const * wellControlsString() { return "wellControlsName"; }
    /// @return String key for the well generator name
    static constexpr char const * wellGeneratorString() { return "wellGeneratorName"; }

    /// ViewKey for the well control name
    dataRepository::ViewKey wellControlsName  = { wellControlsString() };
    /// ViewKey for the well generator name
    dataRepository::ViewKey wellGeneratorName = { wellGeneratorString() };
  }
  /// ViewKey struct for the WellElementRegion class
  viewKeysWellElementRegion;

private:

  /// Name of the (unique) subregion
  const string m_subRegionName;

  /// Name of the well constraint
  string m_wellControlsName;

  /// Name of the well generator
  string m_wellGeneratorName;

};

} /* namespace geos */

#endif /* GEOS_MESH_WELLELEMENTREGION_HPP_ */
