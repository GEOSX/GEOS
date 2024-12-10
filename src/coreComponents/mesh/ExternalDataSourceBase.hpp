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
 * @file ExternalDataSourceBase.hpp
 */

#ifndef GEOS_MESH_EXTERNALDATASOURCEBASE_HPP
#define GEOS_MESH_EXTERNALDATASOURCEBASE_HPP

#include "dataRepository/Group.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

/**
 *  @class ExternalDataSourceBase
 *  @brief The ExternalDataSourceBase class provides an abstract base class implementation for different mesh types.
 *	   The ExternalDataSourceBase is the Group specialization for different type of mesh handling.
 */
class ExternalDataSourceBase : public dataRepository::Group
{
public:

  /**
   * @brief Main constructor for ExternalDataSourceBase base class.
   * @param[in] name of the ExternalDataSourceBase object
   * @param[in] parent the parent Group pointer for the ExternalDataSourceBase object
   */
  explicit ExternalDataSourceBase( string const & name,
                                   Group * const parent );

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /// using alias for templated Catalog ExternalDataSourceBase type
  using CatalogInterface = dataRepository::CatalogInterface< ExternalDataSourceBase, string const &, Group * const >;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Accessor for the singleton Catalog object
   * @return a static reference to the Catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief This function provides the capability to open an external data repository
   * from another component whatever its format.
   */
  virtual void open() = 0;
};

}

#endif /* GEOS_MESH_EXTERNALDATASOURCEBASE_HPP */
