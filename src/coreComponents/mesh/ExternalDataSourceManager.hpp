/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ExternalDataSourceManager.hpp
 */

#ifndef GEOS_MESH_EXTERNALDATASOURCEMANAGER_HPP_
#define GEOS_MESH_EXTERNALDATASOURCEMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

/**
 * @class ExternalDataSourceManager
 * @brief This class manages a data repository whereof objects can be imported to GEOS (reservoir mesh, well mesh)
 */
class ExternalDataSourceManager : public dataRepository::Group
{
public:

  /**
   * @brief Constructor for the ExternalDataSourceManager object.
   * @param[in] name the name of the ExternalDataSourceManager object in the repository
   * @param[in] parent the parent group of the ExternalDataSourceManager object being constructed
   */
  ExternalDataSourceManager( string const & name,
                             Group * const parent );

  virtual ~ExternalDataSourceManager() override;


  /**
   * @brief Create a new sub data repository.
   * @param[in] childKey the key of the new object in the ObjectCatalog
   * @param[in] childName the name of the new object in the collection of sub-meshes
   * @return A pointer to the Group node in the dataRepository of the new object created
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /**
   * @brief Generate the meshes of the physical DomainPartition.
   * @param[in] domain a reference to the physical domain
   */
  void open( DomainPartition & domain );

private:

  /**
   * @brief Deleted default constructor of the ExternalDataSourceManager
   */
  ExternalDataSourceManager() = delete;

};

} /* namespace geos */

#endif /* GEOS_MESH_EXTERNALDATASOURCEMANAGER_HPP_ */
