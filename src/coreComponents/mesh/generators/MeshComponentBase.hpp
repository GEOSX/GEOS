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

/*
 * @file MeshComponentBase.hpp
 *
 */

#ifndef GEOS_MESH_GENERATORS_MESHCOMPONENTBASE_HPP_
#define GEOS_MESH_GENERATORS_MESHCOMPONENTBASE_HPP_

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

/**
 * @class MeshComponentBase
 *
 * Abstract base class defining the information provided by any the well generator class.
 */
class MeshComponentBase : public dataRepository::Group
{
public:

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  MeshComponentBase( string const & name,
                     Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~MeshComponentBase();

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName() { return "MeshComponentBase"; }

  /**
   * @brief Type alias for catalog interface used by this class. See CatalogInterface.
   */
  using CatalogInterface = dataRepository::CatalogInterface< MeshComponentBase, string const &, Group * const >;

  /// @copydoc dataRepository::Group::getCatalog()
  static CatalogInterface::CatalogType & getCatalog();
};

}
#endif /* GEOS_MESH_GENERATORS_MESHCOMPONENTBASE_HPP_ */
