
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
 * @file Region.hpp
 */

#ifndef GEOS_MESH_GENERATORS_REGION_HPP
#define GEOS_MESH_GENERATORS_REGION_HPP

#include "MeshComponentBase.hpp"

namespace geos
{

/**
 * @brief Region parameters with Group capabilities
 *
 * This class has dataRepository::Group capabilities to allow for XML input.
 *
 */
class Region : public MeshComponentBase
{
public:
  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  Region( const string & name, Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~Region() override;

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName() { return "Region"; }

  ///@cond DO_NOT_DOCUMENT
  /// Keys appearing in XML
  struct viewKeyStruct
  {
    static constexpr char const * idString() { return "id"; }
    static constexpr char const * pathInRepositoryString()  { return "pathInRepository"; }
  };
  /// @endcond

private:

  /// Interval region identifier
  integer m_id = 0;

  /// Path of the dataset in the repository
  string m_pathInRepository = "";
};

} // namespace GEOS

#endif
