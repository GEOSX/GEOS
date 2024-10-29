/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_MESH_PARTITIONDESCRIPTORABC_HPP
#define GEOS_MESH_PARTITIONDESCRIPTORABC_HPP

#include <set>
#include <array>

#include "common/DataTypes.hpp"


namespace geos
{

class PartitionDescriptorABC
{
public:
  ~PartitionDescriptorABC() = default;

  /**
   * @brief Gets a reference to the list of metis neighbor list.
   * @return A reference to the Metis neighbor list.
   */
  [[nodiscard]] virtual std::set< int > const & getMetisNeighborList() const = 0;

  [[nodiscard]] virtual array1d< int > getPartitions() const = 0;

  [[nodiscard]] virtual array1d< int > getPeriodic() const = 0;

  [[nodiscard]] virtual array1d< int > getCoords() const = 0;

  [[nodiscard]] virtual std::array< real64, 9 > getGrid() const = 0;

  [[nodiscard]] virtual std::array< real64, 3 > getBlockSize() const = 0;

  [[nodiscard]] virtual std::array< real64, 6 > getBoundingBox() const = 0;

};

} // end of namespace

#endif // Include guard
