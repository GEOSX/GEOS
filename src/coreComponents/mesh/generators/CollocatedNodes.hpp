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

#ifndef GEOS_COLLOCATEDNODES_HPP
#define GEOS_COLLOCATEDNODES_HPP

#include "common/DataTypes.hpp"

#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>

namespace geos::vtk
{

/**
 * @brief Convenience wrapper around the raw vtk information.
 */
class CollocatedNodes
{
public:
  /**
   * @brief Build a convenience wrapper around the raw vtk collocated nodes information.
   * @param faceBlockName The face block name.
   * @param faceMesh The face mesh for which the collocated nodes structure will be fed.
   * @param isParallel Even if the global simulation is parallel,
   * this structure can be built on one unique rank or on multiple ranks.
   * When built on multiple ranks for the same face mesh, some additional extra checks are performed.
   * Those checks are useless when the data is manipulated on one unique rank.
   */
  CollocatedNodes( string const & faceBlockName,
                   vtkSmartPointer< vtkDataSet > faceMesh,
                   bool isParallel = true );

  /**
   * @brief For node @p i of the face block, returns all the duplicated global node indices in the main 3d mesh.
   * @param i the node in the face block (numbering is local to the face block).
   * @return The list of global node indices in the main 3d mesh.
   */
  std::vector< vtkIdType > const & operator[]( std::size_t i ) const
  {
    return m_collocatedNodes[i];
  }

  /**
   * @brief Number of duplicated nodes buckets.
   * Multiple nodes that are considered to be duplicated one of each other make one bucket.
   * @return The number of duplicated nodes buckets.
   */
  [[nodiscard]] std::size_t size() const
  {
    return m_collocatedNodes.size();
  }

private:

  /**
   * @brief Populates the mapping.
   * @param collocatedNodes The pointer to the vtk data array. Guaranteed not to be null.
   * @details Converts the raw vtk data to STL containers.
   */
  void init( vtkIdTypeArray const * collocatedNodes );

  /// For each node of the face block, lists all the collocated nodes in the main 3d mesh.
  std::vector< std::vector< vtkIdType > > m_collocatedNodes;
};

} // geos::vtk

#endif //GEOS_COLLOCATEDNODES_HPP
