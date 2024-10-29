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

#ifndef GEOS_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_
#define GEOS_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"



namespace geos
{

/**
 * @class ChomboCoupler
 * @brief A class managing data exchange with CHOMBO.
 */

class ChomboCoupler
{
public:

  /**
   * @brief Construct a new ChomboCoupler.
   * @param comm Communicator used in reading/writing from/to file.
   * @param outputPath filename The name of the file to write out to.
   * @param inputPath filename The name of the file to read from.
   * @param mesh The mesh to communicate.
   */
  ChomboCoupler( MPI_Comm const comm, const string & outputPath, const string & inputPath, MeshLevel & mesh );

  /**
   * @brief Write data to file.
   * @param dt the current time step.
   */
  void write( double dt );

  /**
   * @brief Read data from file.
   * @param usePressures If true, pressure are read in from file
   */
  void read( bool usePressures );

private:
  /**
   * @brief Copy nodal data into local arrays.
   * @details The purpose of this method is to de-permute the nodal data since the
   *   coupling format has a fixed data layout. If this copying becomes a bottleneck
   *   the coupling format should be updated so we can do a zero copy write.
   */
  void copyNodalData();

  /// The MPI communicator used to read and write the file.
  MPI_Comm const m_comm;
  /// The path to write the file to.
  string const m_outputPath;
  /// The path to read from.
  string const m_inputPath;
  /// The global face offset.
  std::int64_t m_face_offset;
  /// The local number of faces written.
  std::int64_t m_n_faces_written;
  /// The global node offset.
  std::int64_t m_node_offset;
  /// The local number of nodes written.
  std::int64_t m_n_nodes_written;
  /// The mesh to communicate.
  MeshLevel & m_mesh;
  /// Not sure why this is here.
  int m_counter;
  /// A copy of the nodal reference position.
  array2d< real64 > m_referencePositionCopy;
  /// A copy of the nodal displacement.
  array2d< real64 > m_displacementCopy;
  /// A copy of the nodal velocity.
  array2d< real64 > m_velocityCopy;
};

} /* namespace geos */

#endif /* GEOS_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_ */
