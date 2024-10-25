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

#ifndef GEOS_HDFFILE_HPP
#define GEOS_HDFFILE_HPP

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @class HDFFile
 * A class used to control access to an HDF file target.
 */
class HDFFile
{
public:
  /**
   * @brief Constructor -- this creates/opens the target file for read/write.
   * @param[in] fnm The filename.
   * @param[in] deleteExisting Whether to remove/recreate if a file with the same name exists.
   * @param[in] parallelAccess Whether to access one file in parallel or one file per rank in the comm.
   * @param[in] comm An MPI communicator where each rank in the communicator will be accessing the target file.
   */
  HDFFile( string const & fnm, bool deleteExisting, bool parallelAccess, MPI_Comm comm );

  /**
   * @brief Closes the file and accessors.
   */
  ~HDFFile();

  /**
   * @brief Whether a dataset/group with the specified name exists in the target.
   * @param[in] name The dataset/group name to check for.
   * @return Whether the dataset/group exists in the target.
   */
  bool hasDataset( const string & name ) const;

  /**
   * @brief Get the HDF hid_t file identifier.
   * @return the HDF hid_t file id.
   */
  operator int64_t() const { return m_fileId; }
private:
  /// The filename
  string m_filename;
  /// The hdf file id
  int64_t m_fileId;
  /// The hdf file access properties list id
  int64_t m_faplId;
  /// Whether to open the same file in parallel, or open a file per process
  bool m_mpioFapl;
  /// The comminator to operate on the file collectively over
  MPI_Comm m_comm;
};

}

#endif //GEOS_HDFFILE_HPP
