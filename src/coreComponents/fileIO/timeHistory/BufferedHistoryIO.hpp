/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_FILEIO_TIMEHISTORY_BUFFEREDHISTORYIO_HPP_
#define GEOS_FILEIO_TIMEHISTORY_BUFFEREDHISTORYIO_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @class BufferedHistoryIO
 * @brief An abstract class for performing buffered history output.
 */
class BufferedHistoryIO
{
public:

  /// Destructor
  virtual ~BufferedHistoryIO() {}

  /**
   * @brief Get the head of the internal history storage buffer.
   * @return The head of the internal history storage buffer to be written to.
   */
  virtual buffer_unit_type * getBufferHead() = 0;

  /**
   * @brief Perform and intialization needed for time-history output.
   * @param[in] existsOkay Whether it is acceptable for the intended output target to already exist
   *                       ( false on start from scratch, true on restart ).
   */
  virtual void init( bool existsOkay ) = 0;

  /**
   * @brief Write the buffered history data to the output target.
   */
  virtual void write() = 0;

  /**
   * @brief Ensure the repressentation of the data in the output target is dense and terse.
   * @note Typically the file will be oversized to receive data writes without constant resizing.
   */
  virtual void compressInFile() = 0;

  /**
   * @brief Update the number of items being stored for IO in this object.
   * @param count [in] The new number of items being collected
   */
  virtual void updateCollectingCount( localIndex count ) = 0;

  /**
   * @brief Query the number of history states currently stored in the internal buffer.
   * @return The number of discrete time history records buffered to be written.
   * @note Since the size of each discrete time history can change, this should not be used
   *       to calculate size, but is useful to check for consistency between collectors
   *       that should be operating at the same cadence (ie the time collector and the data collector).
   * @deprecated Also, note that this member function is related to restarting HDF5 files.
   *             Which means that getting this information popping through the abstract interface
   *             is a shortcut that should eventually be fixed by keeping this information in the HDF5 buffer.
   */
  virtual localIndex getBufferedCount() = 0;

  /**
   * @brief Get the log-level for BufferedHistoryIO classes
   * @return the current log-level
   */
  int getLogLevel() const { return m_logLevel; }

  /**
   * @brief Set the log-level for BufferedHistoryIO classes
   * @param[in] logLevel the log-level to set
   */
  void setLogLevel( int logLevel ) { m_logLevel = logLevel; }

private:

  /// the log-level
  int m_logLevel = 0;
};

}
#endif
