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

/**
 * @file RestartOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_RESTARTOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_RESTARTOUTPUT_HPP_

#include "OutputBase.hpp"


namespace geos
{

/**
 * @class RestartOutput
 *
 * A class for creating restart-based outputs
 * Note: this class is redundant, and is expected to be retired
 */
class RestartOutput : public OutputBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group(string const & name, Group * const parent)
  RestartOutput( string const & name,
                 Group * const parent );

  /// Destructor
  virtual ~RestartOutput() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "Restart"; }

  /**
   * @brief Writes out a restart file.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Write one final restart file as the code exits
   * @copydetails ExecutableGroup::cleanup()
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain );
  }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    dataRepository::ViewKey writeFEMFaces = { "writeFEMFaces" };
  } viewKeys;
  /// @endcond
};


} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_RESTARTOUTPUT_HPP_ */
