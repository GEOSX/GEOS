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
 * @file PythonOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_PYTHONOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_PYTHONOUTPUT_HPP_

#include "fileIO/Outputs/OutputBase.hpp"

namespace geos
{

/**
 * @class PythonOutput
 * @brief Performs no actual output but returns control to the Python interpreter.
 */
class PythonOutput : public OutputBase
{
public:

  /**
   * @brief Constructor.
   * @param name The name of this PythonOutput in the data repository.
   * @param parent The parent group.
   */
  PythonOutput( string const & name,
                Group * const parent ):
    OutputBase( name, parent )
  {}

  /**
   * @brief Destructor.
   */
  virtual ~PythonOutput() override
  {}

  /**
   * @brief Get the name used to register this object in an XML file.
   * @return The string "Python".
   */
  static string catalogName() { return "Python"; }

  /**
   * @brief Signals the EventManager to exit the loop early, and return to Python.
   * @note It is an error to use this event when not running through pygeosx.
   * @copydetails EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    GEOS_UNUSED_VAR( time_n );
    GEOS_UNUSED_VAR( dt );
    GEOS_UNUSED_VAR( cycleNumber );
    GEOS_UNUSED_VAR( eventCounter );
    GEOS_UNUSED_VAR( eventProgress );
    GEOS_UNUSED_VAR( domain );
    return true;
  }
};


} // namespace geos

#endif // GEOS_FILEIO_OUTPUTS_PYTHONOUTPUT_HPP_
