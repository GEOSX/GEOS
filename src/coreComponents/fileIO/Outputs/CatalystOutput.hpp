/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CatalystOutput.hpp
 */

#ifndef GEOSX_FILEIO_OUTPUTS_CATALYSTOUTPUT_HPP_
#define GEOSX_FILEIO_OUTPUTS_CATALYSTOUTPUT_HPP_

#include "BlueprintOutput.hpp"

#include <memory> // for unique_ptr

namespace geos
{

/**
 * @class CatalystOutput
 * @brief A class for outputing to a catalyst based in-situ pipeline.
 */
class CatalystOutput : public BlueprintOutput
{
public:

  /**
   * @brief Construct a new CatalystOutput object.
   * @param name The name of the CatalystObject in the data repository.
   * @param parent The parent Group.
   */
  CatalystOutput( string const & name,
                   Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~CatalystOutput() override;

  /**
   * @brief Get the name used to register this object in an XML file.
   * @return The string "Catalyst".
   */
  static string catalogName() { return "Catalyst"; }

  /**
   * @brief Launches a catalyst execution
   * @copydetails EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Launches a last catalyst execution for the data present.
   * @copydetails ExecutableGroup::cleanup()
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;
  /**
   * @brief Return PyCatalystOutput type.
   * @return Return PyCatalystOutput type.
   */
#if defined(GEOSX_USE_PYGEOSX)
  virtual PyTypeObject * getPythonType() const override;
#endif

private:
  static const std::string getEnv( const char* );
  ///@{
  /**
   */
  struct CatalystInternals;
  std::unique_ptr<CatalystInternals> internal;
  ///@}
};

}

#endif // GEOSX_FILEIO_OUTPUTS_CATALYSTOUTPUT_HPP_
