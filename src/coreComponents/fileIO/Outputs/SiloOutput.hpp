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
 * @file SiloOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_SILOOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_SILOOUTPUT_HPP_

#include "OutputBase.hpp"


namespace geos
{

/**
 * @class SiloOutput
 *
 * A class for creating silo-based outputs
 */
class SiloOutput : public OutputBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group(string const & name, Group * const parent)
  SiloOutput( string const & name,
              Group * const parent );

  /// Destructor
  virtual ~SiloOutput() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "Silo"; }

  /**
   * @brief Writes out a Silo plot file.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Writes out a Silo plot file at the end of the simulation.
   * @copydoc ExecutableGroup::cleanup()
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
  struct viewKeysStruct : OutputBase::viewKeysStruct
  {
    static constexpr auto plotFileRoot = "plotFileRoot";
    static constexpr auto writeEdgeMesh = "writeEdgeMesh";
    static constexpr auto writeFaceMesh = "writeFEMFaces";
    static constexpr auto writeCellElementMesh = "writeCellElementMesh";
    static constexpr auto writeFaceElementMesh = "writeFaceElementMesh";
    static constexpr auto plotLevel = "plotLevel";
    static constexpr auto onlyPlotSpecifiedFieldNames = "onlyPlotSpecifiedFieldNames";
    static constexpr auto fieldNames = "fieldNames";
  } siloOutputViewKeys;
  /// @endcond

private:

  void postInputInitialization() override;

  string m_plotFileRoot;
  integer m_writeEdgeMesh;
  integer m_writeFaceMesh;
  integer m_writeCellElementMesh;
  integer m_writeFaceElementMesh;
  integer m_plotLevel;

  /// flag to decide whether we only plot the specified field names
  integer m_onlyPlotSpecifiedFieldNames;

  /// array of names of the fields to output
  array1d< string > m_fieldNames;

};


} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_SILOOUTPUT_HPP_ */
