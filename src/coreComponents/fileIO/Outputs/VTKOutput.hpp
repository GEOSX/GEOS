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
 * @file VTKOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_VTKOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_VTKOUTPUT_HPP_

#include "OutputBase.hpp"
#include "fileIO/vtk/VTKPolyDataWriterInterface.hpp"


namespace geos
{

/**
 * @brief A class for creating vtk outputs
 */
class VTKOutput : public OutputBase
{
public:

  /// @copydoc geos::dataRepository::Group::Group(string const & name, Group * const parent)
  VTKOutput( string const & name, Group * const parent );

  /// Destructor
  virtual ~VTKOutput() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "VTK"; }

  virtual void postInputInitialization() override;

  /**
   * @brief Set the plotFileRoot name for the output
   *
   * @param root The string name
   */
  void setPlotFileRoot( string const & root );

  /**
   * @brief Writes out a set of vtk files.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Write one final set of vtk files as the code exits
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

  /**
   * @brief Performs re-initialization of the datasets accumulated in the PVD writer.
   */
  virtual void reinit() override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeysStruct : OutputBase::viewKeysStruct
  {
    static constexpr auto plotFileRoot = "plotFileRoot";
    static constexpr auto writeFEMFaces = "writeFEMFaces";
    static constexpr auto writeGhostCells = "writeGhostCells";
    static constexpr auto writeFaceElementsAs3D = "writeFaceElementsAs3D";
    static constexpr auto plotLevel = "plotLevel";
    static constexpr auto binaryString = "format";
    static constexpr auto outputRegionTypeString = "outputRegionType";
    static constexpr auto onlyPlotSpecifiedFieldNames = "onlyPlotSpecifiedFieldNames";
    static constexpr auto fieldNames = "fieldNames";
    static constexpr auto levelNames = "levelNames";
  } vtkOutputViewKeys;
  /// @endcond

  /**
   * @brief Return PyVTKOutput type.
   * @return Return PyVTKOutput type.
   */
#if defined(GEOS_USE_PYGEOSX)
  virtual PyTypeObject * getPythonType() const override;
#endif

private:

  string m_plotFileRoot;
  integer m_writeFaceMesh;
  integer m_plotLevel;

  /// Should the vtk files contain the ghost cells or not.
  integer m_writeGhostCells;

  /// Should the face elements be written as 3d volumes or not.
  integer m_writeFaceElementsAs3D;

  /// flag to decide whether we only plot the specified field names
  integer m_onlyPlotSpecifiedFieldNames;

  /// array of names of the fields to output
  array1d< string > m_fieldNames;

  /// array of names of the mesh levels to output (an empty array means all levels are saved)
  array1d< string > m_levelNames;

  /// VTK output mode
  vtk::VTKOutputMode m_writeBinaryData = vtk::VTKOutputMode::BINARY;

  /// VTK output region filter
  vtk::VTKRegionTypes m_outputRegionType = vtk::VTKRegionTypes::ALL;

  vtk::VTKPolyDataWriterInterface m_writer;

};


} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_VTKOUTPUT_HPP_ */
