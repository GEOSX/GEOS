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

#ifndef GEOS_FILEIO_VTK_VTKPOLYDATAWRITERINTERFACE_HPP_
#define GEOS_FILEIO_VTK_VTKPOLYDATAWRITERINTERFACE_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "dataRepository/Wrapper.hpp"
#include "fileIO/vtk/VTKPVDWriter.hpp"
#include "fileIO/vtk/VTKVTMWriter.hpp"
#include "codingUtilities/EnumStrings.hpp"

class vtkUnstructuredGrid;
class vtkPointData;
class vtkCellData;

namespace geos
{

class DomainPartition;
class ElementRegionBase;
class ParticleRegionBase;
class EmbeddedSurfaceNodeManager;
class ElementRegionManager;
class NodeManager;
class ParticleManager;
class FaceManager;

namespace vtk
{

enum struct VTKOutputMode
{
  BINARY,
  ASCII
};

enum struct VTKRegionTypes
{
  CELL,
  WELL,
  SURFACE,
  PARTICLE,
  ALL
};

/// Declare strings associated with output enumeration values.
ENUM_STRINGS( VTKOutputMode,
              "binary",
              "ascii" );

/// Declare strings associated with region type enumeration values.
ENUM_STRINGS( VTKRegionTypes,
              "cell",
              "well",
              "surface",
              "particle",
              "all" );

/**
 * @brief Encapsulate output methods for vtk
 */
class VTKPolyDataWriterInterface
{
public:
  /**
   * @brief Constructor
   * @param[in] outputName folder name in which all the files will be written
   */
  explicit VTKPolyDataWriterInterface( string outputName );

  /**
   * @brief Defines if the vtk outputs should contain the ghost cells.
   * @param writeGhostCells The boolean flag.
   */
  void setWriteGhostCells( bool writeGhostCells )
  {
    m_writeGhostCells = writeGhostCells;
  }

  /**
   * @brief Defines whether in the vtk output facelements should be 2D or 3D
   * @param writeFaceElementsAs3D The boolean flag.
   */
  void setWriteFaceElementsAs3D( bool writeFaceElementsAs3D )
  {
    m_writeFaceElementsAs3D = writeFaceElementsAs3D;
  }

  /**
   * @brief Sets the plot level
   * @details All fields have an associated plot level. If it is <= to \p plotLevel,
   * the field will be output.
   * @param[in] plotLevel the limit plotlevel
   */
  void setPlotLevel( integer plotLevel )
  {
    m_plotLevel = dataRepository::toPlotLevel( plotLevel );
  }

  /**
   * @brief Set the binary mode
   * @param[in] mode output mode to be set
   */
  void setOutputMode( VTKOutputMode mode )
  {
    m_outputMode = mode;
  }

  /**
   * @brief Set the output region type
   * @param[in] regionType output region type to be set
   */
  void setOutputRegionType( VTKRegionTypes regionType )
  {
    m_outputRegionType = regionType;
  }

  /**
   * @brief Set the output directory name
   * @param[in] outputDir global output directory location
   * @param[in] outputName name of the VTK output subdirectory and corresponding PVD file
   */
  void setOutputLocation( string outputDir, string outputName )
  {
    m_outputDir = std::move( outputDir );
    m_outputName = std::move( outputName );
    m_pvd.setFileName( joinPath( m_outputDir, m_outputName ) + ".pvd" );
  }

  /**
   * @brief Set the flag to decide whether we only plot the fields specified by fieldNames, or if we also plot fields based on plotLevel
   * @param[in] onlyPlotSpecifiedFieldNames the flag
   */
  void setOnlyPlotSpecifiedFieldNamesFlag( integer const onlyPlotSpecifiedFieldNames )
  {
    m_onlyPlotSpecifiedFieldNames = onlyPlotSpecifiedFieldNames;
  }

  /**
   * @brief Set the names of the fields to output
   * @param[in] fieldNames the fields to output
   */
  void setFieldNames( arrayView1d< string const > const & fieldNames )
  {
    m_fieldNames.insert( fieldNames.begin(), fieldNames.end() );
  }

  /**
   * @brief Set the names of the mesh levels to output
   * @param[in] levelNames the mesh levels to output (an empty array means all levels are saved)
   */
  void setLevelNames( arrayView1d< string const > const & levelNames )
  {
    m_levelNames.insert( levelNames.begin(), levelNames.end() );
  }

  /**
   * @brief Main method of this class. Write all the files for one time step.
   * @details This method writes a .pvd file (if a previous one was created from a precedent time step,
   * it is overwritten). The .pvd file contains relative path to every .vtm files (one vtm file per time step).
   * This method triggers also the writing of a .vtm file. A .vtm file containts relative paths to blocks
   * with the following hierarchy :
   *  - CellElementRegion
   *    - CellElementRegion1
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - CellElementRegion2
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - ...
   *  -WellElementRegion
   *    - Well1
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - Well2
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   * @param[in] time the time step to be written
   * @param[in] cycle the current cycle of event
   * @param[in] domain the computation domain of this rank
   */
  void write( real64 time, integer cycle, DomainPartition const & domain );

  /**
   * @brief Clears the datasets accumulated in the pvd writer
   *
   */
  void clearData();


private:

  /**
   * @brief Check if plotting is enabled for this field
   * @param[in] wrapper the wrapper
   * @return true if this wrapper should be plot, false otherwise
   */
  bool isFieldPlotEnabled( dataRepository::WrapperBase const & wrapper ) const;

  /**
   * @brief Writes the files for all the CellElementRegions.
   * @details There will be one file written per CellElementRegion and per rank.
   * @param[in] time the time-step
   * @param[in] cycle the current cycle number
   * @param[in] elemManager the ElementRegionManager containing the CellElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   * @param[in] meshLevelName the name of the MeshLevel containing the nodes and elements to be output
   * @param[in] meshBodyName the name of the MeshBody containing the nodes and elements to be output
   */
  void writeCellElementRegions( real64 time,
                                ElementRegionManager const & elemManager,
                                NodeManager const & nodeManager,
                                string const & path ) const;

  void writeParticleRegions( real64 const time,
                             ParticleManager const & particleManager,
                             string const & path ) const;

  /**
   * @brief Writes the files containing the well representation
   * @details There will be one file written per WellElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] cycle the current cycle number
   * @param[in] elemManager the ElementRegionManager containing the WellElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   */
  void writeWellElementRegions( real64 time,
                                ElementRegionManager const & elemManager,
                                NodeManager const & nodeManager,
                                string const & path ) const;

  /**
   * @brief Writes the files containing the faces elements
   * @details There will be one file written per FaceElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] elemManager the ElementRegionManager containing the FaceElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   * @param[in] embSurfNodeManager the embedded surface node Manager.
   * @param[in] faceManager the faceManager.
   * @param[in] path the path to the output file.
   */
  void writeSurfaceElementRegions( real64 time,
                                   ElementRegionManager const & elemManager,
                                   NodeManager const & nodeManager,
                                   EmbeddedSurfaceNodeManager const & embSurfNodeManager,
                                   FaceManager const & faceManager,
                                   string const & path ) const;

  /**
   * @brief Writes a VTM file for the time-step \p time.
   * @details a VTM file is a VTK Multiblock file. It contains relative path to different files organized in blocks.
   * @param[in] cycle the current cycle number
   * @param[in] elemManager the ElementRegionManager containing all the regions to be output and referred to in the VTM file
   * @param[in] vtmWriter a writer specialized for the VTM file format
   */

  void writeVtmFile( integer const cycle,
                     DomainPartition const & domain,
                     VTKVTMWriter const & vtmWriter ) const;
  /**
   * @brief Write all the fields associated to the nodes of \p nodeManager if their plotlevel is <= m_plotLevel
   * @param[in] pointData a VTK object containing all the fields associated with the nodes
   * @param[in] nodeIndices list of local node indices to write
   * @param[in] nodeManager the NodeManager associated with the domain being written
   */
  void writeNodeFields( NodeManager const & nodeManager,
                        arrayView1d< localIndex const > const & nodeIndices,
                        vtkPointData * pointData ) const;

  /**
   * @brief Writes all the fields associated to the elements of \p er if their plotlevel is <= m_plotLevel
   * @param[in] subRegion ElementRegion being written
   * @param[in] cellData a VTK object containing all the fields associated with the elements
   */
  void writeElementFields( ElementRegionBase const & subRegion,
                           vtkCellData * cellData ) const;

  void writeParticleFields( ParticleRegionBase const & region,
                            vtkCellData * cellData ) const;

  /**
   * @brief Writes an unstructured grid
   * @details The unstructured grid is the last element in the hierarchy of the output,
   * it contains the cells connectivities and the vertices coordinates as long as the
   * data fields associated with it
   * @param[in] ug a VTK SmartPointer to the VTK unstructured grid.
   * @param[in] path directory path for the grid file
   */
  void writeUnstructuredGrid( string const & path,
                              vtkUnstructuredGrid * ug ) const;

private:

  /// Output directory name
  string m_outputDir;

  /// Name of the PVD file and associated directory containing the outputs
  string m_outputName;

  /// A writter specialized for PVD files. There is one PVD file per simulation. It is the root
  /// file containing all the paths to the VTM files.
  VTKPVDWriter m_pvd;

  /// Should the vtk files contain the ghost cells or not.
  bool m_writeGhostCells;

  /// Maximum plot level to be written.
  dataRepository::PlotLevel m_plotLevel;

  /// Flag to decide whether we only plot the fields specified by fieldNames, or if we also plot fields based on plotLevel
  integer m_onlyPlotSpecifiedFieldNames;

  /// Flag to decide whether we check that the specified fieldNames are actually registered
  bool m_requireFieldRegistrationCheck;

  /// Names of the fields to output
  std::set< string > m_fieldNames;

  /// Names of the mesh levels to output (an empty array means all levels are saved)
  std::set< string > m_levelNames;

  /// The previousCycle
  integer m_previousCycle;

  /// Output mode, could be ASCII or BINARAY
  VTKOutputMode m_outputMode;

  /// Region output type, could be CELL, WELL, SURFACE, or ALL
  VTKRegionTypes m_outputRegionType;

  /// Defines whether to plot a faceElement as a 3D volumetric element or not.
  bool m_writeFaceElementsAs3D;
};

} // namespace vtk
} // namespace geos

#endif
