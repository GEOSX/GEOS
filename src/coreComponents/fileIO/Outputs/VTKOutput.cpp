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

/**
 * @file VTKOutput.cpp
 */

#include "VTKOutput.hpp"


#if defined(GEOS_USE_PYGEOSX)
#include "fileIO/python/PyVTKOutputType.hpp"
#endif

namespace geos
{

using namespace dataRepository;

VTKOutput::VTKOutput( string const & name,
                      Group * const parent ):
  OutputBase( name, parent ),
  m_plotFileRoot( name ),
  m_writeFaceMesh(),
  m_plotLevel(),
  m_onlyPlotSpecifiedFieldNames(),
  m_fieldNames(),
  m_levelNames(),
  m_writer( getOutputDirectory() + '/' + m_plotFileRoot )
{
  enableLogLevelInput();

  registerWrapper( viewKeysStruct::plotFileRoot, &m_plotFileRoot ).
    setDefaultValue( m_plotFileRoot ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the root file for this output." );

  registerWrapper( viewKeysStruct::writeFEMFaces, &m_writeFaceMesh ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "" );

  registerWrapper( viewKeysStruct::plotLevel, &m_plotLevel ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Level detail plot. Only fields with lower of equal plot level will be output." );

  registerWrapper( viewKeysStruct::writeGhostCells, &m_writeGhostCells ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Should the vtk files contain the ghost cells or not." );

  registerWrapper( viewKeysStruct::writeFaceElementsAs3D, &m_writeFaceElementsAs3D ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Should the face elements be written as 3d volumes or not." );

  registerWrapper( viewKeysStruct::onlyPlotSpecifiedFieldNames, &m_onlyPlotSpecifiedFieldNames ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "If this flag is equal to 1, then we only plot the fields listed in `fieldNames`. Otherwise, we plot all the fields with the required `plotLevel`, plus the fields listed in `fieldNames`" );

  registerWrapper( viewKeysStruct::fieldNames, &m_fieldNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the fields to output. If this attribute is specified, GEOSX outputs all the fields specified by the user, regardless of their `plotLevel`" );

  registerWrapper( viewKeysStruct::levelNames, &m_levelNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of mesh levels to output." );

  registerWrapper( viewKeysStruct::binaryString, &m_writeBinaryData ).
    setApplyDefaultValue( m_writeBinaryData ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output data format.  Valid options: ``" + EnumStrings< vtk::VTKOutputMode >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeysStruct::outputRegionTypeString, &m_outputRegionType ).
    setApplyDefaultValue( m_outputRegionType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output region types.  Valid options: ``" + EnumStrings< vtk::VTKRegionTypes >::concat( "``, ``" ) + "``" );
}

VTKOutput::~VTKOutput()
{}

void VTKOutput::postInputInitialization()
{
  m_writer.setOutputLocation( getOutputDirectory(), m_plotFileRoot );
  m_writer.setFieldNames( m_fieldNames.toViewConst() );
  m_writer.setLevelNames( m_levelNames.toViewConst() );
  m_writer.setOnlyPlotSpecifiedFieldNamesFlag( m_onlyPlotSpecifiedFieldNames );

  string const fieldNamesString = viewKeysStruct::fieldNames;
  string const onlyPlotSpecifiedFieldNamesString = viewKeysStruct::onlyPlotSpecifiedFieldNames;

  GEOS_THROW_IF( ( m_onlyPlotSpecifiedFieldNames != 0 ) && m_fieldNames.empty(),
                 GEOS_FMT( "{} `{}`: the flag `{}` is different from zero, but `{}` is empty, which is inconsistent",
                           catalogName(), getDataContext(),
                           onlyPlotSpecifiedFieldNamesString, fieldNamesString ),
                 InputError );

  GEOS_LOG_RANK_0_IF( !m_fieldNames.empty() && ( m_onlyPlotSpecifiedFieldNames != 0 ),
                      GEOS_FMT(
                        "{} `{}`: found {} fields to plot in `{}`. These fields will be output regardless of the `plotLevel` specified by the user. No other field will be output.",
                        catalogName(), getDataContext(),
                        std::to_string( m_fieldNames.size() ), fieldNamesString ) );

  GEOS_LOG_RANK_0_IF( !m_fieldNames.empty() && ( m_onlyPlotSpecifiedFieldNames == 0 ),
                      GEOS_FMT(
                        "{} `{}`: found {} fields to plot in `{}`, in addition to all fields with `plotLevel` smaller or equal to {}.",
                        catalogName(), getDataContext(),
                        std::to_string( m_fieldNames.size() ), fieldNamesString, m_plotLevel ) );

  GEOS_ERROR_IF( m_writeFaceElementsAs3D, GEOS_FMT( "{} `{}`: 3D vtk plot of faceElements is not yet supported.",
                                                    catalogName(), getDataContext() ) );
}


void VTKOutput::setPlotFileRoot( string const & root )
{
  m_plotFileRoot = root;
}


void VTKOutput::reinit()
{
  m_writer.clearData();
}

bool VTKOutput::execute( real64 const time_n,
                         real64 const dt,
                         integer const cycleNumber,
                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                         real64 const GEOS_UNUSED_PARAM ( eventProgress ),
                         DomainPartition & domain )
{
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: writing {} at time {} s (cycle number {})", getName(), m_fieldNames, time_n + dt, cycleNumber ));

  m_writer.setWriteGhostCells( m_writeGhostCells );
  m_writer.setWriteFaceElementsAs3D ( m_writeFaceElementsAs3D );
  m_writer.setOutputMode( m_writeBinaryData );
  m_writer.setOutputRegionType( m_outputRegionType );
  m_writer.setPlotLevel( m_plotLevel );
  m_writer.write( time_n, cycleNumber, domain );

  return false;
}

#if defined(GEOS_USE_PYGEOSX)
PyTypeObject * VTKOutput::getPythonType() const
{
  return python::getPyVTKOutputType();
}
#endif

REGISTER_CATALOG_ENTRY( OutputBase, VTKOutput, string const &, Group * const )
} /* namespace geos */
