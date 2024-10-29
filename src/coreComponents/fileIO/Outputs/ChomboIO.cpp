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
 * @file ChomboIO.cpp
 */
#include "ChomboIO.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/DomainPartition.hpp"
#include "fileIO/coupling/ChomboCoupler.hpp"

#include <fstream>
#include <chrono>

namespace geos
{

using namespace dataRepository;

ChomboIO::ChomboIO( string const & name, Group * const parent ):
  OutputBase( name, parent ),
  m_coupler( nullptr ),
  m_outputPath(),
  m_beginCycle( 0 ),
  m_inputPath( "/INVALID_INPUT_PATH" ),
  m_waitForInput(),
  m_useChomboPressures()
{
  registerWrapper( viewKeyStruct::outputPathString(), &m_outputPath ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Path at which the geos to chombo file will be written." );

  registerWrapper( viewKeyStruct::beginCycleString(), &m_beginCycle ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Cycle at which the coupling will commence." );

  registerWrapper( viewKeyStruct::inputPathString(), &m_inputPath ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( "/INVALID_INPUT_PATH" ).
    setDescription( "Path at which the chombo to geos file will be written." );

  registerWrapper( viewKeyStruct::waitForInputString(), &m_waitForInput ).
    setInputFlag( InputFlags::REQUIRED ).
    setDefaultValue( 0 ).
    setDescription( "True iff geos should wait for chombo to write out a file. When true the inputPath must be set." );

  registerWrapper( viewKeyStruct::useChomboPressuresString(), &m_useChomboPressures ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "True iff geos should use the pressures chombo writes out." );
}

ChomboIO::~ChomboIO()
{
  delete m_coupler;
  m_coupler = nullptr;
}

bool ChomboIO::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                        real64 const dt,
                        integer const cycleNumber,
                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                        DomainPartition & domain )
{
  if( cycleNumber < m_beginCycle )
  {
    return false;
  }
  GEOS_LOG_RANK_0( "Executing chombo coupling" );

  if( m_coupler == nullptr )
  {
    GEOS_ERROR_IF( m_waitForInput && m_inputPath == "/INVALID_INPUT_PATH", "Waiting for input but no input path was specified." );

    GEOS_LOG_LEVEL_RANK_0( 1, "Initializing chombo coupling" );

    m_coupler = new ChomboCoupler( MPI_COMM_GEOSX, m_outputPath, m_inputPath, domain.getMeshBody( 0 ).getBaseDiscretization() );

  }

  m_coupler->write( dt );

  if( m_waitForInput )
  {
    m_coupler->read( m_useChomboPressures );
  }

  return false;
}

REGISTER_CATALOG_ENTRY( OutputBase, ChomboIO, string const &, Group * const )
} /* namespace geos */
