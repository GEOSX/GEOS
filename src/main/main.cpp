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

// Source includes
#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "common/TimingMacros.hpp"
#include "common/Units.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/version.hpp"


using namespace geos;


int main( int argc, char *argv[] )
{
  try
  {
    std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();

    std::unique_ptr< CommandLineOptions > commandLineOptions = basicSetup( argc, argv, true );

    outputVersionInfo();

    GEOS_LOG_RANK_0( GEOS_FMT( "Started at {:%Y-%m-%d %H:%M:%S}", startTime ) );

    std::chrono::system_clock::duration initTime;
    std::chrono::system_clock::duration runTime;
    {
      GeosxState state( std::move( commandLineOptions ) );

      bool const problemToRun = state.initializeDataRepository();
      if( problemToRun )
      {
        state.applyInitialConditions();
        state.run();
        LVARRAY_WARNING_IF( state.getState() != State::COMPLETED, "Simulation exited early." );
      }

      initTime = state.getInitTime();
      runTime = state.getRunTime();
    }

    basicCleanup();

    std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();
    std::chrono::system_clock::duration totalTime = endTime - startTime;

    GEOS_LOG_RANK_0( GEOS_FMT( "Finished at {:%Y-%m-%d %H:%M:%S}", endTime ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "total time            {}", units::TimeFormatInfo::fromDuration( totalTime ) ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "initialization time   {}", units::TimeFormatInfo::fromDuration( initTime ) ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "run time              {}", units::TimeFormatInfo::fromDuration( runTime ) ) );

    return 0;
  }
  // A NotAnError is thrown if "-h" or "--help" option is used.
  catch( NotAnError const & )
  {
    basicCleanup();
    return 0;
  }
  catch( std::exception const & e )
  {
    GEOS_LOG( e.what() );
    LvArray::system::callErrorHandler();
    basicCleanup();
    std::abort();
  }
  return 0;
}
