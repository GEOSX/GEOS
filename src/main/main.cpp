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

// Source includes
#include "common/DataTypes.hpp"
#include "common/Format.hpp"
#include "common/TimingMacros.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"


using namespace geos;


int main( int argc, char *argv[] )
{
  try
  {
    std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();

    std::unique_ptr< CommandLineOptions > commandLineOptions = basicSetup( argc, argv, true );

    logger.rank0Log( GEOS_FMT( "Started at {:%Y-%m-%d %H:%M:%S}", startTime ) );

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

    logger.rank0Log( GEOS_FMT( "Finished at {:%Y-%m-%d %H:%M:%S}", endTime ) );
    logger.rank0Log( GEOS_FMT( "total time            {:%H:%M:%S}", totalTime ) );
    logger.rank0Log( GEOS_FMT( "initialization time   {:%H:%M:%S}", initTime ) );
    logger.rank0Log( GEOS_FMT( "run time              {:%H:%M:%S}", runTime ) );

    return 0;
  }
  // A NotAnError is thrown if "-h" or "--help" option is used.
  catch( NotAnError const & )
  {
    return 0;
  }
  catch( std::exception const & e )
  {
    logger.stdLog( e.what() );
    LvArray::system::callErrorHandler();
    std::abort();
  }
  return 0;
}
