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

#include "common/Units.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::units;


struct DurationCase
{
  string m_expectedString;
  SystemClock::duration m_systemDuration;
  double m_simDuration;

  template< class DURATION >
  DurationCase( string_view expectedString, DURATION durationValue ):
    m_expectedString( expectedString ),
    m_systemDuration( std::chrono::duration_cast< SystemClock::duration >( durationValue ) ),
    m_simDuration( std::chrono::duration_cast< std::chrono::duration< double > >( durationValue ).count() )
  {}
};

TEST( Units, SystemDurationFormatTest )
{
  using namespace std::chrono;

  std::vector< DurationCase > durationCases = {

    DurationCase(
      "00h00m00s (1.11e-07 s)",
      nanoseconds( 111 ) ),

    DurationCase(
      "00h00m00s (0.000111 s)",
      microseconds( 111 ) ),

    DurationCase(
      "00h00m00s (0.111 s)",
      milliseconds( 111 ) ),

    DurationCase(
      "00h02m25s (145.016 s)",
      seconds( 145 ) + milliseconds( 16 ) ),

    DurationCase(
      "22h25m45s (80745.016 s)",
      hours( 20 ) + minutes( 145 ) + seconds( 45 ) + milliseconds( 16 ) ),

    DurationCase(
      "20d, 12h02m25s (1771345.016 s)",
      hours( long( 24 * 20.5 ) ) + seconds( 145 ) + milliseconds( 16 ) ),

    DurationCase(
      "20d, 12h02m25s (1771345.016 s)",
      Days( 20 ) + hours( 12 ) + seconds( 145 ) + milliseconds( 16 ) ),

    DurationCase(
      "1y, 0d, 00h00m00s (31556952 s)",
      Years( 1 ) ),

    DurationCase(
      "1y, 0d, 12h00m00s (31600152 s)",
      Years( 1 ) + hours( 12 ) ),

    DurationCase(
      "1y, 20d, 12h02m25s (33328297.016 s)",
      Years( 1 ) + hours( long( 24 * 20.5 ) ) + seconds( 145 ) + milliseconds( 16 ) ),

    DurationCase(
      "1y, 20d, 12h02m25s (33328297.016 s)",
      Years( 1 ) + Days( 20 ) + hours( 12 ) + seconds( 145 ) + milliseconds( 16 ) ),

    DurationCase(
      "12y, 362d, 05h49m12s (409981176 s)",
      Years( 13 ) - Days( 3 ) ),

    DurationCase(
      "13y, 365d, 05h48m12s (441797268 s)",
      Years( 14 ) - minutes( 1 ) ),

    DurationCase(
      "14y, 0d, 00h00m00s (441797328 s)",
      Years( 14 ) ),

    DurationCase(
      "14y, 0d, 12h00m00s (441840528 s)",
      Years( 14 ) + hours( 12 ) ),

    DurationCase(
      "100y, 20d, 12h02m25s (3157466545.016 s)",
      Years( 100 ) + Days( 20 ) + hours( 12 ) + seconds( 145 ) + milliseconds( 16 ) ),

    DurationCase(
      "292y, 20d, 12h02m25s (9216401329.016 s)",
      Years( 292 ) + Days( 20 ) + hours( 12 ) + seconds( 145 ) + milliseconds( 16 ) ),

    DurationCase(
      "5500y, 20d, 12h02m25s (173565007345 s)",
      Years( 5500 ) + Days( 20 ) + hours( 12 ) + seconds( 145 ) ),

    DurationCase(
      "-(00h00m00s) (-1.11e-07 s)",
      -nanoseconds( 111 ) ),

    DurationCase(
      "-(00h00m01s) (-1 s)",
      -seconds( 1 ) ),

    DurationCase(
      "-(00h02m25s) (-145.016 s)",
      -seconds( 145 ) - milliseconds( 16 ) ),

    DurationCase(
      "-(22h25m45s) (-80745.016 s)",
      -hours( 20 ) - minutes( 145 ) - seconds( 45 ) - milliseconds( 16 ) ),

    DurationCase(
      "-(5500y, 20d, 12h02m25s) (-173565007345 s)",
      -Years( 5500 ) - Days( 20 ) - hours( 12 ) - seconds( 145 ) ),

  };

  const SystemClock::duration maxDuration = SystemClock::duration::max();
  const string errorInfo = GEOS_FMT( "(Max possible duration = {} s)",
                                     duration_cast< seconds >( maxDuration ).count() );

  // Duration with more than 292 years are not supported by the SystemClock type.
  double maxSystemTime = duration_cast< seconds, double, std::ratio< 1 > >( SystemClock::duration::max() ).count();

  for( DurationCase const & durationCase : durationCases )
  {
    // testing "double" typed time (which has a limit that is way higher than the tests cases)
    EXPECT_STREQ( durationCase.m_expectedString.c_str(),
                  TimeFormatInfo::fromSeconds( durationCase.m_simDuration ).toString().c_str() ) << errorInfo;

    if( 0.0 < durationCase.m_simDuration && durationCase.m_simDuration <= maxSystemTime )
    {
      EXPECT_STREQ( durationCase.m_expectedString.c_str(),
                    TimeFormatInfo::fromDuration( durationCase.m_systemDuration ).toString().c_str() ) << errorInfo;
    }
  }
}
