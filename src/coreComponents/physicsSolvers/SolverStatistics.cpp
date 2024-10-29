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

/*
 * @file SolverStatistics.cpp
 */

#include "SolverStatistics.hpp"

namespace geos
{

using namespace dataRepository;

SolverStatistics::SolverStatistics( string const & name, Group * const parent )
  : Group( name, parent ),
  m_currentNumOuterLoopIterations( 0 ),
  m_currentNumNonlinearIterations( 0 ),
  m_currentNumLinearIterations( 0 )
{
  registerWrapper( viewKeyStruct::numTimeStepsString(), &m_numTimeSteps ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time steps" );

  registerWrapper( viewKeyStruct::numTimeStepCutsString(), &m_numTimeStepCuts ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time step cuts" );


  registerWrapper( viewKeyStruct::numSuccessfulOuterLoopIterationsString(), &m_numSuccessfulOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful outer loop iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulNonlinearIterationsString(), &m_numSuccessfulNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful nonlinear iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulLinearIterationsString(), &m_numSuccessfulLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful linear iterations" );


  registerWrapper( viewKeyStruct::numDiscardedOuterLoopIterationsString(), &m_numDiscardedOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded outer loop iterations" );

  registerWrapper( viewKeyStruct::numDiscardedNonlinearIterationsString(), &m_numDiscardedNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded nonlinear iterations" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_numDiscardedLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded linear iterations" );
}

void SolverStatistics::initializeTimeStepStatistics()
{
  // the time step begins, we reset the individual-timestep counters
  m_currentNumOuterLoopIterations = 0;
  m_currentNumNonlinearIterations = 0;
  m_currentNumLinearIterations = 0;
}

void SolverStatistics::logNonlinearIteration( integer const numLinearIterations )
{
  // we have just performed a Newton iteration, so we increment the individual-timestep counters
  m_currentNumNonlinearIterations++;
  m_currentNumLinearIterations += numLinearIterations;
}

void SolverStatistics::logNonlinearIteration()
{
  // we have just performed an outer iteration, so we increment the individual-timestep counter (number of outer iteration)
  m_currentNumNonlinearIterations++;
}

void SolverStatistics::logOuterLoopIteration()
{
  // we have just performed an outer loop iteration, so we increment the individual-timestep counter for outer loop iterations
  m_currentNumOuterLoopIterations++;
}


void SolverStatistics::logTimeStepCut()
{
  // we have just cut the time step, so we increment the cumulative counters for discarded timesteps
  m_numDiscardedOuterLoopIterations += m_currentNumOuterLoopIterations;
  m_numDiscardedNonlinearIterations += m_currentNumNonlinearIterations;
  m_numDiscardedLinearIterations += m_currentNumLinearIterations;
  m_numTimeStepCuts++;

  // we are going to restart the timestep from the previous converged time step, so we have to re-initialize the statistics
  initializeTimeStepStatistics();
}

void SolverStatistics::saveTimeStepStatistics()
{
  // the timestep has converged, so we increment the cumulative counters for successful timesteps
  m_numSuccessfulOuterLoopIterations += m_currentNumOuterLoopIterations;
  m_numSuccessfulNonlinearIterations += m_currentNumNonlinearIterations;
  m_numSuccessfulLinearIterations += m_currentNumLinearIterations;
  m_numTimeSteps++;
}

void SolverStatistics::outputStatistics() const
{
  bool const printOuterLoopIterations = !(m_numSuccessfulOuterLoopIterations == 0 && m_numDiscardedOuterLoopIterations == 0);
  bool const printIterations = !(m_numSuccessfulNonlinearIterations == 0 && m_numDiscardedNonlinearIterations == 0);
  bool const printLinearIterations = !(m_numSuccessfulLinearIterations == 0 && m_numDiscardedLinearIterations == 0);

  auto const logStat = [&]( auto const name, auto const value )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "{}, number of {}: {}",
                               getParent().getName(), name, value ) );
  };

  // TODO: the print logic is really convoluted to accomodate the needs of the different solvers, needs simplification

  if( m_numTimeSteps > 0 )
  {
    logStat( "time steps", m_numTimeSteps );
  }

  if( printIterations )
  {
    if( printOuterLoopIterations )
    {
      logStat( "successful outer loop iterations", m_numSuccessfulOuterLoopIterations );
    }
    logStat( "successful nonlinear iterations", m_numSuccessfulNonlinearIterations );
    if( printLinearIterations ) // don't print for the outer iterations in sequential schemes
    {
      logStat( "successful linear iterations", m_numSuccessfulLinearIterations );
    }

    logStat( "time step cuts", m_numTimeStepCuts );
    if( printOuterLoopIterations )
    {
      logStat( "discarded outer loop iterations", m_numDiscardedOuterLoopIterations );
    }
    logStat( "discarded nonlinear iterations", m_numDiscardedNonlinearIterations );
    if( printLinearIterations ) // don't print for the outer iterations in sequential schemes
    {
      logStat( "discarded linear iterations", m_numDiscardedLinearIterations );
    }
  }
}
} // namespace geos
