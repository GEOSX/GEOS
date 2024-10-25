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

/*
 * @file SolverStatistics.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
#define GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP

#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * @class SolverStatistics
 * @brief This class is used to log the solver statistics
 */
class SolverStatistics : public dataRepository::Group
{
public:

  /**
   * @brief Constructor for SolverStatistics Objects.
   * @param[in] name the name of this instantiation of SolverStatistics in the repository
   * @param[in] parent the parent group of this instantiation of SolverStatistics
   */
  SolverStatistics( string const & name,
                    dataRepository::Group * const parent );

  /**
   * @brief Initialize the counters used for an individual time step
   */
  void initializeTimeStepStatistics();

  /**
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @param[in] numLinearIterations the number of linear iterations done by the linear solver
   * @detail This function is well suited for Newton's method, or for single-physics solvers in sequential schemes
   */
  void logNonlinearIteration( integer const numLinearIterations );

  /**
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @detail This function is well suited for the outer loop in sequential schemes
   */
  void logNonlinearIteration();

  /**
   * @brief Tell the solverStatistics that we are doing an outer loop iteration
   */
  void logOuterLoopIteration();

  /**
   * @brief Tell the solverStatistics that there is a time step cut
   */
  void logTimeStepCut();

  /**
   * @brief Save the statistics for the individual time step and increment the cumulative stats
   */
  void saveTimeStepStatistics();

  /**
   * @brief Output the cumulative statistics to the terminal
   */
  void outputStatistics() const;

  /**
   * @return Number of time steps
   */
  integer getNumTimeSteps() const
  { return m_numTimeSteps; }

  /**
   * @return Number of time step cuts
   */
  integer getNumTimeStepCuts() const
  { return m_numTimeStepCuts; }

  /**
   * @return Cumulative number of successful outer loop iterations
   */
  integer getNumSuccessfulOuterLoopIterations() const
  { return m_numSuccessfulOuterLoopIterations; }

  /**
   * @return Cumulative number of successful nonlinear iterations
   */
  integer getNumSuccessfulNonlinearIterations() const
  { return m_numSuccessfulNonlinearIterations; }

  /**
   * @return Cumulative number of successful linear iterations
   */
  integer getNumSuccessfulLinearIterations() const
  { return m_numSuccessfulLinearIterations; }

  /**
   * @return Cumulative number of discarded outer loop iterations
   */
  integer getNumDiscardedOuterLoopIterations() const
  { return m_numDiscardedOuterLoopIterations; }

  /**
   * @return Cumulative number of discarded nonlinear iterations
   */
  integer getNumDiscardedNonlinearIterations() const
  { return m_numDiscardedNonlinearIterations; }

  /**
   * @return Cumulative number of discarded linear iterations
   */
  integer getNumDiscardedLinearIterations() const
  { return m_numDiscardedLinearIterations; }


  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the number of time steps
    static constexpr char const * numTimeStepsString() { return "numTimeSteps"; }
    /// String key for the number of time step cuts
    static constexpr char const * numTimeStepCutsString() { return "numTimeStepCuts"; }

    /// String key for the successful number of outer loop iterations
    static constexpr char const * numSuccessfulOuterLoopIterationsString() { return "numSuccessfulOuterLoopIterations"; }
    /// String key for the successful number of nonlinear iterations
    static constexpr char const * numSuccessfulNonlinearIterationsString() { return "numSuccessfulNonlinearIterations"; }
    /// String key for the successful number of linear iterations
    static constexpr char const * numSuccessfulLinearIterationsString() { return "numSuccessfulLinearIterations"; }

    /// String key for the discarded number of outer loop iterations
    static constexpr char const * numDiscardedOuterLoopIterationsString() { return "numDiscardedOuterLoopIterations"; }
    /// String key for the discarded number of nonlinear iterations
    static constexpr char const * numDiscardedNonlinearIterationsString() { return "numDiscardedNonlinearIterations"; }
    /// String key for the discarded number of linear iterations
    static constexpr char const * numDiscardedLinearIterationsString() { return "numDiscardedLinearIterations"; }
  };

private:
  /// Number of time steps
  integer m_numTimeSteps;

  /// Number of time step cuts
  integer m_numTimeStepCuts;


  /// Number of outer loop iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumOuterLoopIterations;

  /// Number of nonlinear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumNonlinearIterations;

  /// Number of linear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumLinearIterations;


  /// Cumulative number of successful outer loop iterations
  integer m_numSuccessfulOuterLoopIterations;

  /// Cumulative number of successful nonlinear iterations
  integer m_numSuccessfulNonlinearIterations;

  /// Cumulative number of successful linear iterations
  integer m_numSuccessfulLinearIterations;


  /// Cumulative number of discarded outer loop iterations
  integer m_numDiscardedOuterLoopIterations;

  /// Cumulative number of discarded nonlinear iterations
  integer m_numDiscardedNonlinearIterations;

  /// Cumulative number of discarded linear iterations
  integer m_numDiscardedLinearIterations;

};

} //namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
