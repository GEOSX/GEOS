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
 * @file GeosxState.hpp
 */

#ifndef GEOS_MAININTERFACE_GEOSXSTATE_HPP_
#define GEOS_MAININTERFACE_GEOSXSTATE_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"

// System includes
#include <functional>
#include <memory>

// Forward declaration of conduit::Node.
namespace conduit
{
class Node;
}

#if defined( GEOS_USE_CALIPER )
//Forward declaration of cali::ConfigManager.
namespace cali
{
class ConfigManager;
}
#endif


namespace geos
{

// Forward declarations.
namespace dataRepository
{
class Group;
}

class ProblemManager;
class FieldSpecificationManager;
class FunctionManager;
class CommunicationTools;
struct CommandLineOptions;

/**
 * @brief Return a string representing @c duration in seconds.
 * @param duration The duration to convert to a string.
 * @return A string representing @c duration in seconds.
 */
string durationToString( std::chrono::system_clock::duration const duration );

/**
 * @enum State
 * @brief A description of the global state.
 */
enum class State
{
  UNINITIALIZED = 0,
  INITIALIZED = 1,
  READY_TO_RUN = 2,
  COMPLETED = 3
};

/**
 * @brief Stream printer for a @c State.
 * @param os The stream to print to.
 * @param state The state to print.
 * @return @c os.
 */
std::ostream & operator<<( std::ostream & os, State const state );

/**
 * @class GeosxState
 * @brief Holds the global GEOSX state. This is a singleton class, no more than one instance exists at at time.
 *   After construction the instance can be accessed from anywhere by the free function @c getGlobalState().
 */
class GeosxState
{
public:

  /**
   * @brief Construct a new state from command line options.
   * @param commandLineOptions The command line options.
   * @details Loads in the restart file if applicable, and allocates the ProblemManager.
   * @post After this call @c getState() is @c State::UNINITIALIZED.
   */
  explicit GeosxState( std::unique_ptr< CommandLineOptions > && commandLineOptions );

  /**
   * @brief Destructor.
   * @note This is the default constructor but it must be specified in the `cpp`, otherwise the use of
   *   `std::unique_ptr` with forward declared types won't work.
   */
  ~GeosxState();

  GeosxState( GeosxState const & ) = delete;
  GeosxState( GeosxState && ) = delete;
  GeosxState & operator=( GeosxState const & ) = delete;
  GeosxState & operator=( GeosxState && ) = delete;

  /**
   * @brief Initialize the Data Repository from the input file.
   * @return True iff the problem is ready to be run.
   * @pre This must be called after construction when @c getState() is @c State::UNINITIALIZED.
   * @post After this call @c getState() is @c State::INITIALIZED.
   */
  bool initializeDataRepository();

  /**
   * @brief Apply initial conditions and if performing a restart overwrites the data repository.
   * @pre This must be called after @c initializeDataRepository() when @c getState()
   *   is @c State::INITIALIZED.
   * @post After this call @c getState() is @c State::READY_TO_RUN.
   */
  void applyInitialConditions();

  /**
   * @brief Run the simulation.
   * @pre This must be called after @c applyInitialConditions() when @c getState() is @c State::READY_TO_RUN.
   * @post Iff the simulation ran to completion then @c getState() is @c State::COMPLETED, otherwise if there is
   *   still stuff left to do @c getState() is @c State::READY_TO_RUN.
   */
  void run();

  /**
   * @brief Return the current State.
   * @return The current state.
   */
  State getState() const
  { return m_state; }

  /**
   * @brief Return the command line options.
   * @return The command line options.
   */
  CommandLineOptions const & getCommandLineOptions()
  {
    GEOS_ERROR_IF( m_commandLineOptions == nullptr, "Not initialized." );
    return *m_commandLineOptions;
  }

  /**
   * @brief Return the root conduit node.
   * @return The root conduit node.
   */
  conduit::Node & getRootConduitNode()
  {
    GEOS_ERROR_IF( m_rootNode == nullptr, "Not initialized." );
    return *m_rootNode;
  }

  /**
   * @brief Return the ProblemManager.
   * @return The ProblemManager.
   */
  ProblemManager & getProblemManager()
  {
    GEOS_ERROR_IF( m_problemManager == nullptr, "Not initialized." );
    return *m_problemManager;
  }

  /**
   * @brief Return the @c ProblemManager but as a @c dataRepository::Group.
   * @return The @c ProblemManager but as a @c dataRepository::Group.
   * @note This is useful if you only need at @c dataRepository::Group and don't want to
   *   include @c ProblemManager.hpp.
   */
  dataRepository::Group & getProblemManagerAsGroup();

  /**
   * @brief Return the FieldSpecificationManager.
   * @return The FieldSpecificationManager.
   */
  FieldSpecificationManager & getFieldSpecificationManager();

  /**
   * @brief Return the FunctionManager.
   * @return The FunctionManager.
   */
  FunctionManager & getFunctionManager();

  /**
   * @brief Return the CommunicationTools.
   * @return The CommunicationTools.
   */
  CommunicationTools & getCommunicationTools()
  {
    GEOS_ERROR_IF( m_commTools == nullptr, "Not initialized." );
    return *m_commTools;
  }

  /**
   * @brief Return the time taken to setup the problem.
   * @return The time taken to setup the problem.
   */
  std::chrono::system_clock::duration getInitTime() const
  { return m_initTime; }

  /**
   * @brief Return the time taken to run the problem.
   * @return The time taken to run the problem.
   */
  std::chrono::system_clock::duration getRunTime() const
  { return m_runTime; }

private:

  /// The current State.
  State m_state;

  /// The command line options.
  std::unique_ptr< CommandLineOptions > m_commandLineOptions;

  /// The root conduit node for the data repository.
  std::unique_ptr< conduit::Node > m_rootNode;

  /// The ProblemManager.
  std::unique_ptr< ProblemManager > m_problemManager;

  /// The CommunicationTools.
  std::unique_ptr< CommunicationTools > m_commTools;

#if defined( GEOS_USE_CALIPER )
  /// The Caliper ConfigManager.
  std::unique_ptr< cali::ConfigManager > m_caliperManager;
#endif

  /// The initialization time.
  std::chrono::system_clock::duration m_initTime;

  /// The run time.
  std::chrono::system_clock::duration m_runTime;
};

/**
 * @brief Return the current GeosxState.
 * @return The current GeosxState.
 */
GeosxState & getGlobalState();

} // namespace geos

#endif /* GEOS_MAININTERFACE_GEOSXSTATE_HPP_ */
