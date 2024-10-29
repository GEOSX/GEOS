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
 * @file FieldStatisticsBase.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDSTATISTICSBASE_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDSTATISTICSBASE_HPP_

#include "events/tasks/TaskBase.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshLevel.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geos
{

/**
 * @class FieldStatisticsBase
 *
 * Base task class allowing for the computation of aggregate statistics
 */
template< typename SOLVER >
class FieldStatisticsBase : public TaskBase
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  FieldStatisticsBase( const string & name,
                       Group * const parent )
    : TaskBase( name, parent ),
    m_solver( nullptr ),
    m_outputDir( joinPath( OutputBase::getOutputDirectory(), name ) )
  {
    enableLogLevelInput();

    string const key = SOLVER::coupledSolverAttributePrefix() + "SolverName";
    registerWrapper( key, &m_solverName ).
      setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
      setInputFlag( dataRepository::InputFlags::REQUIRED ).
      setDescription( "Name of the " + SOLVER::coupledSolverAttributePrefix() + " solver" );

    this->registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
      setApplyDefaultValue( 0 ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setDescription( "Write statistics into a CSV file" );
  }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override = 0;

  /**@}*/

protected:

  void postInputInitialization() override
  {
    ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
    PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

    m_solver = physicsSolverManager.getGroupPointer< SOLVER >( m_solverName );
    GEOS_THROW_IF( m_solver == nullptr,
                   GEOS_FMT( "{}: Could not find solver '{}' of type {}",
                             getDataContext(),
                             m_solverName, LvArray::system::demangleType< SOLVER >() ),
                   InputError );

    // create dir for output
    if( m_writeCSV > 0 )
    {
      if( MpiWrapper::commRank() == 0 )
      {
        makeDirsForPath( m_outputDir );
      }
      // wait till the dir is created by rank 0
      MPI_Barrier( MPI_COMM_WORLD );
    }
  }

  struct viewKeyStruct
  {
    static constexpr char const * writeCSVFlagString() { return "writeCSV"; }
  };

  /// Pointer to the physics solver
  SOLVER * m_solver;

  // Output directory
  string const m_outputDir;

  // Flag to enable writing CSV output
  integer m_writeCSV;

private:

  /// Name of the solver
  string m_solverName;

};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDSTATISTICSBASE_HPP_ */
