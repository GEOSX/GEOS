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
 * @file PhysicsSolverBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_PHYSICSSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_PHYSICSSOLVERBASE_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mesh/MeshBody.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"
#include "physicsSolvers/LinearSolverParameters.hpp"
#include "physicsSolvers/SolverStatistics.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"

#include <limits>

namespace geos
{

class DomainPartition;

/**
 * @class PhysicsSolverBase
 * @brief Base class for all physics solvers
 *
 * This class provides the base interface for all physics solvers. It provides the basic
 * functionality for setting up and solving a linear system, as well as the interface for
 * performing a timestep.
 */
class PhysicsSolverBase : public ExecutableGroup
{
public:

  /**
   * @brief Constructor for PhysicsSolverBase
   * @param name the name of this instantiation of PhysicsSolverBase
   * @param parent the parent group of this instantiation of PhysicsSolverBase
   */
  explicit PhysicsSolverBase( string const & name,
                              Group * const parent );

  /**
   * @brief Move constructor for PhysicsSolverBase
   */
  PhysicsSolverBase( PhysicsSolverBase && ) = default;

  /**
   * @brief Destructor for PhysicsSolverBase
   */
  virtual ~PhysicsSolverBase() override;

  /**
   * @brief Deleted constructor
   */
  PhysicsSolverBase() = delete;

  /**
   * @brief Deleted copy constructor
   */
  PhysicsSolverBase( PhysicsSolverBase const & ) = delete;

  /**
   * @brief Deleted copy assignment operator
   */
  PhysicsSolverBase & operator=( PhysicsSolverBase const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   */
  PhysicsSolverBase & operator=( PhysicsSolverBase && ) = delete;

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const = 0;


  /**
   * @brief Register wrappers that contain data on the mesh objects
   * @param MeshBodies the group of mesh bodies
   */
  virtual void registerDataOnMesh( Group & MeshBodies ) override;

  /**
   * @brief Initialization tasks after mesh generation is completed.
   */
  virtual void initialize_postMeshGeneration() override;

  /**
   * @brief Generate mesh targets from target regions
   * @param meshBodies the group of mesh bodies
   */
  void generateMeshTargetsFromTargetRegions( Group const & meshBodies );

  /**
   * @copydoc ExecutableGroup::cleanup
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @copydoc ExecutableGroup::execute
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Getter for system matrix
   * @return a reference to linear system matrix of this solver
   */
  ParallelMatrix & getSystemMatrix() { return m_matrix; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelMatrix const & getSystemMatrix() const { return m_matrix; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelVector & getSystemRhs() { return m_rhs; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelVector const & getSystemRhs() const { return m_rhs; }

  /**
   * @brief Getter for system solution vector
   * @return a reference to solution vector of this solver
   */
  ParallelVector & getSystemSolution() { return m_solution; }

  /**
   * @brief Getter for system solution vector
   * @return a reference to solution vector of this solver
   */
  ParallelVector const & getSystemSolution() const { return m_solution; }

  /**
   * @brief Getter for degree-of-freedom manager
   * @return a reference to degree-of-freedom manager of this solver
   */
  DofManager & getDofManager() { return m_dofManager; }

  /**
   * @brief Getter for degree-of-freedom manager
   * @return a reference to degree-of-freedom manager of this solver
   */
  DofManager const & getDofManager() const { return m_dofManager; }

  /**
   * @brief Getter for local matrix
   * @return a reference to linear system matrix of this solver
   */
  CRSMatrix< real64, globalIndex > & getLocalMatrix() { return m_localMatrix; }

  /**
   * @brief Getter for local matrix
   * @return a reference to linear system matrix of this solver
   */
  CRSMatrixView< real64 const, globalIndex const > getLocalMatrix() const { return m_localMatrix.toViewConst(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  /**
   * @brief entry function to perform a solver step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   *
   * This function is the entry point to perform a solver step. The choice of time integration
   * method is determined in this function, and the appropriate step function is called.
   */
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain );

  /**
   * @brief function to set the next time step size
   * @param[in] currentTime the current time
   * @param[in] lastDt the last time step size
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */
  virtual real64 setNextDt( real64 const & currentTime,
                            real64 const & lastDt,
                            DomainPartition & domain );

  /**
   * @brief function to set the next time step size based on Newton convergence
   * @param[in] lastDt the last time step size
   * @return the prescribed time step size
   */
  virtual real64 setNextDtBasedOnNewtonIter( real64 const & lastDt );

  /**
   * @brief function to set the next dt based on state change
   * @param [in]  lastDt the last time step size
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */
  virtual real64 setNextDtBasedOnStateChange( real64 const & lastDt,
                                              DomainPartition & domain );

  /**
   * @brief Entry function for an explicit time integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   */
  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain );

  /**
   * @brief Function for a nonlinear implicit integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   *
   * This function implements a nonlinear newton method for implicit problems. It requires that the
   * other functions in the solver interface are implemented in the derived physics solver. The
   * nonlinear loop includes a simple line search algorithm, and will cut the timestep if
   * convergence is not achieved according to the parameters in linearSolverParameters member.
   */
  virtual real64 nonlinearImplicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition & domain );

  /**
   * @brief Function to perform line search
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param scaleFactor the scaling factor to apply to the solution
   * @param lastResidual (in) target value below which to reduce residual norm, (out) achieved residual norm
   * @return return true if line search succeeded, false otherwise
   *
   * This function implements a nonlinear newton method for implicit problems. It requires that the
   * other functions in the solver interface are implemented in the derived physics solver. The
   * nonlinear loop includes a simple line search algorithm, and will cut the timestep if
   * convergence is not achieved according to the parameters in linearSolverParameters member.
   */
  virtual bool
  lineSearch( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain,
              DofManager const & dofManager,
              CRSMatrixView< real64, globalIndex const > const & localMatrix,
              ParallelVector & rhs,
              ParallelVector & solution,
              real64 const scaleFactor,
              real64 & lastResidual );

  /**
   * @brief Function to perform line search using a parabolic interpolation to find the scaling factor.
   * @param time_n time at the beginning of the step
   * @param dt the prescribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param scaleFactor the scaling factor to apply to the solution
   * @param lastResidual (in) target value below which to reduce residual norm, (out) achieved residual norm
   * @param residualNormT the residual norm at the end of the line search
   * @return return true if line search succeeded, false otherwise
   *
   */
  virtual bool
  lineSearchWithParabolicInterpolation ( real64 const & time_n,
                                         real64 const & dt,
                                         integer const cycleNumber,
                                         DomainPartition & domain,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         ParallelVector & rhs,
                                         ParallelVector & solution,
                                         real64 const scaleFactor,
                                         real64 & lastResidual,
                                         real64 & residualNormT );

  /**
   * @brief Function for a linear implicit integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   *
   * This function implements a single linear step. Similar to the nonlinear step, however it
   * assumes that the solution is achieved in a iteration. The use of this function requires that
   * the other functions in the solver interface are implemented in the derived physics solver. The
   * nonlinear loop includes a simple line search algorithm, and will cut the timestep if
   * convergence is not achieved according to the parameters in linearSolverParameters member.
   */
  virtual real64 linearImplicitStep( real64 const & time_n,
                                     real64 const & dt,
                                     integer const cycleNumber,
                                     DomainPartition & domain );

  /**
   * @brief function to perform setup for implicit timestep
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   *
   * This function should contain any step level initialization required to perform an implicit
   * step.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain );

  /**
   * @brief Populate degree-of-freedom manager with fields relevant to this solver
   * @param domain the domain containing the mesh and fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   */
  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const;

  /**
   * @brief Set up the linear system (DOF indices and sparsity patterns)
   * @param domain the domain containing the mesh and fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param setSparsity flag to indicate if the sparsity pattern should be set
   *
   * @note While the function is virtual, the base class implementation should be
   *       sufficient for most single-physics solvers.
   */
  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true );

  /**
   * @brief function to assemble the linear system matrix and rhs
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   *
   * This function assembles the residual and the jacobian of the residual wrt the primary
   * variables. In a stand alone physics solver, this function will fill a single block in the
   * block system. However the capability to query the block system structure for any coupled blocks
   * may be implemented to fill in off diagonal blocks of the system to enable coupling between
   * solvers.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs );

  /**
   * @brief apply boundary condition to system
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   *
   * This function applies all boundary conditions to the linear system. This is essentially a
   * completion of the system assembly, but is separated for use in coupled solvers.
   */
  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs );

  /**
   * @brief Output the assembled linear system for debug purposes.
   * @param time beginning-of-step time
   * @param cycleNumber event cycle number
   * @param nonlinearIteration current nonlinear iteration number
   * @param matrix system matrix
   * @param rhs system right-hand side vector
   */
  void
  debugOutputSystem( real64 const & time,
                     integer const cycleNumber,
                     integer const nonlinearIteration,
                     ParallelMatrix const & matrix,
                     ParallelVector const & rhs ) const;

  /**
   * @brief Output the linear system solution for debug purposes.
   * @param time beginning-of-step time
   * @param cycleNumber event cycle number
   * @param nonlinearIteration current nonlinear iteration number
   * @param solution system solution vector
   */
  void
  debugOutputSolution( real64 const & time,
                       integer const cycleNumber,
                       integer const nonlinearIteration,
                       ParallelVector const & solution ) const;

  /**
   * @brief calculate the norm of the global system residual
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localRhs the system right-hand side vector
   * @return norm of the residual
   *
   * This function returns the norm of global residual vector, which is suitable for comparison with
   * a tolerance.
   */
  virtual real64
  calculateResidualNorm( real64 const & time,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs );

  /**
   * @brief function to apply a linear system solver to the assembled system.
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   *
   * This function calls the linear solver package to perform a single linear solve on the block
   * system. The derived physics solver is required to specify the call, as no default is provided.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  solveLinearSystem( DofManager const & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution );

  /**
   * @brief Function to check system solution for physical consistency and constraint violation
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localSolution the solution vector
   * @param scalingFactor factor to scale the solution prior to application
   * @return true if solution can be safely applied without violating physical constraints, false otherwise
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   *
   */
  virtual bool
  checkSystemSolution( DomainPartition & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor );

  /**
   * @brief Function to determine if the solution vector should be scaled back in order to maintain a known constraint.
   * @param[in] domain The domain partition.
   * @param[in] dofManager degree-of-freedom manager associated with the linear system
   * @param[in] localSolution the solution vector
   * @return The factor that should be used to scale the solution vector values when they are being applied.
   */
  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution );

  /**
   * @brief Function to apply the solution vector to the state
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localSolution the solution vector
   * @param scalingFactor factor to scale the solution prior to application
   * @param dt the timestep
   * @param domain the domain partition
   *
   * This function performs 2 operations:
   * 1) extract the solution vector for the "blockSystem" parameter, and applies the
   *    contents of the solution vector to the primary variable field data,
   * 2) perform a synchronization of the primary field variable such that all ghosts are updated,
   *
   * The "scalingFactor" parameter allows for the scaled application of the solution vector. For
   * instance, a line search may apply a negative scaling factor to remove part of the previously
   * applied solution.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   *
   */
  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain );

  /**
   * @brief updates the configuration (if needed) based on the state after a converged Newton loop.
   * @param domain the domain containing the mesh and fields
   * @return a bool that states whether the configuration used to solve the nonlinear loop is still valid or not.
   */
  virtual bool updateConfiguration( DomainPartition & domain );

  /**
   * @brief
   * @param domain the domain containing the mesh and fields
   */
  virtual void outputConfigurationStatistics( DomainPartition const & domain ) const;

  /**
   * @brief resets the configuration to the beginning of the time-step.
   * @param domain the domain containing the mesh and fields
   */
  virtual void resetConfigurationToBeginningOfStep( DomainPartition & domain );

  /**
   * @brief resets the configuration to the default value.
   * @param domain the domain containing the mesh and fields
   * @return a bool that states whether the configuration was reset or not.
   */
  virtual bool resetConfigurationToDefault( DomainPartition & domain ) const;


  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  virtual void updateState( DomainPartition & domain );

  /**
   * @brief reset state of physics back to the beginning of the step.
   * @param domain
   *
   * This function essentially abandons the step, and resets all variables back to thier values at
   * the beginning of the step. This is useful for cases where convergence was not achieved, and
   * a cut in timestep was required.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain );

  /**
   * @brief perform cleanup for implicit timestep
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   *
   * This function performs whatever tasks are required to complete an implicit timestep. For
   * example, the acceptance of the solution will occur during this step, and deallocation of
   * temporaries will be be performed in this function.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain );


  /**
   * @brief getter for the next timestep size
   * @return the next timestep size m_nextDt
   */
  virtual real64 getTimestepRequest( real64 const ) override
  {return m_nextDt;};
  /**@}*/

  /**
   * @brief getter for the next timestep size
   * @return the next timestep size m_nextDt
   */
  real64 getTimestepRequest()
  {return m_nextDt;};

  /**
   * @brief creates a child group of of this PhysicsSolverBase instantiation
   * @param childKey the key of the child type
   * @param childName the name of the child
   * @return a pointer to the child group
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Type alias for catalog interface used by this class. See CatalogInterface.
   */
  using CatalogInterface = dataRepository::CatalogInterface< PhysicsSolverBase, string const &, Group * const >;

  /**
   * @brief Get the singleton catalog for PhysicsSolverBase.
   * @return reference to the catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Structure to hold scoped key names
   */
  struct viewKeyStruct
  {
    /// @return string for the cflFactor wrapper
    static constexpr char const * cflFactorString() { return "cflFactor"; }

    /// @return string for the initialDt wrapper
    static constexpr char const * initialDtString() { return "initialDt"; }

    /// @return string for the minDtIncreaseInterval wrapper
    static constexpr char const * minDtIncreaseIntervalString() { return "minDtIncreaseInterval"; }

    /// @return string for the maxStableDt wrapper
    static constexpr char const * maxStableDtString() { return "maxStableDt"; }

    /// @return string for the discretization wrapper
    static constexpr char const * discretizationString() { return "discretization"; }

    /// @return string for the nextDt targetRegions wrapper
    static constexpr char const * targetRegionsString() { return "targetRegions"; }

    /// @return string for the meshTargets wrapper
    static constexpr char const * meshTargetsString() { return "meshTargets"; }

    /// @return string for the writeLinearSystem wrapper
    static constexpr char const * writeLinearSystemString() { return "writeLinearSystem"; }
  };

  /**
   * @brief Structure to hold scoped key names
   */
  struct groupKeyStruct
  {
    /// @return string for the linearSolverParameters wrapper
    static constexpr char const * linearSolverParametersString() { return "LinearSolverParameters"; }

    /// @return string for the nonlinearSolverParameters wrapper
    static constexpr char const * nonlinearSolverParametersString() { return "NonlinearSolverParameters"; }

    /// @return string for the solverStatistics wrapper
    static constexpr char const * solverStatisticsString() { return "SolverStatistics"; }
  };

  /**
   * @brief getter for the timestamp of the system setup
   * @return the timestamp of the last time systemSetup was called
   */
  Timestamp getSystemSetupTimestamp() const { return m_systemSetupTimestamp; }

  /**
   * @brief getter for the timestamp of the mesh modification on the mesh levels
   * @param[in] domain the domain partition (cannot be const because we use forDiscretizationsInMeshTargets inside the function)
   * @return the timestamp of the last time at which one of the mesh levels was modified
   */
  Timestamp getMeshModificationTimestamp( DomainPartition & domain ) const;

  /**
   * @brief set the timestamp of the system setup
   * @param[in] timestamp the new timestamp of system setup
   */
  void setSystemSetupTimestamp( Timestamp timestamp ) { m_systemSetupTimestamp = timestamp; }

  /**
   * @brief return the value of the gravity vector specified in PhysicsSolverManager
   * @return the value of the gravity vector
   *
   * @note if the solver is instantiated outside of a simulation (for instance for a unit test)
   *       and therefore does not have a parent of type PhysicsSolverManager, this function returns
   *       {0.0,0.0,-9.81}
   */
  R1Tensor const gravityVector() const;

  /**
   * @brief Check if the solution increments are ok to use
   * @param domain the domain partition
   * @return true if the solution increments are ok to use, false otherwise
   */
  virtual bool checkSequentialSolutionIncrements( DomainPartition & domain ) const;

  /**
   * @brief Save the state of the solver for sequential iteration
   * @param domain the domain partition
   */
  virtual void saveSequentialIterationState( DomainPartition & domain );

  /**
   * @brief accessor for the linear solver parameters.
   * @return the linear solver parameter list
   */
  LinearSolverParameters & getLinearSolverParameters()
  {
    return m_linearSolverParameters.get();
  }

  /**
   * @brief const accessor for the linear solver parameters.
   * @return the linear solver parameter list
   */
  LinearSolverParameters const & getLinearSolverParameters() const
  {
    return m_linearSolverParameters.get();
  }

  /**
   * @brief accessor for the nonlinear solver parameters.
   * @return the nonlinear solver parameter list
   */
  NonlinearSolverParameters & getNonlinearSolverParameters()
  {
    return m_nonlinearSolverParameters;
  }

  /**
   * @brief const accessor for the nonlinear solver parameters.
   * @return the nonlinear solver parameter list
   */
  NonlinearSolverParameters const & getNonlinearSolverParameters() const
  {
    return m_nonlinearSolverParameters;
  }

  /**
   * @brief syncronize the nonlinear solver parameters.
   */
  virtual void
  synchronizeNonlinearSolverParameters()
  { /* empty here, overriden in CoupledSolver */ }

  /**
   * @brief Get position of a given region within solver's target region list
   * @param regionName the region name to find
   * @return index within target regions list
   */
  localIndex targetRegionIndex( string const & regionName ) const;



  /**
   * @brief Loop over the target discretization on all mesh targets and apply callback.
   * @tparam LAMBDA The callback function type
   * @param meshBodies The group of MeshBodies
   * @param lambda The callback function. Takes the name of the meshBody,
   * reference to the MeshLevel, and a list of regionNames.
   */
  template< typename LAMBDA >
  void forDiscretizationOnMeshTargets( Group const & meshBodies, LAMBDA && lambda ) const
  {
    for( auto const & target: m_meshTargets )
    {
      string const meshBodyName = target.first.first;
      string const meshLevelName = target.first.second;
      arrayView1d< string const > const & regionNames = target.second.toViewConst();
      MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );

      MeshLevel const * meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( meshLevelName );
      if( meshLevelPtr==nullptr )
      {
        meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( MeshBody::groupStructKeys::baseDiscretizationString() );
      }
      lambda( meshBodyName, *meshLevelPtr, regionNames );
    }
  }

  /**
   * @brief Loop over the target discretization on all mesh targets and apply callback.
   * @tparam LAMBDA The callback function type
   * @param meshBodies The group of MeshBodies
   * @param lambda The callback function. Takes the name of the meshBody,
   * reference to the MeshLevel, and a list of regionNames.
   */
  template< typename LAMBDA >
  void forDiscretizationOnMeshTargets( Group & meshBodies, LAMBDA && lambda ) const
  {
    for( auto const & target: m_meshTargets )
    {
      string const meshBodyName = target.first.first;
      string const meshLevelName = target.first.second;
      arrayView1d< string const > const & regionNames = target.second.toViewConst();
      MeshBody & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );

      MeshLevel * meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( meshLevelName );
      if( meshLevelPtr==nullptr )
      {
        meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( MeshBody::groupStructKeys::baseDiscretizationString() );
      }
      lambda( meshBodyName, *meshLevelPtr, regionNames );
    }
  }

  /**
   * @brief return the name of the discretization object
   * @return the name of the discretization object
   */
  string getDiscretizationName() const {return m_discretizationName;}

  /**
   * @brief function to set the value of m_assemblyCallback
   * @param func the function to set m_assemblyCallback to
   * @param funcType the type of the function
   * @return true if the function was successfully set, false otherwise
   *
   * This is used to provide a callback function for to be called in the assembly step.
   */
  virtual bool registerCallback( void * func, const std::type_info & funcType ) final override;

  /**
   * @brief accessor for the solver statistics.
   * @return reference to m_solverStatistics
   */
  SolverStatistics & getSolverStatistics() { return m_solverStatistics; }

  /**
   * @brief const accessor for the solver statistics.
   * @return reference to m_solverStatistics
   */
  SolverStatistics const & getSolverStatistics() const { return m_solverStatistics; }

#if defined(GEOS_USE_PYGEOSX)
  /**
   * @brief Return PySolver type.
   * @return Return PySolver type.
   */
  virtual PyTypeObject * getPythonType() const override;
#endif

  /**
   * @brief accessor for m_meshTargets
   * @return reference to m_meshTargets
   */
  map< std::pair< string, string >, array1d< string > > const & getMeshTargets() const
  {
    return m_meshTargets;
  }
protected:

  /**
   * @brief Eisenstat-Walker adaptive tolerance
   *
   * This method enables an inexact-Newton method is which the linear solver
   * tolerance is chosen based on the nonlinear solver convergence behavior.
   * In early Newton iterations, the search direction is usually imprecise, and
   * therefore a weak linear convergence tolerance can be chosen to minimize
   * computational cost.  As the search gets closer to the true solution, however,
   * more stringent linear tolerances are necessary to maintain quadratic convergence
   * behavior.
   *
   * The user can set the weakest tolerance allowed, with a default of 1e-3.
   * Even weaker values (e.g. 1e-2,1e-1) can be used for further speedup, but may
   * occasionally cause convergence problems.  Use this parameter with caution.  The
   * most stringent tolerance is hardcoded to 1e-8, which is sufficient for
   * most problems.
   *
   * See Eisenstat, S.C. and Walker, H.F., 1996. Choosing the forcing terms in an
   * inexact Newton method. SIAM Journal on Scientific Computing, 17(1), pp.16-32.
   *
   * @param newNewtonNorm Residual norm at current iteration
   * @param oldNewtonNorm Residual norm at previous iteration
   * @param krylovParams Linear solver parameters
   * @return Adaptive tolerance recommendation
   */
  real64 eisenstatWalker( real64 const newNewtonNorm,
                          real64 const oldNewtonNorm,
                          LinearSolverParameters::Krylov const & krylovParams );

  /**
   * @brief Get the Constitutive Name object
   *
   * @tparam CONSTITUTIVE_BASE_TYPE the base type of the constitutive model.
   * @param subRegion the element subregion on which the constitutive model is registered
   * @return the name name of the constitutive model of type CONSTITUTIVE_BASE_TYPE registered on the subregion.
   */
  template< typename CONSTITUTIVE_BASE_TYPE >
  static string getConstitutiveName( ElementSubRegionBase const & subRegion );

  /**
   * @brief Get the Constitutive Name object
   *
   * @tparam CONSTITUTIVE_BASE_TYPE the base type of the constitutive model.
   * @param subRegion the particle subregion on which the constitutive model is registered
   * @return the name name of the constitutive model of type CONSTITUTIVE_BASE_TYPE registered on the subregion.
   */
  template< typename CONSTITUTIVE_BASE_TYPE >
  static string getConstitutiveName( ParticleSubRegionBase const & subRegion ); // particle overload

  /**
   * @brief This function sets constitutive name fields on an
   *  ElementSubRegionBase, and calls the base function it overrides.
   * @param subRegion The ElementSubRegionBase that will have constitutive
   *  names set.
   */
  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const { GEOS_UNUSED_VAR( subRegion ); }


  /**
   * @brief Get the Constitutive Model object
   * @tparam BASETYPE the base type of the constitutive model.
   * @tparam LOOKUP_TYPE the type of the key used to look up the constitutive model.
   * @param dataGroup the data group containing the constitutive models.
   * @param key the key used to look up the constitutive model.
   * @return the constitutive model of type @p BASETYPE registered on the @p dataGroup with the key @p key.
   */
  template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
  static BASETYPE const & getConstitutiveModel( dataRepository::Group const & dataGroup, LOOKUP_TYPE const & key )
  {
    dataRepository::Group const & constitutiveModels = dataGroup.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
    return constitutiveModels.getGroup< BASETYPE >( key );
  }

  /**
   * @brief Get the Constitutive Model object
   * @tparam BASETYPE the base type of the constitutive model.
   * @tparam LOOKUP_TYPE the type of the key used to look up the constitutive model.
   * @param dataGroup the data group containing the constitutive models.
   * @param key the key used to look up the constitutive model.
   * @return the constitutive model of type @p BASETYPE registered on the @p dataGroup with the key @p key.
   */
  template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
  static BASETYPE & getConstitutiveModel( dataRepository::Group & dataGroup, LOOKUP_TYPE const & key )
  {
    dataRepository::Group & constitutiveModels = dataGroup.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
    return constitutiveModels.getGroup< BASETYPE >( key );
  }



  /// Courant–Friedrichs–Lewy factor for the timestep
  real64 m_cflFactor;

  /// maximum stable time step
  real64 m_maxStableDt;

  /// timestep of the next cycle
  real64 m_nextDt;

  /// Number of cycles since last timestep cut
  integer m_numTimestepsSinceLastDtCut;

  /// name of the FV discretization object in the data repository
  string m_discretizationName;

  /// Data structure to handle degrees of freedom
  DofManager m_dofManager;

  /// System matrix
  ParallelMatrix m_matrix;

  /// System right-hand side vector
  ParallelVector m_rhs;

  /// System solution vector
  ParallelVector m_solution;

  /// Local system matrix and rhs
  CRSMatrix< real64, globalIndex > m_localMatrix;

  /// Custom preconditioner for the "native" iterative solver
  std::unique_ptr< PreconditionerBase< LAInterface > > m_precond;

  /// flag for debug output of matrix, rhs, and solution
  integer m_writeLinearSystem;

  /// Linear solver parameters
  LinearSolverParametersInput m_linearSolverParameters;

  /// Result of the last linear solve
  LinearSolverResult m_linearSolverResult;

  /// Nonlinear solver parameters
  NonlinearSolverParameters m_nonlinearSolverParameters;

  /// Solver statistics
  SolverStatistics m_solverStatistics;

  /// Timestamp of the last call to setup system
  Timestamp m_systemSetupTimestamp;

  /// Callback function for assembly step
  std::function< void( CRSMatrix< real64, globalIndex >, array1d< real64 > ) > m_assemblyCallback;

  /// Timers for the aggregate profiling of the solver
  std::map< std::string, std::chrono::system_clock::duration > m_timers;

private:
  /// List of names of regions the solver will be applied to
  array1d< string > m_targetRegionNames;

  /// Map containing the array of target regions (value) for each MeshBody (key).
  map< std::pair< string, string >, array1d< string > > m_meshTargets;

  /**
   * @brief This function sets constitutive name fields on an
   *  ElementSubRegionBase, and DOES NOT call the base function it overrides.
   * @param subRegion The ElementSubRegionBase that will have constitutive
   *  names set.
   */
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const { GEOS_UNUSED_VAR( subRegion ); }

  /**
   * @brief Solve a nonlinear system using a Newton method
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain partition
   * @return true if the nonlinear system was solved, false otherwise
   */
  bool solveNonlinearSystem( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain );

  /**
   * @brief output information about the cycle to the log
   * @param cycleNumber the current cycle number
   * @param numOfSubSteps the number of substeps taken
   * @param subStepDt the time step size for each substep
   */
  void logEndOfCycleInformation( integer const cycleNumber,
                                 integer const numOfSubSteps,
                                 std::vector< real64 > const & subStepDt ) const;

};

template< typename CONSTITUTIVE_BASE_TYPE >
string PhysicsSolverBase::getConstitutiveName( ElementSubRegionBase const & subRegion )
{
  string validName;
  dataRepository::Group const & constitutiveModels = subRegion.getConstitutiveModels();

  constitutiveModels.forSubGroups< CONSTITUTIVE_BASE_TYPE >( [&]( dataRepository::Group const & model )
  {
    GEOS_ERROR_IF( !validName.empty(), "A valid constitutive model was already found." );
    validName = model.getName();
  } );
  return validName;
}

template< typename CONSTITUTIVE_BASE_TYPE >
string PhysicsSolverBase::getConstitutiveName( ParticleSubRegionBase const & subRegion ) // particle overload
{
  string validName;
  dataRepository::Group const & constitutiveModels = subRegion.getConstitutiveModels();

  constitutiveModels.forSubGroups< CONSTITUTIVE_BASE_TYPE >( [&]( dataRepository::Group const & model )
  {
    GEOS_ERROR_IF( !validName.empty(), "A valid constitutive model was already found." );
    validName = model.getName();
  } );
  return validName;
}


} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_PHYSICSSOLVERBASE_HPP_ */
