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

#ifndef GEOS_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_
#define GEOS_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "physicsSolvers/SolverBaseKernels.hpp"

namespace geos
{

/**
 * @class Holds parameters and status for execution of nonlinear solution schemes.
 */
class NonlinearSolverParameters : public dataRepository::Group
{
public:

  /**
   * @brief Default constructor.
   */
  NonlinearSolverParameters() = delete;

  /**
   * @brief Constructor
   * @param[in] name The name of the new instantiation of this Group.
   * @param[in] parent A pointer to the parent of this Group.
   */
  NonlinearSolverParameters( string const & name,
                             Group * const parent );

  /**
   * @brief Default destructor
   */
  virtual ~NonlinearSolverParameters() = default;

  /**
   * @brief Default Move Constructor
   * @param The source object of the move.
   */
  NonlinearSolverParameters( NonlinearSolverParameters && ) = default;

  /**
   * @brief Copy Constructor
   * @param The source object.
   */
  NonlinearSolverParameters & operator=( const NonlinearSolverParameters & params )
  {
    m_lineSearchAction = params.m_lineSearchAction;
    m_lineSearchInterpType = params.m_lineSearchInterpType;
    m_lineSearchMaxCuts = params.m_lineSearchMaxCuts;
    m_lineSearchCutFactor = params.m_lineSearchCutFactor;
    m_lineSearchStartingIteration = params.m_lineSearchStartingIteration;
    m_lineSearchResidualFactor = params.m_lineSearchResidualFactor;

    m_newtonTol = params.m_newtonTol;
    m_maxIterNewton = params.m_maxIterNewton;
    m_minIterNewton = params.m_minIterNewton;
    m_numNewtonIterations = params.m_numNewtonIterations;

    m_maxAllowedResidualNorm = params.m_maxAllowedResidualNorm;
    m_allowNonConverged = params.m_allowNonConverged;

    m_timeStepDecreaseIterLimit = params.m_timeStepDecreaseIterLimit;
    m_timeStepIncreaseIterLimit = params.m_timeStepIncreaseIterLimit;
    m_timeStepDecreaseFactor = params.m_timeStepDecreaseFactor;
    m_timeStepIncreaseFactor = params.m_timeStepIncreaseFactor;
    m_minTimeStepIncreaseInterval = params.m_minTimeStepIncreaseInterval;
    m_maxSubSteps = params.m_maxSubSteps;
    m_maxTimeStepCuts = params.m_maxTimeStepCuts;
    m_timeStepCutFactor = params.m_timeStepCutFactor;
    m_maxNumConfigurationAttempts = params.m_maxNumConfigurationAttempts;
    m_configurationTolerance = params.m_configurationTolerance;

    setLogLevel( params.getLogLevel());

    return *this;
  }

  /**
   * @brief The name of this object in the catalog.
   * @return A string containing the name of this object in the catalog.
   * The CatalogName is the string that will result in the creation of a new
   * NonlinearSolverParameters object when calling
   * Group::getCatalog()::Allocate().
   */
  static string catalogName() { return "NonlinearSolverParameters"; }

  virtual void postInputInitialization() override;

  void print() const;

  struct viewKeysStruct
  {
    static constexpr char const * lineSearchActionString()        { return "lineSearchAction"; }
    static constexpr char const * lineSearchMaxCutsString()       { return "lineSearchMaxCuts"; }
    static constexpr char const * lineSearchCutFactorString()     { return "lineSearchCutFactor"; }
    static constexpr char const * lineSearchInterpolationTypeString() { return "lineSearchInterpolationType"; }
    static constexpr char const * lineSearchStartingIterationString() { return "lineSearchStartingIteration"; }
    static constexpr char const * lineSearchResidualFactorString() { return "lineSearchResidualFactor"; }

    static constexpr char const * normTypeString()                { return "normType"; }
    static constexpr char const * minNormalizerString()           { return "minNormalizer"; }
    static constexpr char const * newtonTolString()               { return "newtonTol"; }
    static constexpr char const * newtonMaxIterString()           { return "newtonMaxIter"; }
    static constexpr char const * newtonMinIterString()           { return "newtonMinIter"; }
    static constexpr char const * newtonNumIterationsString()     { return "newtonNumberOfIterations"; }
    static constexpr char const * newtonSplitOperMaxIterString()  { return "newtonSplitOperMaxIter"; }

    static constexpr char const * allowNonConvergedString()       { return "allowNonConverged"; }
    static constexpr char const * timeStepDecreaseIterLimString() { return "timeStepDecreaseIterLimit"; }
    static constexpr char const * timeStepIncreaseIterLimString() { return "timeStepIncreaseIterLimit"; }
    static constexpr char const * timeStepDecreaseFactorString()  { return "timeStepDecreaseFactor"; }
    static constexpr char const * timeStepIncreaseFactorString()  { return "timeStepIncreaseFactor"; }
    static constexpr char const * minTimeStepIncreaseIntervalString()  { return "minTimeStepIncreaseInterval"; }

    static constexpr char const * maxSubStepsString()             { return "maxSubSteps"; }
    static constexpr char const * maxTimeStepCutsString()         { return "maxTimeStepCuts"; }
    static constexpr char const * timeStepCutFactorString()       { return "timeStepCutFactor"; }
    static constexpr char const * maxAllowedResidualNormString()  { return "maxAllowedResidualNorm"; }

    static constexpr char const * maxNumConfigurationAttemptsString() { return "maxNumConfigurationAttempts"; }
    static constexpr char const * configurationToleranceString() { return "configurationTolerance"; }

    static constexpr char const * couplingTypeString()                   { return "couplingType"; }
    static constexpr char const * sequentialConvergenceCriterionString() { return "sequentialConvergenceCriterion"; }
    static constexpr char const * subcyclingOptionString()               { return "subcycling"; }
    static constexpr char const * nonlinearAccelerationTypeString() { return "nonlinearAccelerationType"; }
  } viewKeys;

  /**
   * @brief Indicates the handling of line search in a Newton loop.
   */
  enum class LineSearchAction : integer
  {
    None,    ///< Do not use line search
    Attempt, ///< Use line search. Allow exit from line search without achieving smaller residual than starting residual.
    Require, ///< Use line search. If smaller residual than starting residual is not achieved, cut time step.
  };

  /**
   * @brief Indicates the handling of line each interpolation strategy.
   */
  enum class LineSearchInterpolationType : integer
  {
    Linear,    ///< linear decrease of line search scaling factor.
    Parabolic, ///< use parabolic interpolation to define line search scaling factor.
  };

  /**
   * @brief Coupling type.
   */
  enum class CouplingType : integer
  {
    FullyImplicit,      ///< Fully-implicit coupling
    Sequential ///< Sequential coupling
  };

  /**
   * @brief Sequential convergence criterion
   */
  enum class SequentialConvergenceCriterion : integer
  {
    ResidualNorm, ///< convergence achieved when the residual drops below a given norm
    NumberOfNonlinearIterations, ///< convergence achieved when the subproblems convergence is achieved in less than minNewtonIteration
    SolutionIncrements ///< convergence achieved when the solution increments are small enough
  };

  /**
   * @brief Nonlinear acceleration type
   */
  enum class NonlinearAccelerationType : integer
  {
    None, ///< no acceleration
    Aitken ///< Aitken acceleration
  };

  /**
   * @brief Calculates the upper limit for the number of iterations to allow a
   * decrease to the next time step.
   * @return The scaled value of the limit (m_timeStepDecreaseIterLimit * m_maxIterNewton)
   */
  integer timeStepDecreaseIterLimit() const
  {
    return std::ceil( m_timeStepDecreaseIterLimit * m_maxIterNewton );
  }


  /**
   * @brief Calculates the lower limit for the number of iterations to force an
   * increase to the next time step.
   * @return The scaled value of the limit (m_timeStepIncreaseIterLimit * m_maxIterNewton)
   */
  integer timeStepIncreaseIterLimit() const
  {
    return std::ceil( m_timeStepIncreaseIterLimit * m_maxIterNewton );
  }

  /**
   * @brief Getter for the factor used to decrease the time step size
   * @return the factor used to decrease the time step size
   */
  real64 timeStepDecreaseFactor() const
  {
    return m_timeStepDecreaseFactor;
  }

  /**
   * @brief Getter for the factor used to increase the time step size
   * @return the factor used to increase the time step size
   */
  real64 timeStepIncreaseFactor() const
  {
    return m_timeStepIncreaseFactor;
  }

  /**
   * @brief Getter for the minimum interval for increasing the time-step
   * @return the minimum interval for increasing the time-step
   */
  integer minTimeStepIncreaseInterval() const
  {
    return m_minTimeStepIncreaseInterval;
  }

  /**
   * @brief Getter for the norm type used to check convergence in the flow/well solvers
   * @return the norm type
   */
  solverBaseKernels::NormType normType() const
  {
    return m_normType;
  }

  /**
   * @brief Getter for the coupling type
   * @return the coupling type
   */
  CouplingType couplingType() const
  {
    return m_couplingType;
  }

  /**
   * @brief Getter for the sequential convergence criterion
   * @return the sequential convergence criterion
   */
  SequentialConvergenceCriterion sequentialConvergenceCriterion() const
  {
    return m_sequentialConvergenceCriterion;
  }

  /// Flag to apply a line search.
  LineSearchAction m_lineSearchAction;

  /// Flag to pick the type of line search
  LineSearchInterpolationType m_lineSearchInterpType;

  /// The maximum number of line search cuts to attempt.
  integer m_lineSearchMaxCuts;

  /// The reduction factor for each line search cut.
  real64 m_lineSearchCutFactor;

  /// Iteration when line search starts
  integer m_lineSearchStartingIteration;

  /// Factor to determine residual increase
  real64 m_lineSearchResidualFactor;

  /// Norm used to check the nonlinear loop convergence
  solverBaseKernels::NormType m_normType;

  /// The tolerance for the nonlinear convergence check.
  real64 m_newtonTol;

  /// The maximum number of nonlinear iterations that are allowed.
  integer m_maxIterNewton;

  /// The minimum number of nonlinear iterations that may be applied.
  integer m_minIterNewton;

  /// The number of nonlinear iterations that have been exectued.
  integer m_numNewtonIterations;

  /// The maximum value of residual norm that we allow (otherwise, we cut the time step)
  real64 m_maxAllowedResidualNorm;

  /// Flag to allow for a non-converged nonlinear solution and continue with the problem.
  integer m_allowNonConverged;

  /// Fraction of the max Newton iterations above which the solver asks for the time-step to be decreased for the next timeStep
  real64 m_timeStepDecreaseIterLimit;

  /// Fraction of the max Newton iterations below which the solver asks for the time-step to be increased for the next timeStep
  real64 m_timeStepIncreaseIterLimit;

  /// Factor used to decrease the time step size
  real64 m_timeStepDecreaseFactor;

  /// Factor used to increase the time step size
  real64 m_timeStepIncreaseFactor;

  /// Minimum interval, since the last time-step cut, for increasing the time-step
  integer m_minTimeStepIncreaseInterval;

  /// Maximum number of time sub-steps allowed for the solver
  integer m_maxSubSteps;

  /// Max number of time step cuts
  integer m_maxTimeStepCuts;

  /// Factor by which the time step will be cut if a timestep cut is required.
  real64 m_timeStepCutFactor;

  /// Number of times that the time-step had to be cut
  integer m_numTimeStepAttempts;

  /// Number of times that the configuration had to be changed
  integer m_numConfigurationAttempts;

  /// Max number of times that the configuration can be changed
  integer m_maxNumConfigurationAttempts;

  /// Configuration tolerance
  double m_configurationTolerance;

  /// Type of coupling
  CouplingType m_couplingType;

  /// Criterion used to check outer-loop convergence in sequential schemes
  SequentialConvergenceCriterion m_sequentialConvergenceCriterion;

  /// Flag to specify whether subcycling is allowed or not in sequential schemes
  integer m_subcyclingOption;

  /// Type of nonlinear acceleration for sequential solver
  NonlinearAccelerationType m_nonlinearAccelerationType;

  /// Value used to make sure that residual normalizers are not too small when computing residual norm
  real64 m_minNormalizer = 1e-12;
};

ENUM_STRINGS( NonlinearSolverParameters::LineSearchAction,
              "None",
              "Attempt",
              "Require" );

ENUM_STRINGS( NonlinearSolverParameters::LineSearchInterpolationType,
              "Linear",
              "Parabolic" );

ENUM_STRINGS( NonlinearSolverParameters::CouplingType,
              "FullyImplicit",
              "Sequential" );

ENUM_STRINGS( NonlinearSolverParameters::SequentialConvergenceCriterion,
              "ResidualNorm",
              "NumberOfNonlinearIterations",
              "SolutionIncrements" );

ENUM_STRINGS( NonlinearSolverParameters::NonlinearAccelerationType,
              "None",
              "Aitken" );

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_ */
