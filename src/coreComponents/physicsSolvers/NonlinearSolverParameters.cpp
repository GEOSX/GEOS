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

#include "NonlinearSolverParameters.hpp"
#include "common/logger/Logger.hpp"
#include "common/format/table/TableFormatter.hpp"

namespace geos
{

using namespace dataRepository;

NonlinearSolverParameters::NonlinearSolverParameters( string const & name,
                                                      Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  // This enables logLevel filtering
  enableLogLevelInput();

  registerWrapper( viewKeysStruct::lineSearchActionString(), &m_lineSearchAction ).
    setApplyDefaultValue( LineSearchAction::Attempt ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "How the line search is to be used. Options are: \n "
                    "* None    - Do not use line search.\n"
                    "* Attempt - Use line search. Allow exit from line search without achieving smaller residual than starting residual.\n"
                    "* Require - Use line search. If smaller residual than starting resdual is not achieved, cut time-step." );

  registerWrapper( viewKeysStruct::lineSearchInterpolationTypeString(), &m_lineSearchInterpType ).
    setApplyDefaultValue( LineSearchInterpolationType::Linear ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strategy to cut the solution update during the line search. Options are: \n "
                    "* Linear\n"
                    "* Parabolic" );

  registerWrapper( viewKeysStruct::lineSearchMaxCutsString(), &m_lineSearchMaxCuts ).
    setApplyDefaultValue( 4 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of line search cuts." );

  registerWrapper( viewKeysStruct::lineSearchCutFactorString(), &m_lineSearchCutFactor ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Line search cut factor. For instance, a value of 0.5 will result in the effective application of"
                    " the last solution by a factor of (0.5, 0.25, 0.125, ...)" );

  registerWrapper( viewKeysStruct::lineSearchStartingIterationString(), &m_lineSearchStartingIteration ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Iteration when line search starts." );

  registerWrapper( viewKeysStruct::lineSearchResidualFactorString(), &m_lineSearchResidualFactor ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor to determine residual increase (recommended values: 1.1 (conservative), 2.0 (relaxed), 10.0 (aggressive))." );

  registerWrapper( viewKeysStruct::normTypeString(), &m_normType ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( physicsSolverBaseKernels::NormType::Linf ).
    setDescription( "Norm used by the flow solver to check nonlinear convergence. "
                    "Valid options:\n* " + EnumStrings< physicsSolverBaseKernels::NormType >::concat( "\n* " ) );

  registerWrapper( viewKeysStruct::minNormalizerString(), &m_minNormalizer ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e-12 ).
    setDescription( "Value used to make sure that residual normalizers are not too small when computing residual norm." );

  registerWrapper( viewKeysStruct::newtonTolString(), &m_newtonTol ).
    setApplyDefaultValue( 1.0e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The required tolerance in order to exit the Newton iteration loop." );

  registerWrapper( viewKeysStruct::newtonMaxIterString(), &m_maxIterNewton ).
    setApplyDefaultValue( 5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of iterations that are allowed in a Newton loop." );

  registerWrapper( viewKeysStruct::newtonMinIterString(), &m_minIterNewton ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum number of iterations that are required before exiting the Newton loop." );

  registerWrapper( viewKeysStruct::newtonNumIterationsString(), &m_numNewtonIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of Newton's iterations." );

  registerWrapper( viewKeysStruct::maxAllowedResidualNormString(), &m_maxAllowedResidualNorm ).
    setApplyDefaultValue( 1e9 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum value of residual norm that is allowed in a Newton loop" );


  registerWrapper( viewKeysStruct::allowNonConvergedString(), &m_allowNonConverged ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Allow non-converged solution to be accepted. "
                    "(i.e. exit from the Newton loop without achieving the desired tolerance)" );

  registerWrapper( viewKeysStruct::timeStepDecreaseIterLimString(), &m_timeStepDecreaseIterLimit ).
    setApplyDefaultValue( 0.7 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fraction of the max Newton iterations above which the solver asks for the time-step to be decreased for the next time-step." );

  registerWrapper( viewKeysStruct::timeStepIncreaseIterLimString(), &m_timeStepIncreaseIterLimit ).
    setApplyDefaultValue( 0.4 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fraction of the max Newton iterations below which the solver asks for the time-step to be increased for the next time-step." );

  registerWrapper( viewKeysStruct::timeStepDecreaseFactorString(), &m_timeStepDecreaseFactor ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which the time-step is decreased when the number of Newton iterations is large." );

  registerWrapper( viewKeysStruct::timeStepIncreaseFactorString(), &m_timeStepIncreaseFactor ).
    setApplyDefaultValue( 2.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which the time-step is increased when the number of Newton iterations is small." );

  registerWrapper( viewKeysStruct::minTimeStepIncreaseIntervalString(), &m_minTimeStepIncreaseInterval ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum number of cycles since the last time-step cut for increasing the time-step again." );

  registerWrapper( viewKeysStruct::timeStepCutFactorString(), &m_timeStepCutFactor ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which the time-step will be cut if a time-step cut is required." );

  registerWrapper( viewKeysStruct::maxTimeStepCutsString(), &m_maxTimeStepCuts ).
    setApplyDefaultValue( 2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max number of time-step cuts" );

  registerWrapper( viewKeysStruct::maxSubStepsString(), &m_maxSubSteps ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of time sub-steps allowed for the solver" );

  registerWrapper( viewKeysStruct::maxNumConfigurationAttemptsString(), &m_maxNumConfigurationAttempts ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max number of times that the configuration can be changed" );

  registerWrapper( viewKeysStruct::configurationToleranceString(), &m_configurationTolerance ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Configuration tolerance" );

  /// GEOS mainly uses FIM coupling so let's define FIM as the default.
  registerWrapper( viewKeysStruct::couplingTypeString(), &m_couplingType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( CouplingType::FullyImplicit ).
    setDescription( "Type of coupling. "
                    "Valid options:\n* " + EnumStrings< CouplingType >::concat( "\n* " ) );

  registerWrapper( viewKeysStruct::sequentialConvergenceCriterionString(), &m_sequentialConvergenceCriterion ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( SequentialConvergenceCriterion::ResidualNorm ).
    setDescription( "Criterion used to check outer-loop convergence in sequential schemes. "
                    "Valid options:\n* " + EnumStrings< SequentialConvergenceCriterion >::concat( "\n* " ) );

  registerWrapper( viewKeysStruct::subcyclingOptionString(), &m_subcyclingOption ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to decide whether to iterate between sequentially coupled solvers or not." );

  registerWrapper( viewKeysStruct::nonlinearAccelerationTypeString(), &m_nonlinearAccelerationType ).
    setApplyDefaultValue( NonlinearAccelerationType::None ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Nonlinear acceleration type for sequential solver." );

}

void NonlinearSolverParameters::postInputInitialization()
{
  GEOS_ERROR_IF_LE_MSG( m_timeStepDecreaseIterLimit, m_timeStepIncreaseIterLimit,
                        getWrapperDataContext( viewKeysStruct::timeStepIncreaseIterLimString() ) <<
                        ": should be smaller than " << viewKeysStruct::timeStepDecreaseIterLimString() );

  GEOS_ERROR_IF_LE_MSG( m_lineSearchResidualFactor, 0.0,
                        getWrapperDataContext( viewKeysStruct::lineSearchResidualFactorString() ) << ": should be positive" );

  if( getLogLevel() > 0 )
  {
    print();
  }
}

void NonlinearSolverParameters::print() const
{
  TableData tableData;
  tableData.addRow( "Log level", getLogLevel());
  tableData.addRow( "Line search", "" );
  tableData.addRow( "  Action", m_lineSearchAction );
  if( m_lineSearchAction != LineSearchAction::None )
  {
    tableData.addRow( "  Interpolation type", m_lineSearchInterpType );
    tableData.addRow( "  Maximum number of cuts", m_lineSearchMaxCuts );
    tableData.addRow( "  Cut factor", m_lineSearchCutFactor );
    tableData.addRow( "  Starting iteration", m_lineSearchStartingIteration );
    tableData.addRow( "  Residual increase factor", m_lineSearchResidualFactor );
  }
  tableData.addRow( "Norm type (flow solver)", m_normType );
  tableData.addRow( "Minimum residual normalizer", m_minNormalizer );
  tableData.addRow( "Convergence tolerance", m_newtonTol );
  tableData.addRow( "Maximum iterations", m_maxIterNewton );
  tableData.addRow( "Minimum iterations", m_minIterNewton );
  tableData.addRow( "Maximum allowed residual norm", m_maxAllowedResidualNorm );
  tableData.addRow( "Allow non-converged", m_allowNonConverged );
  tableData.addRow( "Time-step decrease iterations limit", m_timeStepDecreaseIterLimit );
  tableData.addRow( "Time-step increase iterations limit", m_timeStepIncreaseIterLimit );
  tableData.addRow( "Time-step decrease factor", m_timeStepDecreaseFactor );
  tableData.addRow( "Time-step increase factor", m_timeStepDecreaseFactor );
  tableData.addRow( "Time-step cut factor", m_timeStepCutFactor );
  tableData.addRow( "Minimum time-step increase interval", m_minTimeStepIncreaseInterval );
  tableData.addRow( "Maximum time-step cuts", m_maxTimeStepCuts );
  tableData.addRow( "Maximum sub time-steps", m_maxSubSteps );
  tableData.addRow( "Maximum number of configuration attempts", m_maxNumConfigurationAttempts );
  tableData.addRow( "Coupling type", m_couplingType );
  if( m_couplingType == CouplingType::Sequential )
  {
    tableData.addRow( "Sequential convergence criterion", m_sequentialConvergenceCriterion );
    tableData.addRow( "Subcycling", m_subcyclingOption );
  }
  TableLayout const tableLayout = TableLayout( GEOS_FMT( "{}: nonlinear solver", getParent().getName() ),
                                               { TableLayout::Column{"Parameter", TableLayout::Alignment::left},
                                                 TableLayout::Column{"Value", TableLayout::Alignment::left} }
                                               );
  TableTextFormatter const tableFormatter( tableLayout );
  GEOS_LOG_RANK_0( tableFormatter.toString( tableData ));
}

REGISTER_CATALOG_ENTRY( Group, NonlinearSolverParameters, string const &, Group * const )

} /* namespace geos */
