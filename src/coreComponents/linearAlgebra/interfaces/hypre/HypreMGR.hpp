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
 * @file HypreMGR.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_

#include "linearAlgebra/common/common.hpp"

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <_hypre_utilities.h>

namespace geos
{

/**
 * @brief Container for hypre preconditioner auxiliary data for MGR.
 */
struct HypreMGRData
{
  array1d< HYPRE_Int > pointMarkers;  ///< array1d of unique tags for local degrees of freedom
  HyprePrecWrapper coarseSolver;      ///< MGR coarse solver pointer and functions
  HyprePrecWrapper mechSolver;        ///< MGR mechanics fine solver pointer and functions
};

namespace hypre
{

namespace mgr
{

/**
 * @brief Helper to simplify MGR setup
 * @tparam STRATEGY strategy class (one of structs defined below)
 * @param params MGR parameters
 * @param numComponentsPerField array containing number of components in each field
 * @param precond the preconditioner wrapper
 * @param mgrData additional MGR data struct with marker array populated
 */
template< typename STRATEGY >
void setStrategy( LinearSolverParameters::MGR const & params,
                  arrayView1d< int const > const & numComponentsPerField,
                  HyprePrecWrapper & precond,
                  HypreMGRData & mgrData )
{
  STRATEGY strategy( numComponentsPerField );
  strategy.setup( params, precond, mgrData );
}

/**
 * @brief Helper struct for strategies that provides some basic parameter arrays needed by MGR.
 * @tparam NLEVEL number of reduction levels (not including the coarsest level)
 */
template< int NLEVEL >
class MGRStrategyBase
{
public:

  static constexpr HYPRE_Int numLevels = NLEVEL;              ///< Number of levels

protected:

  HYPRE_Int m_numBlocks{ 0 };                                 ///< Number of different matrix blocks treated separately

  std::vector< HYPRE_Int > m_labels[numLevels]{};             ///< Dof labels kept at each level
  HYPRE_Int m_numLabels[numLevels]{ -1 };                     ///< Number of dof labels kept
  HYPRE_Int * m_ptrLabels[numLevels]{ nullptr };              ///< Pointers to each level's labels, as consumed by MGR

  MGRFRelaxationType m_levelFRelaxType[numLevels];            ///< F-relaxation type for each level
  HYPRE_Int m_levelFRelaxIters[numLevels]{ -1 };              ///< Number of F-relaxation iterations for each level
  MGRInterpolationType m_levelInterpType[numLevels];          ///< Interpolation type for each level
  MGRRestrictionType m_levelRestrictType[numLevels];          ///< Restriction type for each level
  MGRCoarseGridMethod m_levelCoarseGridMethod[numLevels];     ///< Coarse grid method for each level
  MGRGlobalSmootherType m_levelGlobalSmootherType[numLevels]; ///< Global smoother type for each level
  HYPRE_Int m_levelGlobalSmootherIters[numLevels]{ -1 };      ///< Number of global smoother iterations for each level
  HYPRE_Real m_coarseGridThreshold{ 1.0e-20 };                ///< Coarse grid truncation threshold

  // TODO: the following options are currently commented out in MGR's code.
  //       Let's consider their use when re-enable in hypre
  // HYPRE_Int m_numRestrictSweeps{ -1 }; ///< Number of restrict sweeps
  // HYPRE_Int m_numInterpSweeps{ -1 };   ///< Number of interpolation sweeps

  /**
   * @brief Constructor.
   * @param numBlocks number of blocks
   */
  explicit MGRStrategyBase( HYPRE_Int const numBlocks )
    : m_numBlocks( numBlocks )
  {
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      m_levelFRelaxType[i]         = MGRFRelaxationType::jacobi;
      m_levelInterpType[i]         = MGRInterpolationType::jacobi;
      m_levelRestrictType[i]       = MGRRestrictionType::injection;
      m_levelCoarseGridMethod[i]   = MGRCoarseGridMethod::galerkin;
      m_levelGlobalSmootherType[i] = MGRGlobalSmootherType::none;
    }
  }

  /**
   * @brief Call this after populating lv_cindexes.
   */
  void setupLabels()
  {
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      m_numLabels[i] = m_labels[i].size();
      m_ptrLabels[i] = m_labels[i].data();
    }
  }

  /**
   * @brief Helper function that sets the reduction features common to all mgr strategies
   * @param precond the preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setReduction( HyprePrecWrapper & precond,
                     HypreMGRData & mgrData )

  {
    // Ensure that if no F-relaxation or global smoothing is chosen the corresponding number
    // of iteration is set to 0
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      if( m_levelFRelaxType[i] == MGRFRelaxationType::none )
      {
        m_levelFRelaxIters[i] = 0;
      }
      if( m_levelGlobalSmootherType[i] == MGRGlobalSmootherType::none )
      {
        m_levelGlobalSmootherIters[i] = 0;
      }
    }

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                 m_numBlocks, numLevels,
                                                                 m_numLabels, m_ptrLabels,
                                                                 mgrData.pointMarkers.data() ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxType( precond.ptr, toUnderlyingPtr( m_levelFRelaxType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelNumRelaxSweeps( precond.ptr, m_levelFRelaxIters ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( precond.ptr, toUnderlyingPtr( m_levelGlobalSmootherType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( precond.ptr, m_levelGlobalSmootherIters ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( precond.ptr, m_coarseGridThreshold ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonGalerkinMaxElmts( precond.ptr, 1 ));
  }

  /**
   * @brief Set up BoomerAMG to perform the solve for the displacement system
   * @param solver the solver wrapper
   * @param separateComponents flag controlling the use of the separate displacement component (SDC) approximation
   */
  void setDisplacementAMG( HyprePrecWrapper & solver,
                           integer const & separateComponents )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &solver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( solver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( solver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( solver.ptr, 1.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( solver.ptr, 0.6 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( solver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpType( solver.ptr, hypre::getAMGInterpolationType( LinearSolverParameters::AMG::InterpType::modifiedExtendedE )) );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( solver.ptr, hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( solver.ptr, hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::chebyshev ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( solver.ptr, 1 ) );
#else
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( solver.ptr, 1 ) );
#endif

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( solver.ptr, 3 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetFilterFunctions( solver.ptr, separateComponents ) );

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up BoomerAMG to perform the solve for the pressure system
   * @param solver the solver wrapper
   */
  void setPressureAMG( HyprePrecWrapper & solver )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &solver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( solver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( solver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( solver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggPMaxElmts( solver.ptr, 16 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggInterpType( solver.ptr, hypre::getAMGAggressiveInterpolationType( LinearSolverParameters::AMG::AggInterpType::multipass ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( solver.ptr, 0.0 ) );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( solver.ptr, hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( solver.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( solver.ptr, 2 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( solver.ptr, 1.0 ) );
#else
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( solver.ptr, 1 ) );
#endif

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up BoomerAMG to perform the solve for the pressure/temperature system
   * @param solver the solver wrapper
   */
  void setPressureTemperatureAMG( HyprePrecWrapper & solver )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &solver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( solver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( solver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( solver.ptr, 1 ) ); // TODO: keep or not 1 aggressive level?
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggPMaxElmts( solver.ptr, 16 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( solver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( solver.ptr, 2 ) ); // pressure and temperature (CPTR)
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( solver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( solver.ptr, toUnderlying( AMGCoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( solver.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( solver.ptr, 2 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( solver.ptr, 1.0 ) );
#else
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( solver.ptr, 1 ) );
#endif

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up BoomerAMG to perform the mechanics F-solve for the first F-relaxation
   * @param precond the preconditioner wrapper
   * @param mgrData auxiliary MGR data
   * @param separateComponents flag controlling the use of the separate displacement component (SDC) approximation
   *
   * @note This function should be rethought once MGR allows for customizing boomerAMG (or
   *       any other solver) for F-relaxation at any level
   */
  void setMechanicsFSolver( HyprePrecWrapper & precond,
                            HypreMGRData & mgrData,
                            integer const & separateComponents )
  {
    setDisplacementAMG( mgrData.mechSolver, separateComponents );
    HYPRE_MGRSetFSolver( precond.ptr, mgrData.mechSolver.solve, mgrData.mechSolver.setup, mgrData.mechSolver.ptr );
  }

};

/**
 * @brief Create the MGR preconditioner object.
 * @param params preconditioner parameters
 * @param dofManager pointer to DofManager for the linear system
 * @param precond the preconditioner
 * @param mgrData auxiliary data for MGR
 */
void createMGR( LinearSolverParameters const & params,
                DofManager const * const dofManager,
                HyprePrecWrapper & precond,
                HypreMGRData & mgrData );

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_*/
