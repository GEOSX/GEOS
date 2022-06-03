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

/**
 * @file CompositionalMultiphaseFVM.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief CompositionalMultiphaseFVM strategy.
 *
 * Labels description stored in point_marker_array
 *               0 = pressure
 *               1 = density
 *             ... = densities
 *   numLabels - 1 = density
 *
 * 2-level MGR reduction strategy which seems to work well for 2 components:
 *   - 1st level: eliminate the reservoir density associated with the volume constraint
 *   - 2nd level: eliminate the other reservoir density
 *   - The coarse grid (pressure system) is solved with BoomerAMG.
 *
 * @todo:
 *   - Experiment with block Jacobi for F-relaxation/interpolation of the reservoir densities
 */
class CompositionalMultiphaseFVM : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] ) )
  {
    // Level 0: eliminate last density which corresponds to the volume constraint equation
    m_labels[0].resize( m_numBlocks - 1 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    // Level 1: eliminate the other density
    m_labels[1].push_back( 0 );

    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[0]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::galerkin;

    // Level 1
    m_levelFRelaxMethod[1]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[1]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = MGRCoarseGridMethod::galerkin;

    m_globalSmoothType = MGRGlobalSmootherType::ilu0; // TODO switch to levelGlobalSmoother when hypre allows
    m_numGlobalSmoothSweeps = 0; // No global smoother
  }

  /**
   * @brief Setup the MGR strategy.
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const &,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                  m_numBlocks, numLevels,
                                                                  m_numLabels, m_ptrLabels,
                                                                  mgrData.pointMarkers.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, toUnderlyingPtr( m_levelFRelaxMethod ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalSmoothType( precond.ptr, toUnderlying( m_globalSmoothType ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalSmoothIters( precond.ptr, m_numGlobalSmoothSweeps ) );
#ifdef GEOSX_USE_HYPRE_CUDA // GEOSX_USE_HYPRE_ROCM
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( precond.ptr, 1e-20 )); // truncate intermediate/coarse grids
#endif

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
#ifdef GEOSX_USE_HYPRE_CUDA // GEOSX_USE_HYPRE_ROCM
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.coarseSolver.ptr, toUnderlying( AMGCoarseningType::PMIS ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.coarseSolver.ptr, 2 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.coarseSolver.ptr, 1.0 ) );
#else
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ) );
#endif

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_*/
