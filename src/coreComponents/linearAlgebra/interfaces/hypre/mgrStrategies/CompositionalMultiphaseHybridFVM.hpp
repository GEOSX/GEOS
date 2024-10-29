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
 * @file CompositionalMultiphaseHybridFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief
 *
 * Labels description stored in point_marker_array
 *                         0 = pressure
 *                         1 = density
 *                       ... = ... (densities)
 * numCellCenteredLabels - 1 = density
 * numLabels - 1             = face pressure
 *
 * 3-level MGR reduction strategy inspired from CompositionalMultiphaseReservoir:
 *   - 1st level: eliminate the density associated with the volume constraint
 *   - 2nd level: eliminate the rest of the densities
 *   - 3rd level: eliminate the cell-centered pressure
 *   - The coarse grid is the interface pressure system and is solved with BoomerAMG
 */
class CompositionalMultiphaseHybridFVM : public MGRStrategyBase< 3 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseHybridFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    // Level 0: eliminate the last density of the cell-centered block
    m_labels[0].resize( m_numBlocks - 1 );
    std::iota( m_labels[0].begin(), m_labels[0].begin() + m_numBlocks - 2, 0 );
    m_labels[0][m_numBlocks - 2] = m_numBlocks - 1;
    // Level 1: eliminate remaining densities of the cell-centered block
    m_labels[1].resize( 2 );
    m_labels[1][0] = 0;
    m_labels[1][1] = m_numBlocks - 1;
    // Level 2: eliminate reservoir cell-centered pressure
    m_labels[2].resize( 1 );
    m_labels[2][0] = m_numBlocks - 1;

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]          = MGRFRelaxationType::none;
    m_levelInterpType[1]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::blockColLumped; // True-IMPES
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;

    // Level 2
    m_levelFRelaxType[2]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[2]         = 1;
    m_levelInterpType[2]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[2]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::blockGaussSeidel;
    m_levelGlobalSmootherIters[2] = 1;
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
    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_*/
