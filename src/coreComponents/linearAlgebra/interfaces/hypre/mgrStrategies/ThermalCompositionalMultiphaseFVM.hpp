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
 * @file ThermalCompositionalMultiphaseFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALCOMPOSITIONALMULTIPHASEFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALCOMPOSITIONALMULTIPHASEFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief ThermalCompositionalMultiphaseFVM strategy.
 *
 * Labels description stored in point_marker_array
 *               0 = pressure
 *               1 = density
 *             ... = densities
 *   numLabels - 2 = last density
 *   numLabels - 1 = temperature
 *
 * 2-level MGR reduction strategy
 *   - 1st level: eliminate the reservoir density associated with the volume constraint
 *   - 2nd level: eliminate the other reservoir densities
 *   - The coarse grid (pressure and temperature system) is solved with BoomerAMG with numFunctions==2.
 *
 */
class ThermalCompositionalMultiphaseFVM : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit ThermalCompositionalMultiphaseFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] ) )
  {
    // Level 0: eliminate last density which corresponds to the volume constraint equation
    m_labels[0].resize( m_numBlocks - 2 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    m_labels[0].push_back( m_numBlocks-1 ); // keep temperature
    // Level 1: eliminate the other densities
    m_labels[1].push_back( 0 ); // keep pressure
    m_labels[1].push_back( m_numBlocks-1 ); // keep temperature

    setupLabels();

    m_levelFRelaxType[0]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::jacobi; // Diagonal scaling (Jacobi)
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::blockGaussSeidel;
    m_levelGlobalSmootherIters[0] = 1;

    m_levelFRelaxType[1]          = MGRFRelaxationType::none;
    m_levelInterpType[1]          = MGRInterpolationType::injection;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::cprLikeBlockDiag; // Non-Galerkin Quasi-IMPES CPR
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[1] = 1;
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

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure/temperature reduced system
    setPressureTemperatureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_*/
