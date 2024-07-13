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
 * @file CompositionalMultiphaseReservoirFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief CompositionalMultiphaseReservoirFVM strategy.
 *
 * Labels description stored in point_marker_array
 *                0 = reservoir pressure
 *                1 = reservoir density
 *              ... = ... (reservoir densities)
 * numResLabels - 1 = reservoir density
 *     numResLabels = well pressure
 * numResLabels + 1 = well density
 *              ... = ... (well densities)
 * numResLabels + numWellLabels - 2 = well density
 * numResLabels + numWellLabels - 1 = well rate
 *
 * 3-level MGR reduction strategy
 *   - 1st level: eliminate the well block
 *   - 2nd level: eliminate the reservoir density associated with the volume constraint
 *   - 3rd level: eliminate the remaining the reservoir densities
 *   - The coarse grid is the pressure system and is solved with BoomerAMG
 */
class CompositionalMultiphaseReservoirFVM : public MGRStrategyBase< 3 >
{
public:

  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseReservoirFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    HYPRE_Int const numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );

    // Level 0: eliminate the well block
    m_labels[0].resize( numResLabels );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    // Level 1: eliminate the last density of the reservoir block
    m_labels[1].resize( numResLabels - 1 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 0 );
    // Level 2: eliminate the rest of the densities
    m_labels[2].push_back( 0 );

    setupLabels();

    // level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    // level 1
    m_levelFRelaxType[1]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;

    // level 2
    m_levelFRelaxType[2]          = MGRFRelaxationType::none;
    m_levelInterpType[2]          = MGRInterpolationType::injection;
    m_levelRestrictType[2]        = MGRRestrictionType::blockColLumped; // True-IMPES
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[2] = 1;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams parameters for the configuration of the MGR recipe
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    // if the wells are shut, using Gaussian elimination as F-relaxation for the well block is an overkill
    // in that case, we just use Jacobi
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[0] = MGRFRelaxationType::jacobi;
    }

    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_*/
