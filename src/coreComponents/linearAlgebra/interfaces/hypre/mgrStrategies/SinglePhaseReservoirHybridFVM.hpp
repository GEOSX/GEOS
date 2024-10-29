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
 * @file SinglePhaseReservoirHybridFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRHYBRIDFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRHYBRIDFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhaseReservoirHybridFVM strategy.
 *
 * Labels description stored in point_marker_array
 *   dofLabel: 0 = cell-centered reservoir pressure
 *   dofLabel: 1 = face-centered reservoir pressure
 *   dofLabel: 2 = well pressure
 *   dofLabel: 3 = well rate
 *
 * Ingredients
 *
 * 2-level MGR reduction strategy
 * 1st level: eliminate the well block
 * 2nd level: eliminate the cell-centered reservoir pressure
 * The coarse grid is the face-centered reservoir pressure (Lagrange multiplier) system
 * The coarse grid solved with BoomerAMG
 */
class SinglePhaseReservoirHybridFVM : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit SinglePhaseReservoirHybridFVM( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 4 ) )
  {
    // Level 0: eliminate the cell-centered pressure
    m_labels[0].push_back( 0 );
    m_labels[0].push_back( 1 );
    // Level 1: eliminate the face-centered pressure
    m_labels[1].push_back( 1 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]         = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[1]        = 1;
    m_levelInterpType[1]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
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
    if( !mgrParams.areWellsShut )
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

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRHYBRIDFVM_HPP_*/
