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
 * @file ImmiscibleMultiphaseFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRIMMISCIBLEMULTIPHASEFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRIMMISCIBLEMULTIPHASEFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief ImmiscibleMultiphaseFVM strategy.
 *
 * Labels description stored in point_marker_array
 *               0 = pressure
 *               1 = phase volume fraction
 *             ... = phase volume fractions
 *   numLabels - 1 = phase volume fraction
 *
 * 1-level MGR reduction strategy:
 *   - 1st level: eliminate the reservoir phase volume fractions
 *   - The coarse grid (pressure system) is solved with BoomerAMG.
 *
 */
class ImmiscibleMultiphaseFVM : public MGRStrategyBase< 1 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit ImmiscibleMultiphaseFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] ) )
  {
    // Level 0: eliminate last density which corresponds to the volume constraint equation
    // m_labels[0].resize( m_numBlocks - 1 );
    // std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );

    // Level 0: eliminate the phase volume fractions
    m_labels[0].push_back( 0 );

    setupLabels();

    m_levelFRelaxType[0]          = MGRFRelaxationType::none;
    m_levelInterpType[0]          = MGRInterpolationType::injection;
    m_levelRestrictType[0]        = MGRRestrictionType::blockColLumped; // True-IMPES
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[0] = 1;
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

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_*/
