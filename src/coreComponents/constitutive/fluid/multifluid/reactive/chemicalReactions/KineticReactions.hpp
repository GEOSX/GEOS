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
 * @file EquilibriumReaction.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_CHEMICALREACTIONS_REACTIONBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_CHEMICALREACTIONS_REACTIONBASE_HPP_

#include "ReactionsBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "common/PhysicsConstants.hpp"

namespace geos
{

namespace constitutive
{

namespace chemicalReactions
{

class KineticReactions : public ReactionsBase
{
public:

  KineticReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies, integer const numKineticReactions );

  class KernelWrapper final : public ReactionsBase::KernelWrapper
  {
public:

    static constexpr real64 RConst = constants::gasConstant;

    KernelWrapper( integer const numPrimarySpecies,
                   integer const numSecondarySpecies,
                   integer const numKineticReactions,
                   arrayView1d< real64 > const & log10EqConst,
                   arrayView2d< real64 > const & stoichMatrix,
                   arrayView1d< integer > const & chargePrimary,
                   arrayView1d< integer > const & chargeSec,
                   arrayView1d< real64 > const & ionSizePrimary,
                   arrayView1d< real64 > const & ionSizeSec,
                   real64 const DebyeHuckelA,
                   real64 const DebyeHuckelB,
                   real64 const WATEQBDot,
                   arrayView1d< real64 > const & reactionRateConstant,
                   real64 const specificSurfaceArea ):
      ReactionsBase::KernelWrapper( numPrimarySpecies,
                                    numSecondarySpecies,
                                    log10EqConst,
                                    stoichMatrix,
                                    chargePrimary,
                                    chargeSec,
                                    ionSizePrimary,
                                    ionSizeSec,
                                    DebyeHuckelA,
                                    DebyeHuckelB,
                                    WATEQBDot ),
      m_reactionRateConstant( reactionRateConstant ),
      m_numKineticReactions( numKineticReactions ),
      m_specificSurfaceArea( specificSurfaceArea )
    {}

    /**
     * @brief Compute kinetic reaction rates.
     *
     * @param temperature
     * @param primarySpeciesConcentration concentration of the primary species
     * @param log10PrimaryActCoeff
     * @param specificSurfaceArea the surface area available per unit volume
     * @param reactionRates
     */
    void computeReactionRates( real64 const & temperature,
                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                               arraySlice1d< real64, compflow::USD_COMP - 1 > const & reactionRates ) const;

private:

    arrayView1d< real64 > m_reactionRateConstant;

    integer m_numKineticReactions;

    real64 m_specificSurfaceArea;

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  array1d< real64 > m_reactionRateConstant;

  integer m_numKineticReactions;

  real64 m_specificSurfaceArea;
};

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_
