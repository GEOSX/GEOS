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
 * @file EquilibriumReactions.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_CHEMICALREACTIONS_EQUILIBRIUMREACTIONS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_CHEMICALREACTIONS_EQUILIBRIUMREACTIONS_HPP_

#include "ReactionsBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"

namespace geos
{

namespace constitutive
{

namespace chemicalReactions
{

class EquilibriumReactions : public ReactionsBase
{
public:

  EquilibriumReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies );

  class KernelWrapper final : public ReactionsBase::KernelWrapper
  {
public:

    KernelWrapper( integer const numPrimarySpecies,
                   integer const numSecondarySpecies,
                   arrayView1d< real64 > const & log10EqConst,
                   arrayView2d< real64 > const & stoichMatrix,
                   arrayView1d< integer > const & chargePrimary,
                   arrayView1d< integer > const & chargeSec,
                   arrayView1d< real64 > const & ionSizePrimary,
                   arrayView1d< real64 > const & ionSizeSec,
                   real64 const DebyeHuckelA,
                   real64 const DebyeHuckelB,
                   real64 const WATEQBDot ):
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
                                    WATEQBDot )
    {}

    /**
     * @brief Construct a new update Concentrations object
     *
     * @param temperature
     * @param totalConc
     * @param dLog10PrimaryConc_dTotalConc
     */
    void updateConcentrations( real64 const temperature,
                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                               arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesContentration,
                               arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;
private:

    void assembleEquilibriumReactionSystem( real64 const temperature,
                                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                            arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                            arraySlice2d< real64 > const & matrix,
                                            arraySlice1d< real64 > const & rhs ) const;

    void computeSeondarySpeciesConcAndDerivative( real64 const temperature,
                                                  arraySlice1d< real64 const > const & log10PrimaryActCoeff,
                                                  arraySlice1d< real64 const > const & dLog10PrimaryActCoeff_dIonicStrength,
                                                  arraySlice1d< real64 const > const & log10SecActCoeff,
                                                  arraySlice1d< real64 const > const & dLog10SecActCoeff_dIonicStrength,
                                                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                  arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                  arraySlice2d< real64 > const & dLog10SecConc_dLog10PrimaryConc ) const;

    void computeTotalConcAndDerivative( real64 const temperature,
                                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                        arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                        arraySlice1d< real64 > const & totalConc,
                                        arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc ) const;

    void updatePrimarySpeciesConcentrations( arraySlice1d< real64 const > const solution,
                                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const;

    void setInitialGuess( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                          arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const;

    static constexpr integer m_maxNumIterations = MultiFluidConstants::maxNewtonIterations;
    static constexpr real64 m_newtonTol = 1e-6;
  };

/**
 * @brief Create an update kernel wrapper.
 * @return the wrapper
 */
  KernelWrapper createKernelWrapper() const;

};


} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REQUILIBRIUMREACTIONS_HPP_
