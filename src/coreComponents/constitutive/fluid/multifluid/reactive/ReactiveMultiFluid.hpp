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
 * @file ReactiveMultiFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_REACTIVEMULTIFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_REACTIVEMULTIFLUID_HPP_


#include "common/format/EnumStrings.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/reactive/chemicalReactions/EquilibriumReactions.hpp"
#include "constitutive/fluid/multifluid/reactive/chemicalReactions/KineticReactions.hpp"

#include <memory>

namespace geos
{

namespace constitutive
{

class ReactiveMultiFluid : public MultiFluidBase
{
public:

  using exec_policy = serialPolicy;

  ReactiveMultiFluid( string const & name,
                      Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  virtual bool isThermal() const override;

  arrayView2d< real64 const, compflow::USD_COMP > primarySpeciesConcentration() const
  { return m_primarySpeciesConcentration; }

  arrayView2d< real64 const, compflow::USD_COMP > secondarySpeciesConcentration() const
  { return m_secondarySpeciesConcentration; }

  arrayView2d< real64 const, compflow::USD_COMP > kineticReactionRates() const
  { return m_kineticReactionRates; }

  integer numPrimarySpecies() const { return m_numPrimarySpecies; }

  integer numSecondarySpecies() const { return m_numSecondarySpecies; }

  integer numKineticReactions() const { return m_numKineticReactions; }

  /**
   * @brief Kernel wrapper class for ReactiveMultiFluid.
   */
  class KernelWrapper : public MultiFluidBase::KernelWrapper
  {

public:

    void computeChemistry( real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & kineticReactionRates ) const;

    virtual void updateChemistry( localIndex const k,
                                  localIndex const q,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const = 0;

    /**
     * @brief Construct a new Kernel Wrapper object
     *
     * @param componentMolarWeight
     * @param useMass
     * @param isThermal
     * @param phaseFraction
     * @param phaseDensity
     * @param phaseMassDensity
     * @param phaseViscosity
     * @param phaseEnthalpy
     * @param phaseInternalEnergy
     * @param phaseCompFraction
     * @param totalDensity
     * @param numPrimarySpecies
     * @param equilibriumReactions
     * @param kineticReactions
     * @param primarySpeciesConcentration
     * @param secondarySpeciesConcentration
     * @param primarySpeciesTotalConcentration
     * @param kineticReactionRates
     */
    KernelWrapper( arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity,
                   integer const numPrimarySpecies,
                   chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                   chemicalReactions::KineticReactions const & kineticReactions,
                   arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesConcentration,
                   arrayView2d< real64, compflow::USD_COMP > const & secondarySpeciesConcentration,
                   arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesTotalConcentration,
                   arrayView2d< real64, compflow::USD_COMP > const & kineticReactionRates ):
      MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                     useMass,
                                     std::move( phaseFraction ),
                                     std::move( phaseDensity ),
                                     std::move( phaseMassDensity ),
                                     std::move( phaseViscosity ),
                                     std::move( phaseEnthalpy ),
                                     std::move( phaseInternalEnergy ),
                                     std::move( phaseCompFraction ),
                                     std::move( totalDensity )  ),
      m_numPrimarySpecies( numPrimarySpecies ),
      m_equilibriumReactions( equilibriumReactions.createKernelWrapper() ),
      m_kineticReactions( kineticReactions.createKernelWrapper() ),
      m_primarySpeciesConcentration( primarySpeciesConcentration ),
      m_secondarySpeciesConcentration( secondarySpeciesConcentration ),
      m_primarySpeciesTotalConcentration( primarySpeciesTotalConcentration ),
      m_kineticReactionRates( kineticReactionRates )
    {}

protected:

    void convertMoleFractionToMolarity( real64 const totalDensity,
                                        real64 const totalMolecularWeight,
                                        arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition,
                                        arraySlice1d< geos::real64, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration ) const;

    friend class ReactiveMultiFluid;
    /// Reaction related terms
    integer m_numPrimarySpecies;

    chemicalReactions::EquilibriumReactions::KernelWrapper m_equilibriumReactions;

    chemicalReactions::KineticReactions::KernelWrapper m_kineticReactions;

    arrayView2d< real64, compflow::USD_COMP >  m_primarySpeciesConcentration;

    arrayView2d< real64, compflow::USD_COMP >  m_secondarySpeciesConcentration;

    arrayView2d< real64, compflow::USD_COMP >  m_primarySpeciesTotalConcentration;

    arrayView2d< real64, compflow::USD_COMP >  m_kineticReactionRates;
  };

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {};

protected:

  virtual void postInputInitialization() override;

  void createChemicalReactions();

  virtual void resizeFields( localIndex const size, localIndex const numPts ) override;

  /// Reaction related terms
  integer m_numPrimarySpecies;

  integer m_numSecondarySpecies;

  integer m_numKineticReactions;

  std::unique_ptr< chemicalReactions::EquilibriumReactions > m_equilibriumReactions;

  std::unique_ptr< chemicalReactions::KineticReactions > m_kineticReactions;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_primarySpeciesConcentration;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_secondarySpeciesConcentration;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_primarySpeciesTotalConcentration;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_kineticReactionRates;
};

inline void
ReactiveMultiFluid::KernelWrapper::
  computeChemistry( real64 const pressure,
                    real64 const temperature,
                    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & kineticReactionRates ) const
{
  GEOS_UNUSED_VAR( pressure );

  // 2. solve for equilibrium
  m_equilibriumReactions.updateConcentrations( temperature,
                                               primarySpeciesTotalConcentration,
                                               primarySpeciesConcentration,
                                               secondarySpeciesConcentration );

  // 3. compute kinetic reaction rates
  m_kineticReactions.computeReactionRates( temperature,
                                           primarySpeciesConcentration,
                                           secondarySpeciesConcentration,
                                           kineticReactionRates );
}


inline void
ReactiveMultiFluid::KernelWrapper::
  convertMoleFractionToMolarity( real64 const totalDensity,
                                 real64 const totalMolecularWeight,
                                 arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition,
                                 arraySlice1d< geos::real64, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration ) const
{
  // 1. Convert from mole fraction to molarity ( mol/L )
  real64 const conversionFactor = totalDensity / totalMolecularWeight * 1e-3;  //conversion to L instead of cubic meters
  for( int i=0; i < m_numPrimarySpecies; i++ )
  {
    primarySpeciesTotalConcentration[i] = composition[i] * conversionFactor;
  }
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP
