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
 * @file ReactiveMultiFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP_


#include "codingUtilities/EnumStrings.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/chemicalReactions/EquilibriumReactions.hpp"
#include "constitutive/fluid/chemicalReactions/KineticReactions.hpp"

#include <memory>

namespace geosx
{

namespace constitutive
{

class ReactiveMultiFluid : public MultiFluidBase
{
public:

  using exec_policy = parallelDevicePolicy<>;

  ReactiveMultiFluid( string const & name,
                      Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "ReactiveMultiFluid";}

  virtual string getCatalogName() const override { return catalogName(); }

  virtual bool isThermal() const override;

  virtual integer getWaterPhaseIndex() const override { return 0; }

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
  class KernelWrapper final : public MultiFluidBase::KernelWrapper
  {

public:

    GEOSX_HOST_DEVICE
    void computeChemistry( real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & kineticReactionRates ) const;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class ReactiveMultiFluid;

    KernelWrapper( arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   bool const isThermal,
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
                   arrayView2d< real64 > const & primarySpeciesConcentration,
                   arrayView2d< real64 > const & secondarySpeciesConcentration,
                   arrayView2d< real64 > const & primarySpeciesTotalConcentration,
                   arrayView2d< real64 > const & kineticReactionRates );

    /// Reaction related terms
    integer m_numPrimarySpecies;

    chemicalReactions::EquilibriumReactions::KernelWrapper m_equilibriumReactions;

    chemicalReactions::KineticReactions::KernelWrapper m_kineticReactions;

    arrayView2d< real64 >  m_primarySpeciesConcentration;

    arrayView2d< real64 >  m_secondarySpeciesConcentration;

    arrayView2d< real64 >  m_primarySpeciesTotalConcentration;

    arrayView2d< real64 >  m_kineticReactionRates;
  };

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {};

protected:

  virtual void postProcessInput() override;

  void createChemicalReactions();

  virtual void resizeFields( localIndex const size, localIndex const numPts ) override;

  /// Reaction related terms
  integer m_numPrimarySpecies;

  integer m_numSecondarySpecies;

  integer m_numKineticReactions;

  std::unique_ptr< chemicalReactions::EquilibriumReactions > m_equilibriumReactions;

  std::unique_ptr< chemicalReactions::KineticReactions > m_kineticReactions;

  array2d< real64 >  m_primarySpeciesConcentration;

  array2d< real64 >  m_secondarySpeciesConcentration;

  array2d< real64 >  m_primarySpeciesTotalConcentration;

  array2d< real64 >  m_kineticReactionRates;
};

GEOSX_HOST_DEVICE
inline void
ReactiveMultiFluid::KernelWrapper::
  computeChemistry( real64 const pressure,
                    real64 const temperature,
                    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & kineticReactionRates ) const
{
  // I am assuming that the primary variable is the concentration of the primary species.
  for( int i=0; i < m_numPrimarySpecies; i++ )
  {
    primarySpeciesTotalConcentration[i] = composition[i];
  }

  m_equilibriumReactions.updateConcentrations( temperature,
                                               primarySpeciesTotalConcentration,
                                               primarySpeciesConcentration,
                                               secondarySpeciesConcentration );

  m_kineticReactions.computeReactionRates( temperature,
                                           primarySpeciesConcentration,
                                           secondarySpeciesConcentration,
                                           kineticReactionRates );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP
