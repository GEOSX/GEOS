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
 * @file JFunctionCapillaryPressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_JFUNCTIONCAPILLARYPRESSURE_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_JFUNCTIONCAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

#include "codingUtilities/EnumStrings.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{

class JFunctionCapillaryPressure : public CapillaryPressureBase
{
public:

  /// order of the phase properties for three-phase flow
  struct ThreePhasePairPhaseType
  {
    enum : integer
    {
      INTERMEDIATE_WETTING = 0,   ///< index for intermediate-wetting
      INTERMEDIATE_NONWETTING = 1 ///< index for intermediate-non-wetting
    };
  };

  JFunctionCapillaryPressure( std::string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void initializeRockState( arrayView2d< real64 const > const & initialPorosity,
                                    arrayView3d< real64 const > const & initialPermeability ) const override;

  virtual void saveConvergedRockState( arrayView2d< real64 const > const & initialPorosity,
                                       arrayView3d< real64 const > const & initialPermeability ) const override;

  static std::string catalogName() { return "JFunctionCapillaryPressure"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  class KernelWrapper final : public CapillaryPressureBaseUpdate
  {
public:

    KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & jFuncKernelWrappers,
                   arrayView2d< real64 const > const & jFuncMultiplier,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPres,
                   arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPres_dPhaseVolFrac );

    GEOS_HOST_DEVICE
    void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                  arraySlice1d< real64 const > const & jFuncMultiplier,
                  arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                  arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const;

    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;

private:

    /// Array of kernel wrappers for the J-function
    /// Is of size 1 for two-phase flow, and of size 2 for three-phase flow
    arrayView1d< TableFunction::KernelWrapper const > const m_jFuncKernelWrappers;

    /// Array of cell-wise J-function multipliers
    /// The second dimension is of size 1 for two-phase flow, and of size 2 for three-phase flow
    arrayView2d< real64 const > const m_jFuncMultiplier;

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * wettingNonWettingJFuncTableNameString() { return "wettingNonWettingJFunctionTableName"; }
    static constexpr char const * wettingIntermediateJFuncTableNameString() { return "wettingIntermediateJFunctionTableName"; }
    static constexpr char const * nonWettingIntermediateJFuncTableNameString() { return "nonWettingIntermediateJFunctionTableName"; }
    static constexpr char const * wettingNonWettingSurfaceTensionString() { return "wettingNonWettingSurfaceTension"; }
    static constexpr char const * wettingIntermediateSurfaceTensionString() { return "wettingIntermediateSurfaceTension"; }
    static constexpr char const * nonWettingIntermediateSurfaceTensionString() { return "nonWettingIntermediateSurfaceTension"; }
    static constexpr char const * porosityExponentString() { return "porosityExponent"; }
    static constexpr char const * permeabilityExponentString() { return "permeabilityExponent"; }
    static constexpr char const * permeabilityDirectionString() { return "permeabilityDirection"; }
    static constexpr char const * jFunctionWrappersString() { return "jFunctionWrappers"; }
  };

  /**
   * @brief Type of permeability directions
   */
  enum class PermeabilityDirection : integer
  {
    XY, ///< use the average of permx and permy
    X,  ///< use permx only
    Y,  ///< use permy only
    Z,  ///< use permz only
  };


private:

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /// J-function table names (one for each phase in the wetting-non-wetting pair)
  string m_wettingNonWettingJFuncTableName;

  /// J-function table names (one for each phase in the wetting-intermediate pair)
  string m_wettingIntermediateJFuncTableName;

  /// J-function table names (one for each phase in the non-wetting-intermediate pair)
  string m_nonWettingIntermediateJFuncTableName;

  /// Surface tension for the pair wetting-non-wetting phase
  real64 m_wettingNonWettingSurfaceTension;

  /// Surface tension for the pair wetting-intermediate phase
  real64 m_wettingIntermediateSurfaceTension;

  /// Surface tension for the pair non-wetting-intermediate phase
  real64 m_nonWettingIntermediateSurfaceTension;

  /// Porosity exponent used in the multiplier
  real64 m_porosityExponent;

  /// Permeability exponent used in the multiplier
  real64 m_permeabilityExponent;

  /// Permeability direction (XY, X, Y, Z)
  PermeabilityDirection m_permeabilityDirection;

  /// Cell-wise array of J-function multipliers
  array2d< real64 > m_jFuncMultiplier;

  /// J-function kernel wrapper for the first pair (wetting-intermediate if NP=3, wetting-non-wetting otherwise)
  array1d< TableFunction::KernelWrapper > m_jFuncKernelWrappers;

};

GEOS_HOST_DEVICE
inline void
JFunctionCapillaryPressure::KernelWrapper::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64 const > const & jFuncMultiplier,
           arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
           arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseCapPres_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  using PT = CapillaryPressureBase::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];

  if( ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    using TPT = JFunctionCapillaryPressure::ThreePhasePairPhaseType;

    // Step 1: water-oil

    // water-oil J-function
    phaseCapPres[ipWater] =
      m_jFuncKernelWrappers[TPT::INTERMEDIATE_WETTING].compute( &(phaseVolFraction)[ipWater],
                                                                &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );
    // apply water-oil multiplier
    phaseCapPres[ipWater] *= jFuncMultiplier[TPT::INTERMEDIATE_WETTING];
    dPhaseCapPres_dPhaseVolFrac[ipWater][ipWater] *= jFuncMultiplier[TPT::INTERMEDIATE_WETTING];

    // Step 2: gas-oil

    // gas-oil J-function
    phaseCapPres[ipGas] =
      m_jFuncKernelWrappers[TPT::INTERMEDIATE_NONWETTING].compute( &(phaseVolFraction)[ipGas],
                                                                   &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );
    // apply gas-oil multiplier
    phaseCapPres[ipGas] *= jFuncMultiplier[TPT::INTERMEDIATE_NONWETTING];
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= jFuncMultiplier[TPT::INTERMEDIATE_NONWETTING];

    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( ipWater < 0 )
  {
    // put J-function on the non-wetting phase
    phaseCapPres[ipGas] =
      m_jFuncKernelWrappers[0].compute( &(phaseVolFraction)[ipGas],
                                        &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );
    // apply multiplier
    phaseCapPres[ipGas] *= jFuncMultiplier[0];
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= jFuncMultiplier[0];

    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( ipOil < 0 || ipGas < 0 )
  {
    // put J-function on the wetting phase
    phaseCapPres[ipWater] =
      m_jFuncKernelWrappers[0].compute( &(phaseVolFraction)[ipWater],
                                        &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );
    // apply multiplier
    phaseCapPres[ipWater] *= jFuncMultiplier[0];
    dPhaseCapPres_dPhaseVolFrac[ipWater][ipWater] *= jFuncMultiplier[0];
  }
}

GEOS_HOST_DEVICE
inline void
JFunctionCapillaryPressure::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          arraySlice1d< geos::real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const
{
  compute( phaseVolFraction,
           m_jFuncMultiplier[k],
           m_phaseCapPressure[k][q],
           m_dPhaseCapPressure_dPhaseVolFrac[k][q] );
}

/// Declare strings associated with enumeration values.
ENUM_STRINGS( JFunctionCapillaryPressure::PermeabilityDirection,
              "XY",
              "X",
              "Y",
              "Z" );


} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_JFUNCTIONCAPILLARYPRESSURE_HPP
