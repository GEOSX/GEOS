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
 * @file TableCapillaryPressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSURE_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{

class TableCapillaryPressure : public CapillaryPressureBase
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

  TableCapillaryPressure( std::string const & name, dataRepository::Group * const parent );

  static std::string catalogName() { return "TableCapillaryPressure"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  class KernelWrapper final : public CapillaryPressureBaseUpdate
  {
public:

    KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & capPresKernelWrappers,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPres,
                   arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPres_dPhaseVolFrac );

    GEOS_HOST_DEVICE
    void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                  arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                  arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const;

    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;

private:

    /// Array of kernel wrappers for the capillary pressures
    /// Is of size 1 for two-phase flow, and of size 2 for three-phase flow
    arrayView1d< TableFunction::KernelWrapper const > const m_capPresKernelWrappers;

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * wettingNonWettingCapPresTableNameString() { return "wettingNonWettingCapPressureTableName"; }
    static constexpr char const * wettingIntermediateCapPresTableNameString() { return "wettingIntermediateCapPressureTableName"; }
    static constexpr char const * nonWettingIntermediateCapPresTableNameString() { return "nonWettingIntermediateCapPressureTableName"; }
    static constexpr char const * capPresWrappersString() { return "capPresWrappers"; }

  };

private:

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /// Capillary pressure table names (one for each phase in the wetting-non-wetting pair)
  string m_wettingNonWettingCapPresTableName;

  /// Capillary pressure table names (one for each phase in the wetting-intermediate pair)
  string m_wettingIntermediateCapPresTableName;

  /// Capillary pressure table names (one for each phase in the non-wetting-intermediate pair)
  string m_nonWettingIntermediateCapPresTableName;

  /// Capillary pressure kernel wrapper for the first pair (wetting-intermediate if NP=3, wetting-non-wetting otherwise)
  array1d< TableFunction::KernelWrapper > m_capPresKernelWrappers;

};

GEOS_HOST_DEVICE
inline void
TableCapillaryPressure::KernelWrapper::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
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
    using TPT = TableCapillaryPressure::ThreePhasePairPhaseType;

    // water-oil capillary pressure
    phaseCapPres[ipWater] =
      m_capPresKernelWrappers[TPT::INTERMEDIATE_WETTING].compute( &(phaseVolFraction)[ipWater],
                                                                  &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );

    // gas-oil capillary pressure
    phaseCapPres[ipGas] =
      m_capPresKernelWrappers[TPT::INTERMEDIATE_NONWETTING].compute( &(phaseVolFraction)[ipGas],
                                                                     &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );

    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( ipWater < 0 )
  {
    // put capillary pressure on the non-wetting phase
    phaseCapPres[ipGas] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipGas],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );

    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( ipOil < 0 || ipGas < 0 )
  {
    // put capillary pressure on the wetting phase
    phaseCapPres[ipWater] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipWater],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );
  }
}

GEOS_HOST_DEVICE
inline void
TableCapillaryPressure::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          arraySlice1d< geos::real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const
{
  compute( phaseVolFraction,
           m_phaseCapPressure[k][q],
           m_dPhaseCapPressure_dPhaseVolFrac[k][q] );
}

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSURE_HPP
