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
 * @file MultiPhaseVolumeWeightedThermalConductivity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_VOLUMEWEIGHTEDTHERMALCONDUCTIVITY_HPP_
#define GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_VOLUMEWEIGHTEDTHERMALCONDUCTIVITY_HPP_

#include "constitutive/thermalConductivity/MultiPhaseThermalConductivityBase.hpp"


namespace geos
{
namespace constitutive
{

/**
 * @brief The update class for volume-weighted thermal conductivity
 */
class MultiPhaseVolumeWeightedThermalConductivityUpdate : public MultiPhaseThermalConductivityBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the thermal conductivity updates
   * @param effectiveConductivity the array of cell-wise effective conductivities in the subregion
   * @param dEffectiveConductivity_dPhaseVolFrac the array of cell-wise derivatives of effective conductivities wrt phase vol fractions in
   * the subregion
   * @param rockThermalConductivity the array of cell-wise rock thermal conductivities
   * @param phaseThermalConductivity the array of fluid phase thermal conductivities
   */
  MultiPhaseVolumeWeightedThermalConductivityUpdate( arrayView3d< real64 > const & effectiveConductivity,
                                                     arrayView4d< real64 > const & dEffectiveConductivity_dPhaseVolFrac,
                                                     arrayView3d< real64 const > const & rockThermalConductivity,
                                                     arrayView1d< real64 const > const & phaseThermalConductivity )
    : MultiPhaseThermalConductivityBaseUpdate( effectiveConductivity, dEffectiveConductivity_dPhaseVolFrac ),
    m_rockThermalConductivity( rockThermalConductivity ),
    m_phaseThermalConductivity( phaseThermalConductivity )
  {}

  /**
   * @brief Cell-wise compute function for the volume-weighted effective conductivity
   * @param[in] laggedPorosity the cell lagged porosity
   * @param[in] phaseVolFrac the cell phase volume fraction
   * @param[in] rockThermalConductivity the cell rock thermal conductivity
   * @param[out] effectiveConductivity the cell effective conductivity
   * @param[out] dEffectiveConductivity_dPhaseVolFrac the cell derivative of effective conductivity wrt phase volume fractions
   */
  GEOS_HOST_DEVICE
  void compute( real64 const & laggedPorosity,
                arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
                arraySlice1d< real64 const > const rockThermalConductivity,
                arraySlice1d< real64 > const effectiveConductivity,
                arraySlice2d< real64 > const dEffectiveConductivity_dPhaseVolFrac ) const
  {
    LvArray::forValuesInSlice( effectiveConductivity, []( real64 & val ){ val = 0.0; } );
    for( integer ip = 0; ip < numPhases(); ++ip )
    {
      real64 const phaseVolFracTimesConductivity = phaseVolFrac[ip] * m_phaseThermalConductivity[ip];
      real64 const porosityTimesConductivity = m_phaseThermalConductivity[ip] * laggedPorosity;
      for( integer dir = 0; dir < 3; ++dir )
      {
        effectiveConductivity[dir] += phaseVolFracTimesConductivity;
        dEffectiveConductivity_dPhaseVolFrac[dir][ip] = porosityTimesConductivity;
      }
    }
    for( integer dir = 0; dir < 3; ++dir )
    {
      effectiveConductivity[dir] *= laggedPorosity;
      effectiveConductivity[dir] += ( 1.0 - laggedPorosity ) * rockThermalConductivity[dir];
    }
  }

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & laggedPorosity,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac ) const override
  {
    compute( laggedPorosity,
             phaseVolFrac,
             m_rockThermalConductivity[k][q],
             m_effectiveConductivity[k][q],
             m_dEffectiveConductivity_dPhaseVolFrac[k][q] );
  }

private:

  /// View on the cell-wise rock thermal conductivity
  arrayView3d< real64 const > m_rockThermalConductivity;

  /// View on the phase thermal conductivity
  arrayView1d< real64 const > m_phaseThermalConductivity;

};

/**
 * @brief The class for volume-weighted thermal conductivity
 */
class MultiPhaseVolumeWeightedThermalConductivity : public MultiPhaseThermalConductivityBase
{
public:

  /**
   * @brief Constructor for the class implementing volume-weighted conductivity
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  MultiPhaseVolumeWeightedThermalConductivity( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void initializeRockFluidState( arrayView2d< real64 const > const & initialPorosity,
                                         arrayView2d< real64 const, compflow::USD_PHASE > const & initialPhaseVolumeFraction ) const override;

  virtual void saveConvergedRockFluidState( arrayView2d< real64 const > const & convergedPorosity,
                                            arrayView2d< real64 const, compflow::USD_PHASE > const & convergedPhaseVolumeFraction ) const override;

  static string catalogName() { return "MultiPhaseVolumeWeightedThermalConductivity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = MultiPhaseVolumeWeightedThermalConductivityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_effectiveConductivity,
                          m_dEffectiveConductivity_dPhaseVolFrac,
                          m_rockThermalConductivity,
                          m_phaseThermalConductivity );
  }

  struct viewKeyStruct : public MultiPhaseThermalConductivityBase::viewKeyStruct
  {
    static constexpr char const * rockThermalConductivityComponentsString() { return "rockThermalConductivityComponents"; }
    static constexpr char const * phaseThermalConductivityString() { return "phaseThermalConductivity"; }
  } viewKeys;

protected:

  virtual void postInputInitialization() override;

private:

  /// default rock thermal conductivity in the subRegion
  R1Tensor m_rockThermalConductivityComponents;

  /// cell-wise rock thermal conductivity in the subRegion
  array3d< real64 > m_rockThermalConductivity;

  /// fluid phase thermal conductivity
  array1d< real64 > m_phaseThermalConductivity;

};

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_VOLUMEWEIGHTEDTHERMALCONDUCTIVITY_HPP_
