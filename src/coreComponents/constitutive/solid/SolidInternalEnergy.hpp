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
 * @file SolidInternalEnergy.hpp
 */


#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDINTERNALENERGY_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDINTERNALENERGY_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{


class SolidInternalEnergyUpdates
{
public:

  SolidInternalEnergyUpdates( arrayView2d< real64 > const & internalEnergy,
                              arrayView2d< real64 > const & dInternalEnergy_dTemperature,
                              real64 const & referenceVolumetricHeatCapacity,
                              real64 const & dVolumetricHeatCapacity_dTemperature,
                              real64 const & referenceTemperature,
                              real64 const & referenceInternalEnergy ):
    m_internalEnergy( internalEnergy ),
    m_dInternalEnergy_dTemperature( dInternalEnergy_dTemperature ),
    m_referenceVolumetricHeatCapacity( referenceVolumetricHeatCapacity ),
    m_dVolumetricHeatCapacity_dTemperature( dVolumetricHeatCapacity_dTemperature ),
    m_referenceTemperature( referenceTemperature ),
    m_referenceInternalEnergy( referenceInternalEnergy )
  {}

  GEOS_DEVICE
  void update( localIndex const k,
               real64 const & temperature ) const
  {
    compute( temperature,
             m_internalEnergy[k][0],
             m_dInternalEnergy_dTemperature[k][0] );
  }

  GEOS_DEVICE
  void compute( real64 const & temperature,
                real64 & internalEnergy,
                real64 & dInternalEnergy_dTemperature ) const
  {
    real64 volumetricHeatCapacity = m_referenceVolumetricHeatCapacity + m_dVolumetricHeatCapacity_dTemperature * ( temperature - m_referenceTemperature );

    internalEnergy = m_referenceInternalEnergy + volumetricHeatCapacity * ( temperature - m_referenceTemperature );
    dInternalEnergy_dTemperature =  volumetricHeatCapacity + m_dVolumetricHeatCapacity_dTemperature * ( temperature - m_referenceTemperature );
  }

private:

  /// Solid internal energy
  arrayView2d< real64 > m_internalEnergy;

  /// Derivative of the solid internal energy w.r.t. the temperature
  arrayView2d< real64 > m_dInternalEnergy_dTemperature;

  /// Solid volumetric heat capacity at the reference tempearture
  real64 m_referenceVolumetricHeatCapacity;

  /// Derivative of the solid volumetric heat capacity w.r.t. the temperature
  real64 m_dVolumetricHeatCapacity_dTemperature;

  /// Reference temperature
  real64 m_referenceTemperature;

  /// Internal energy at the reference temperature
  real64 m_referenceInternalEnergy;
};

class SolidInternalEnergy : public ConstitutiveBase
{
public:

  SolidInternalEnergy( string const & name, Group * const parent );

  static string catalogName() { return "SolidInternalEnergy"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * internalEnergyString() { return "internalEnergy"; }
    static constexpr char const * oldInternalEnergyString() { return "internalEnergy_n"; }
    static constexpr char const * dInternalEnergy_dTemperatureString() { return "dInternalEnergy_dTemperature"; }
    static constexpr char const * referenceVolumetricHeatCapacityString() { return "referenceVolumetricHeatCapacity"; }
    static constexpr char const * dVolumetricHeatCapacity_dTemperatureString() { return "dVolumetricHeatCapacity_dTemperature"; }
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }
    static constexpr char const * referenceInternalEnergyString() { return "referenceInternalEnergy"; }
  } viewKeys;

  using KernelWrapper = SolidInternalEnergyUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper( m_internalEnergy,
                          m_dInternalEnergy_dTemperature,
                          m_referenceVolumetricHeatCapacity,
                          m_dVolumetricHeatCapacity_dTemperature,
                          m_referenceTemperature,
                          m_referenceInternalEnergy );
  }


  /**
   * @brief Const accessor for internalEnergy.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getInternalEnergy() const { return m_internalEnergy; }

  /**
   * @brief Const/non-mutable accessor for internalEnergy_n.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getInternalEnergy_n() const { return m_internalEnergy_n; }

  /**
   * @brief Const/non-mutable accessor for dInternalEnergy_dTemperature.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getDinternalEnergy_dTemperature() const { return m_dInternalEnergy_dTemperature; }

  /// Save state data in preparation for next timestep
  virtual void saveConvergedState() const override;

private:

  /// Solid internal energy
  array2d< real64 > m_internalEnergy;

  /// Old solid internal energy
  array2d< real64 > m_internalEnergy_n;

  /// Derivative of the solid internal energy w.r.t. the temperature
  array2d< real64 > m_dInternalEnergy_dTemperature;

  /// Solid volumetric heat capacity at the reference temperature
  real64 m_referenceVolumetricHeatCapacity;

  /// Derivative of the solid volumetric heat capacity w.r.t. the temperature
  real64 m_dVolumetricHeatCapacity_dTemperature;

  /// Reference temperature
  real64 m_referenceTemperature;

  /// Reference internal energy
  real64 m_referenceInternalEnergy;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDINTERNALENERGY_HPP_ */
