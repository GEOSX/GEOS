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
 * @file SinglePhaseThermalConductivity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITY_HPP_
#define GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITY_HPP_

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"


namespace geos
{
namespace constitutive
{

/**
 * @brief The update class for constant thermal conductivity (does not do anything)
 */
class SinglePhaseThermalConductivityUpdate : public SinglePhaseThermalConductivityBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the thermal conductivity updates
   * @param effectiveConductivity the array of cell-wise effective conductivities in the subregion
   * the subregion
   */
  SinglePhaseThermalConductivityUpdate( arrayView3d< real64 > const & effectiveConductivity,
                                        arrayView3d< real64 > const & dEffectiveConductivity_dT )
    : SinglePhaseThermalConductivityBaseUpdate( effectiveConductivity,
                                                dEffectiveConductivity_dT )
  {}

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & laggedPorosity ) const override
  { GEOS_UNUSED_VAR( k, q, laggedPorosity ); }

};

/**
 * @brief The class for constant thermal conductivity
 */
class SinglePhaseThermalConductivity : public SinglePhaseThermalConductivityBase
{
public:

  /**
   * @brief Constructor for the class storing constant conductivity
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  SinglePhaseThermalConductivity( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "SinglePhaseThermalConductivity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void initializeRockFluidState( arrayView2d< real64 const > const & initialPorosity ) const override final;

  virtual void updateFromTemperature( arrayView1d< real64 const > const & temperature ) const override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = SinglePhaseThermalConductivityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_effectiveConductivity,
                          m_dEffectiveConductivity_dT );
  }

  struct viewKeyStruct : public SinglePhaseThermalConductivityBase::viewKeyStruct
  {
    static constexpr char const * defaultThermalConductivityComponentsString() { return "defaultThermalConductivityComponents"; }
    static constexpr char const * thermalConductivityGradientComponentsString() { return "thermalConductivityGradientComponents"; }
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }
  } viewKeys;

protected:

  virtual void postInputInitialization() override;

private:

  /// Default thermal conductivity components in the subRegion
  R1Tensor m_defaultThermalConductivityComponents;

  /// Thermal conductivity gradient components in the subRegion
  R1Tensor m_thermalConductivityGradientComponents;

  /// Reference temperature
  real64 m_referenceTemperature;

};

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_CONSTANTTHERMALCONDUCTIVITY_HPP_
