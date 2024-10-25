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
 * @file MultiPhaseConstantThermalConductivity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_CONSTANTTHERMALCONDUCTIVITY_HPP_
#define GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_CONSTANTTHERMALCONDUCTIVITY_HPP_

#include "constitutive/thermalConductivity/MultiPhaseThermalConductivityBase.hpp"


namespace geos
{
namespace constitutive
{

/**
 * @brief The update class for constant thermal conductivity (does not do anything)
 */
class MultiPhaseConstantThermalConductivityUpdate : public MultiPhaseThermalConductivityBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the thermal conductivity updates
   * @param effectiveConductivity the array of cell-wise effective conductivities in the subregion
   * @param dEffectiveConductivity_dPhaseVolFrac the array of cell-wise derivatives of effective conductivities wrt phase vol fractions in
   * the subregion
   */
  MultiPhaseConstantThermalConductivityUpdate( arrayView3d< real64 > const & effectiveConductivity,
                                               arrayView4d< real64 > const & dEffectiveConductivity_dPhaseVolFrac )
    : MultiPhaseThermalConductivityBaseUpdate( effectiveConductivity, dEffectiveConductivity_dPhaseVolFrac )
  {}

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & laggedPorosity,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override
  { GEOS_UNUSED_VAR( k, q, laggedPorosity, phaseVolFraction ); }

};

/**
 * @brief The class for constant thermal conductivity
 */
class MultiPhaseConstantThermalConductivity : public MultiPhaseThermalConductivityBase
{
public:

  /**
   * @brief Constructor for the class storing constant conductivity
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  MultiPhaseConstantThermalConductivity( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "MultiPhaseConstantThermalConductivity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = MultiPhaseConstantThermalConductivityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_effectiveConductivity,
                          m_dEffectiveConductivity_dPhaseVolFrac );
  }

  struct viewKeyStruct : public MultiPhaseThermalConductivityBase::viewKeyStruct
  {
    static constexpr char const * thermalConductivityComponentsString() { return "thermalConductivityComponents"; }
  } viewKeys;

protected:

  virtual void postInputInitialization() override;

private:

  /// default thermal conductivity in the subRegion
  R1Tensor m_thermalConductivityComponents;

};

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_CONSTANTTHERMALCONDUCTIVITY_HPP_
