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
 * @file ConstantDiffusion.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DIFFUSION_CONSTANTDIFFUSION_HPP_
#define GEOS_CONSTITUTIVE_DIFFUSION_CONSTANTDIFFUSION_HPP_

#include "constitutive/diffusion/DiffusionBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief The update class for constant diffusion (does not do anything)
 */
class ConstantDiffusionUpdate : public DiffusionBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the diffusion updates
   * @param diffusivity the array of cell-wise diffusivities in the subregion
   * @param dDiffusivity_dTemperature the array of cell-wise derivatives of diffusivities wrt temperature in the subregion
   */
  ConstantDiffusionUpdate( arrayView3d< real64 > const & diffusivity,
                           arrayView3d< real64 > const & dDiffusivity_dTemperature )
    : DiffusionBaseUpdate( diffusivity, dDiffusivity_dTemperature )
  {}

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & temperature ) const override
  { GEOS_UNUSED_VAR( k, q, temperature ); }

};

/**
 * @brief The class for constant diffusion
 */
class ConstantDiffusion : public DiffusionBase
{
public:

  /**
   * @brief Constructor for the class storing constant diffusion
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  ConstantDiffusion( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ConstantDiffusion"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ConstantDiffusionUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_diffusivity,
                          m_dDiffusivity_dTemperature );
  }

  struct viewKeyStruct : public DiffusionBase::viewKeyStruct
  {
    static constexpr char const * diffusivityComponentsString() { return "diffusivityComponents"; }
  } viewKeys;

protected:

  virtual void postInputInitialization() override;

private:

  /// default diffusivity in the subRegion
  array1d< real64 > m_diffusivityComponents;

};

} // namespace constitutive

} // namespace geos


#endif // GEOS_CONSTITUTIVE_DIFFUSION_CONSTANTDIFFUSION_HPP_
