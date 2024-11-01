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
 * @file DiffusionBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DIFFUSION_DIFFUSIONBASE_HPP
#define GEOS_CONSTITUTIVE_DIFFUSION_DIFFUSIONBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief The abstract base class to perform the diffusion update
 */
class DiffusionBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_diffusivity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_diffusivity.size( 1 ); }

protected:

  /**
   * @brief Constructor for the class performing the diffusion updates
   * @param diffusivity the array of cell-wise diffusion in the subregion
   * @param dDiffusivity_dTemperature the array of cell-wise derivatives of diffusion wrt temperature in the subregion
   */
  DiffusionBaseUpdate( arrayView3d< real64 > const & diffusivity,
                       arrayView3d< real64 > const & dDiffusivity_dTemperature )
    : m_diffusivity( diffusivity ),
    m_dDiffusivity_dTemperature( dDiffusivity_dTemperature )
  {}

  /// View on the cell-wise diffusivity
  arrayView3d< real64 > const m_diffusivity;

  /// View on the cell-wise derivatives of diffusivity wrt temperature
  arrayView3d< real64 > const m_dDiffusivity_dTemperature;

private:

  /**
   * @brief Pointwise update function called from the solver
   * @param[in] k index of the cell in the subRegion
   * @param[in] q constitutive index (equal to one in this class)
   * @param[in] temperature temperature in the cell
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & temperature ) const = 0;
};

/**
 * @brief The abstract base class for diffusion
 */
class DiffusionBase : public ConstitutiveBase
{
public:

  /// Max number of phases allowed in the class
  static constexpr integer MAX_NUM_PHASES = 3;

  /**
   * @brief Constructor for the abstract base class
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  DiffusionBase( string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Getter for the number of fluid phases
   * @return the number of fluid phases
   */
  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  /**
   * @brief Getter for the phase names
   * @return an arrayView of phase names
   */
  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  /**
   * @brief Getter for the phase diffusivity multipliers
   * @return an arrayView of multipliers
   */
  arrayView3d< real64 const > phaseDiffusivityMultiplier() const { return m_phaseDiffusivityMultiplier; }

  /**
   * @brief Getter for the diffusivities in the subRegion
   * @return an arrayView of diffusivities
   */
  arrayView3d< real64 const > diffusivity() const { return m_diffusivity; }

  /**
   * @brief Getter for the derivatives of diffusivity wrt temperature in the subRegion
   * @return an arrayView of derivatives of diffusivities wrt temperature
   */
  arrayView3d< real64 const > dDiffusivity_dTemperature() const { return m_dDiffusivity_dTemperature; }

  /**
   * @brief Initialize the diffusivity state (needed when diffusion depends on temperature)
   * @param[in] initialTemperature the initial temperature field after reservoir initialization
   *
   * Note: this is needed because for now, the temperature is treated **explicitly** in the diffusion tensor
   */
  virtual void initializeTemperatureState( arrayView1d< real64 const > const & initialTemperature ) const
  { GEOS_UNUSED_VAR( initialTemperature ); }

  /**
   * @brief Save the temperature state (needed when diffusion depends on temperature)
   * @param[in] convergedTemperature the converged temperature field
   *
   * Note: this is needed because for now, the temperature is treated **explicitly** in the diffusion tensor
   */
  virtual void saveConvergedTemperatureState( arrayView1d< real64 const > const & convergedTemperature ) const
  { GEOS_UNUSED_VAR( convergedTemperature ); }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * defaultPhaseDiffusivityMultiplierString() { return "defaultPhaseDiffusivityMultipliers"; }
  };

private:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void resizeFields( localIndex const size, localIndex const numPts );

protected:

  virtual void postInputInitialization() override;

  /// phase names read from input
  string_array m_phaseNames;

  /// array of phase-wise default diffusion multipliers
  array1d< real64 > m_defaultPhaseDiffusivityMultiplier;

  /// array of phase-wise diffusion multipliers
  array3d< real64 > m_phaseDiffusivityMultiplier;

  /// cell-wise diffusivities in the subregion
  array3d< real64 > m_diffusivity;

  /// cell-wise derivatives of diffusivities wrt temperature in the subregion
  array3d< real64 > m_dDiffusivity_dTemperature;

};

} // namespace constitutive

} // namespace geos


#endif // GEOS_CONSTITUTIVE_DIFFUSION_DIFFUSIONBASE_HPP
