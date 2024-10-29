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
 * @file MultiPhaseThermalConductivityBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYBASE_HPP
#define GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief The abstract base class to perform the thermal conductivity
 */
class MultiPhaseThermalConductivityBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_effectiveConductivity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_effectiveConductivity.size( 1 ); }

  /**
   * @brief Get number of fluid phases.
   * @return number of phases
   */
  GEOS_HOST_DEVICE
  integer numPhases() const { return LvArray::integerConversion< integer >( m_dEffectiveConductivity_dPhaseVolFrac.size( 3 ) ); }

protected:

  /**
   * @brief Constructor for the class performing the thermal conductivity updates
   * @param effectiveConductivity the array of cell-wise effective conductivities in the subregion
   * @param dEffectiveConductivity_dPhaseVolFrac the array of cell-wise derivatives of effective conductivities wrt phase vol fractions in
   * the subregion
   */
  MultiPhaseThermalConductivityBaseUpdate( arrayView3d< real64 > const & effectiveConductivity,
                                           arrayView4d< real64 > const & dEffectiveConductivity_dPhaseVolFrac )
    : m_effectiveConductivity( effectiveConductivity ),
    m_dEffectiveConductivity_dPhaseVolFrac( dEffectiveConductivity_dPhaseVolFrac )
  {}

  /// View on the cell-wise effective conductivities
  arrayView3d< real64 > m_effectiveConductivity;

  /// View on the cell-wise derivatives of effective conductivities wrt phase volume fractions
  arrayView4d< real64 > m_dEffectiveConductivity_dPhaseVolFrac;

private:

  /**
   * @brief Pointwise update function called from the solver
   * @param[in] k index of the cell in the subRegion
   * @param[in] q constitutive index (equal to one in this class)
   * @param[in] laggedPorosity lagged porosity in the cell (for fractures, this will be unused)
   * @param[in] phaseVolFrac phase volume fraction in the cell
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & laggedPorosity,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac ) const = 0;
};

/**
 * @brief The abstract base class for thermal conductivity
 */
class MultiPhaseThermalConductivityBase : public ConstitutiveBase
{
public:

  /// Max number of phases allowed in the class
  static constexpr integer MAX_NUM_PHASES = 3;

  /**
   * @brief Constructor for the abstract base class
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  MultiPhaseThermalConductivityBase( string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Initialize the thermal conductivity state (needed when thermal conductivity depends on porosity and phase volume fraction)
   * @param[in] initialPorosity the initial porosity field after reservoir initialization
   * @param[in] initialPhaseVolumeFraction the initial phase volume fraction field
   *
   * Note: this is needed because for now, the porosity and phase volume fraction are treated **explictly**
   */
  virtual void initializeRockFluidState( arrayView2d< real64 const > const & initialPorosity,
                                         arrayView2d< real64 const, compflow::USD_PHASE > const & initialPhaseVolumeFraction ) const
  { GEOS_UNUSED_VAR( initialPorosity, initialPhaseVolumeFraction ); }

  /**
   * @brief Save the thermal conductivity state (needed when thermal conductivity depends on porosity and phase volume fraction)
   * @param[in] convergedPorosity the converged porosity field after reservoir initialization
   * @param[in] convergedPhaseVolumeFraction the converged phase volume fraction field
   *
   * Note: this is needed because for now, the porosity and phase volume fraction are treated **explictly**
   */
  virtual void saveConvergedRockFluidState( arrayView2d< real64 const > const & convergedPorosity,
                                            arrayView2d< real64 const, compflow::USD_PHASE > const & convergedPhaseVolumeFraction ) const
  { GEOS_UNUSED_VAR( convergedPorosity, convergedPhaseVolumeFraction ); }

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
   * @brief Getter for the effective conductivities in the subRegion
   * @return an arrayView of effective conductivities
   */
  arrayView3d< real64 const > effectiveConductivity() const { return m_effectiveConductivity; }

  /**
   * @brief Getter for the derivatives of effective conductivity wrt phase volume fraction in the subRegion
   * @return an arrayView of effective conductivities
   */
  arrayView4d< real64 const > dEffectiveConductivity_dPhaseVolFraction() const { return m_dEffectiveConductivity_dPhaseVolFrac; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
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

  /// cell-wise effective conductivities in the subregion
  array3d< real64 > m_effectiveConductivity;

  /// cell-wise derivatives of effective conductivities wrt phase volume fraction in the subregion
  /// note that for now, the effectiveConductivity is evaluated explicitly, so this field is not used in the solver
  array4d< real64 > m_dEffectiveConductivity_dPhaseVolFrac;

  // the effectiveConductivity is used to compute fluxes, so it is easier to lag porosity for now

};

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_MULTIPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYBASE_HPP
