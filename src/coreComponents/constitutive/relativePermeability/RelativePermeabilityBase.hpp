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
 * @file RelativePermeabilityBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
#define GEOS_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief enum to dispatch thress phase interpolation
 */
enum class ThreePhaseInterpolator : integer
{
  BAKER,
  STONEII
};

/**
 * @brief stringify enum to dispatch three phase interpolators
 */
ENUM_STRINGS( ThreePhaseInterpolator,
              "BAKER",
              "STONEII" );

class RelativePermeabilityBaseUpdate
{
public:

  /**
   * @brief Get handle to relperm
   * @return arrayView for the [elt][gauss][phase] relperm container
   */

  GEOS_HOST_DEVICE
  arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > relperm() const
  { return m_phaseRelPerm; }

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_phaseRelPerm.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_phaseRelPerm.size( 1 ); }

  /**
   * @brief Get number of fluid phases.
   * @return number of phases
   */
  GEOS_HOST_DEVICE
  integer numPhases() const { return LvArray::integerConversion< integer >( m_phaseTypes.size() ); }


protected:

  RelativePermeabilityBaseUpdate( arrayView1d< integer const > const & phaseTypes,
                                  arrayView1d< integer const > const & phaseOrder,
                                  arrayView3d< real64, constitutive::relperm::USD_RELPERM > const & phaseRelPerm,
                                  arrayView4d< real64, constitutive::relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
                                  arrayView3d< real64, constitutive::relperm::USD_RELPERM > const & phaseTrappedVolFrac )
    : m_phaseTypes( phaseTypes ),
    m_phaseOrder( phaseOrder ),
    m_phaseRelPerm( phaseRelPerm ),
    m_dPhaseRelPerm_dPhaseVolFrac( dPhaseRelPerm_dPhaseVolFrac ),
    m_phaseTrappedVolFrac( phaseTrappedVolFrac ) {}

  arrayView1d< integer const > m_phaseTypes;
  arrayView1d< integer const > m_phaseOrder;

  arrayView3d< real64, constitutive::relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64, constitutive::relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  arrayView3d< real64, constitutive::relperm::USD_RELPERM > m_phaseTrappedVolFrac;

private:
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const = 0;
};


class RelativePermeabilityBase : public ConstitutiveBase
{
public:

  static constexpr integer MAX_NUM_PHASES = 3;

  struct PhaseType
  {
    enum : integer
    {
      OIL = 0,
      GAS = 1,
      WATER = 2,
    };
  };

  /// order of the phase properties in the water-oil data
  struct WaterOilPairPhaseType
  {
    enum : integer
    {
      WATER = 0, ///< first water phase property
      OIL = 1,   ///< second oil phase property
    };
  };

  /// order of the phase properties in the gas-oil data
  struct GasOilPairPhaseType
  {
    enum : integer
    {
      GAS = 0, ///< first gas phase property
      OIL = 1, ///< second oil phase property
    };
  };

  RelativePermeabilityBase( string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > phaseTrappedVolFraction() const { return m_phaseTrappedVolFrac; }

  arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > phaseRelPerm() const { return m_phaseRelPerm; }
  arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > dPhaseRelPerm_dPhaseVolFraction() const { return m_dPhaseRelPerm_dPhaseVolFrac; }

  arrayView1d< integer const > getPhaseOrder() const { return m_phaseOrder; }
  virtual arrayView1d< real64 const > getPhaseMinVolumeFraction() const = 0;

  std::tuple< integer, integer > wettingAndNonWettingPhaseIndices() const;
  /**
   * @brief Save converged phase volume fraction at the end of a time step (needed for hysteresis)
   * @param[in] phaseVolFraction an array containing the phase volume fractions at the end of a converged time step
   */
  virtual void saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const
  { GEOS_UNUSED_VAR( phaseVolFraction ); }

  /**
   * @brief Save converged state after the relperm update
   */
  virtual void saveConvergedState() const override;

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * phaseTypesString() { return "phaseTypes"; }
    static constexpr char const * phaseOrderString() { return "phaseOrder"; }
  };

private:

  /**
   * @brief Called internally to set array dim labels.
   */
  void setLabels();

protected:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  virtual void resizeFields( localIndex const size, localIndex const numPts );

  virtual void postInputInitialization() override;

  // phase names read from input
  string_array m_phaseNames;

  // phase ordering info
  array1d< integer > m_phaseTypes;
  array1d< integer > m_phaseOrder;

protected:

  // output quantities
  array3d< real64, constitutive::relperm::LAYOUT_RELPERM > m_phaseRelPerm;
  array3d< real64, constitutive::relperm::LAYOUT_RELPERM >  m_phaseRelPerm_n;
  array4d< real64, constitutive::relperm::LAYOUT_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;
  array3d< real64, constitutive::relperm::LAYOUT_RELPERM > m_phaseTrappedVolFrac;

};

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
