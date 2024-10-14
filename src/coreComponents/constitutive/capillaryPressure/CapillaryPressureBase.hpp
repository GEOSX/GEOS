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
 * @file CapillaryPressureBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREBASE_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREBASE_HPP

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace constitutive
{

class CapillaryPressureBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_phaseCapPressure.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_phaseCapPressure.size( 1 ); }

  /**
   * @brief Get number of fluid phases.
   * @return number of phases
   */
  GEOS_HOST_DEVICE
  integer numPhases() const { return LvArray::integerConversion< integer >( m_phaseTypes.size() ); }

protected:

  CapillaryPressureBaseUpdate( arrayView1d< integer const > const & phaseTypes,
                               arrayView1d< integer const > const & phaseOrder,
                               arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPressure,
                               arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPressure_dPhaseVolFrac )
    : m_phaseTypes( phaseTypes ),
    m_phaseOrder( phaseOrder ),
    m_phaseCapPressure( phaseCapPressure ),
    m_dPhaseCapPressure_dPhaseVolFrac( dPhaseCapPressure_dPhaseVolFrac )
  {}

  arrayView1d< integer const > m_phaseTypes;
  arrayView1d< integer const > m_phaseOrder;

  arrayView3d< real64, cappres::USD_CAPPRES > m_phaseCapPressure;
  arrayView4d< real64, cappres::USD_CAPPRES_DS > m_dPhaseCapPressure_dPhaseVolFrac;

private:

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const = 0;
};

class CapillaryPressureBase : public ConstitutiveBase
{
public:

  static constexpr integer MAX_NUM_PHASES = 3;

  struct PhaseType
  {
    enum : integer
    {
      OIL            = 0,
      GAS            = 1,
      WATER          = 2,
    };
  };

  // choose the reference pressure to be the oil pressure for all models
  static constexpr integer REFERENCE_PHASE = PhaseType::OIL;

  CapillaryPressureBase( string const & name,
                         dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Initialize the capillary pressure state (needed when capillary pressure depends on porosity and permeability)
   * @param[in] initialPorosity the initial porosity field after reservoir initialization
   * @param[in] initialPermeability the initial permeability field
   */
  virtual void initializeRockState( arrayView2d< real64 const > const & initialPorosity,
                                    arrayView3d< real64 const > const & initialPermeability ) const
  { GEOS_UNUSED_VAR( initialPorosity, initialPermeability ); }

  /**
   * @brief Save the capillary pressure state (needed when capillary pressure depends on porosity and permeability)
   * @param[in] convergedPorosity the converged porosity field after reservoir initialization
   * @param[in] convergedPermeability the converged permeability field
   *
   * Note: for now, this function is only used to save the porosity and permeability.
   *       But, once hysteresis is implemented, it will also be used to save the max historical saturation
   *       At this point, it will makes sense to just rename this function to "saveConvergedState"
   *       (this will require changes in ConstitutiveBase, SolidBase, PorosityBase, etc)
   */
  virtual void saveConvergedRockState( arrayView2d< real64 const > const & convergedPorosity,
                                       arrayView3d< real64 const > const & convergedPermeability ) const
  { GEOS_UNUSED_VAR( convergedPorosity, convergedPermeability ); }

  /*
   * @brief Getter for the number of fluid phases
   * @return the number of fluid phases
   */
  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  /*
   * @brief Getter for the phase names
   * @return an array of phase names
   */
  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  /*
   * @brief Getter for the cell-wise phase capillary pressures
   * @return an array of cell-wise phase capillary pressures
   */
  arrayView3d< real64 const, cappres::USD_CAPPRES > phaseCapPressure() const { return m_phaseCapPressure; }

  /*
   * @brief Getter for the cell-wise derivatives of phase capillary pressures wrt phase volume fractions
   * @return an array of cell-wise derivatives of phase capillary pressures wrt phase volume fractions
   */
  arrayView4d< real64 const, cappres::USD_CAPPRES_DS > dPhaseCapPressure_dPhaseVolFraction() const { return m_dPhaseCapPressure_dPhaseVolFrac; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * phaseTypesString() { return "phaseTypes"; }
    static constexpr char const * phaseOrderString() { return "phaseOrder"; }
  };

private:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void resizeFields( localIndex const size, localIndex const numPts );

  /**
   * @brief Called internally to set array dim labels.
   */
  void setLabels();

protected:

  virtual void postInputInitialization() override;

  // phase names read from input
  string_array m_phaseNames;

  // phase ordering info
  array1d< integer > m_phaseTypes;
  array1d< integer > m_phaseOrder;

  // output quantities
  array3d< real64, cappres::LAYOUT_CAPPRES >  m_phaseCapPressure;
  array4d< real64, cappres::LAYOUT_CAPPRES_DS >  m_dPhaseCapPressure_dPhaseVolFrac;

};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREBASE_HPP
