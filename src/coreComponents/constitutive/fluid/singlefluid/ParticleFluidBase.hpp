/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParticleFluidBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Base class for particle fluid model kernel wrappers.
 */
class ParticleFluidBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex numElems() const { return m_settlingFactor.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex numGauss() const { return m_settlingFactor.size( 1 ); };

protected:

  /**
   * @brief Constructor.
   * @param density     fluid density
   * @param dDens_dPres derivative of density w.r.t. pressure
   * @param viscosity   fluid viscosity
   * @param dVisc_dPres derivative of viscosity w.r.t. pressure
   */
  ParticleFluidBaseUpdate( bool const isCollisionalSlip,
                           real64 const maxProppantConcentration,
                           arrayView1d< real64 > const & settlingFactor,
                           arrayView1d< real64 > const & dSettlingFactor_dPressure,
                           arrayView1d< real64 > const & dSettlingFactor_dProppantConcentration,
                           arrayView2d< real64 > const & dSettlingFactor_dComponentConcentration,
                           arrayView1d< real64 > const & collisionFactor,
                           arrayView1d< real64 > const & dCollisionFactor_dProppantConcentration,
                           arrayView1d< real64 > const & proppantPackPermeability )
    : m_isCollisionalSlip( isCollisionalSlip ),
    m_maxProppantConcentration( maxProppantConcentration ),
    m_settlingFactor( settlingFactor ),
    m_dSettlingFactor_dPressure( dSettlingFactor_dPressure ),
    m_dSettlingFactor_dProppantConcentration( dSettlingFactor_dProppantConcentration ),
    m_dSettlingFactor_dComponentConcentration( dSettlingFactor_dComponentConcentration ),
    m_collisionFactor( collisionFactor ),
    m_dCollisionFactor_dProppantConcentration( dCollisionFactor_dProppantConcentration ),
    m_proppantPackPermeability( proppantPackPermeability )
  {}

  /**
   * @brief Copy constructor.
   */
  ParticleFluidBaseUpdate( ParticleFluidBaseUpdate const & ) = default;

  /**
   * @brief Move constructor.
   */
  ParticleFluidBaseUpdate( ParticleFluidBaseUpdate && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return reference to this object
   */
  GEOS_HOST_DEVICE
  ParticleFluidBaseUpdate & operator=( ParticleFluidBaseUpdate const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return reference to this object
   */
  GEOS_HOST_DEVICE
  ParticleFluidBaseUpdate & operator=( ParticleFluidBaseUpdate && ) = delete;

  bool m_isCollisionalSlip;

  real64 m_maxProppantConcentration;

  arrayView1d< real64 > m_settlingFactor;
  arrayView1d< real64 > m_dSettlingFactor_dPressure;
  arrayView1d< real64 > m_dSettlingFactor_dProppantConcentration;
  arrayView2d< real64 > m_dSettlingFactor_dComponentConcentration;

  arrayView1d< real64 > m_collisionFactor;
  arrayView1d< real64 > m_dCollisionFactor_dProppantConcentration;

  arrayView1d< real64 > m_proppantPackPermeability;

private:

  /**
   * @brief Perform a single point constitutive update.
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] proppantConcentration proppant concentration value
   * @param[in] fluidDensity fluid density
   * @param[in] dFluidDensity_dPressure derivatives of the fluid density wrt the pressure
   * @param[in] dFluidDensity_dComponentConcentration derivatives of the fluid density wrt the composition
   * @param[in] fluidViscosity fluid viscosity
   * @param[in] dFluidViscosity_dPressure derivatives of the fluid viscosity wrt the pressure
   * @param[in] dFluidViscosity_dComponentConcentration derivatives of the fluid viscosity wrt the composition
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       real64 const proppantConcentration,
                       real64 const fluidDensity,
                       real64 const dFluidDensity_dPressure,
                       arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                       real64 const fluidViscosity,
                       real64 const dFluidViscosity_dPressure,
                       arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration ) const = 0;
};

/**
 * Base class for models calculating proppant transport and settling properties
 */
class ParticleFluidBase : public ConstitutiveBase
{
public:

  ParticleFluidBase( string const & name, Group * const parent );

  virtual ~ParticleFluidBase() override;

  // *** ConstitutiveBase interface
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr localIndex MAX_NUM_COMPONENTS = 4;

protected:

  virtual void postInputInitialization() override;

  array1d< real64 > m_settlingFactor;
  array1d< real64 > m_dSettlingFactor_dPressure;
  array1d< real64 > m_dSettlingFactor_dProppantConcentration;
  array2d< real64 > m_dSettlingFactor_dComponentConcentration;

  array1d< real64 > m_collisionFactor;
  array1d< real64 > m_dCollisionFactor_dProppantConcentration;

  array1d< real64 > m_proppantPackPermeability;

  integer m_isCollisionalSlip;

  real64 m_maxProppantConcentration;

private:

  // *** Data repository keys
  struct viewKeyStruct
  {
    static constexpr char const * maxProppantConcentrationString() { return "maxProppantConcentration"; }
    static constexpr char const * isCollisionalSlipString() { return "isCollisionalSlip"; }
  };

};

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PARTICLEFLUIDBASE_HPP_
