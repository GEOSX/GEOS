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
 * @file SlurryFluidBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SLURRYFLUIDBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SLURRYFLUIDBASE_HPP_

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Base class for slurry fluid model kernel wrappers.
 */
class SlurryFluidBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex numElems() const { return m_density.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex numGauss() const { return m_density.size( 1 ); };

  /**
   * @brief Get number of fluid components.
   * @return number of fluid components
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex numFluidComponents() const { return m_defaultComponentDensity.size(); }

protected:

  /**
   * @brief
   * @param defaultDensity
   * @param defaultCompressibility
   * @param defaultViscosity
   * @param nIndices
   * @param Ks
   * @param isNewtonianFluid
   * @param density
   * @param dDens_dPres
   * @param dDens_dProppantConc
   * @param dDens_dCompConc
   * @param componentDensity
   * @param dCompDens_dPres
   * @param dCompDens_dCompConc
   * @param fluidDensity
   * @param dFluidDens_dPres
   * @param dFluidDens_dCompConc
   * @param fluidViscosity
   * @param dFluidVisc_dPres
   * @param dFluidVisc_dCompConc
   * @param viscosity
   * @param dVisc_dPres
   * @param dVisc_dProppantConc
   * @param dVisc_dCompConc
   */
  SlurryFluidBaseUpdate( arrayView1d< real64 const > const & defaultDensity,
                         arrayView1d< real64 const > const & defaultCompressibility,
                         arrayView1d< real64 const > const & defaultViscosity,
                         arrayView1d< real64 const > const & nIndices,
                         arrayView1d< real64 const > const & Ks,
                         bool const isNewtonianFluid,
                         arrayView2d< real64 > const & density,
                         arrayView2d< real64 > const & dDens_dPres,
                         arrayView2d< real64 > const & dDens_dProppantConc,
                         arrayView3d< real64 > const & dDens_dCompConc,
                         arrayView3d< real64 > const & componentDensity,
                         arrayView3d< real64 > const & dCompDens_dPres,
                         arrayView4d< real64 > const & dCompDens_dCompConc,
                         arrayView2d< real64 > const & fluidDensity,
                         arrayView2d< real64 > const & dFluidDens_dPres,
                         arrayView3d< real64 > const & dFluidDens_dCompConc,
                         arrayView2d< real64 > const & fluidViscosity,
                         arrayView2d< real64 > const & dFluidVisc_dPres,
                         arrayView3d< real64 > const & dFluidVisc_dCompConc,
                         arrayView2d< real64 > const & viscosity,
                         arrayView2d< real64 > const & dVisc_dPres,
                         arrayView2d< real64 > const & dVisc_dProppantConc,
                         arrayView3d< real64 > const & dVisc_dCompConc )
    : m_defaultComponentDensity( defaultDensity ),
    m_defaultComponentCompressibility( defaultCompressibility ),
    m_defaultComponentViscosity( defaultViscosity ),
    m_nIndices( nIndices ),
    m_Ks( Ks ),
    m_isNewtonianFluid( isNewtonianFluid ),
    m_density( density ),
    m_dDensity_dPressure( dDens_dPres ),
    m_dDensity_dProppantConc( dDens_dProppantConc ),
    m_dDensity_dCompConc( dDens_dCompConc ),
    m_componentDensity( componentDensity ),
    m_dCompDens_dPres( dCompDens_dPres ),
    m_dCompDens_dCompConc( dCompDens_dCompConc ),
    m_fluidDensity( fluidDensity ),
    m_dFluidDens_dPres( dFluidDens_dPres ),
    m_dFluidDens_dCompConc( dFluidDens_dCompConc ),
    m_fluidViscosity( fluidViscosity ),
    m_dFluidVisc_dPres( dFluidVisc_dPres ),
    m_dFluidVisc_dCompConc( dFluidVisc_dCompConc ),
    m_viscosity( viscosity ),
    m_dViscosity_dPressure( dVisc_dPres ),
    m_dViscosity_dProppantConc( dVisc_dProppantConc ),
    m_dViscosity_dCompConc( dVisc_dCompConc )
  {}

  /**
   * @brief Copy constructor.
   */
  SlurryFluidBaseUpdate( SlurryFluidBaseUpdate const & ) = default;

  /**
   * @brief Move constructor.
   */
  SlurryFluidBaseUpdate( SlurryFluidBaseUpdate && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return reference to this object
   */
  GEOS_HOST_DEVICE
  SlurryFluidBaseUpdate & operator=( SlurryFluidBaseUpdate const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return reference to this object
   */
  GEOS_HOST_DEVICE
  SlurryFluidBaseUpdate & operator=( SlurryFluidBaseUpdate && ) = delete;

  arrayView1d< real64 const > m_defaultComponentDensity;
  arrayView1d< real64 const > m_defaultComponentCompressibility;
  arrayView1d< real64 const > m_defaultComponentViscosity;

  arrayView1d< real64 const > m_nIndices;
  arrayView1d< real64 const > m_Ks;

  bool m_isNewtonianFluid;

  arrayView2d< real64 > m_density;
  arrayView2d< real64 > m_dDensity_dPressure;
  arrayView2d< real64 > m_dDensity_dProppantConc;
  arrayView3d< real64 > m_dDensity_dCompConc;

  arrayView3d< real64 > m_componentDensity;
  arrayView3d< real64 > m_dCompDens_dPres;
  arrayView4d< real64 > m_dCompDens_dCompConc;

  arrayView2d< real64 > m_fluidDensity;
  arrayView2d< real64 > m_dFluidDens_dPres;
  arrayView3d< real64 > m_dFluidDens_dCompConc;

  arrayView2d< real64 > m_fluidViscosity;
  arrayView2d< real64 > m_dFluidVisc_dPres;
  arrayView3d< real64 > m_dFluidVisc_dCompConc;

  arrayView2d< real64 > m_viscosity;
  arrayView2d< real64 > m_dViscosity_dPressure;
  arrayView2d< real64 > m_dViscosity_dProppantConc;
  arrayView3d< real64 > m_dViscosity_dCompConc;

private:

  /**
   * @brief Perform a single point constitutive update.
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] q second constitutive index (e.g. quadrature index)
   * @param[in] pressure pressure value
   * @param[in] proppantConcentration proppant concentration value
   * @param[in] Componentconcentration fluid composition array
   * @param[in] shearRate shear rate for power-law fluid calculation
   * @param[in] isProppantBoundary proppant boundary flag
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const proppantConcentration,
                       arraySlice1d< real64 const > const & componentConcentration,
                       real64 const shearRate,
                       integer const isProppantBoundary ) const = 0;

  /**
   * @brief Perform a single point fluid property update.
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] q second constitutive index (e.g. quadrature index)
   * @param[in] pressure pressure value
   * @param[in] componentconcentration fluid composition array
   * @param[in] shearRate shear rate for power-law fluid calculation
   */
  GEOS_HOST_DEVICE
  virtual void updateFluidProperty( localIndex const k,
                                    localIndex const q,
                                    real64 const pressure,
                                    arraySlice1d< real64 const > const & componentConcentration,
                                    real64 const shearRate ) const = 0;

  /**
   * @brief Perform a single point fluid component property update.
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] q second constitutive index (e.g. quadrature index)
   * @param[in] pressure pressure value
   * @param[in] componentconcentration fluid composition array
   */
  GEOS_HOST_DEVICE
  virtual void updateComponentDensity( localIndex const k,
                                       localIndex const q,
                                       real64 const pressure,
                                       arraySlice1d< real64 const > const & componentConcentration ) const = 0;

};

/**
 * Base class for models calculating slurry fluid properties.
 */
class SlurryFluidBase : public SingleFluidBase
{
public:

  SlurryFluidBase( string const & name, Group * const parent );

  virtual ~SlurryFluidBase() override;

  // *** ConstitutiveBase interface
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr localIndex MAX_NUM_COMPONENTS = 3;

  localIndex numFluidComponents() const;

  arrayView1d< real64 > kIndex() const { return m_Ks; }
  arrayView1d< real64 const > nIndex() const { return m_nIndices; }

  arrayView2d< real64 > dDensity_dProppantConcentration() { return m_dDensity_dProppantConc; }
  arrayView2d< real64 const > dDensity_dProppantConcentration() const { return m_dDensity_dProppantConc; }

  arrayView3d< real64 > dDensity_dComponentConcentration() { return m_dDensity_dCompConc; }
  arrayView3d< real64 const > dDensity_dComponentConcentration() const { return m_dDensity_dCompConc; }

  arrayView2d< real64 > fluidDensity() { return m_fluidDensity; }
  arrayView2d< real64 const > fluidDensity() const { return m_fluidDensity; }

  arrayView2d< real64 > dFluidDensity_dPressure() { return m_dFluidDens_dPres; }
  arrayView2d< real64 const > dFluidDensity_dPressure() const { return m_dFluidDens_dPres; }

  arrayView3d< real64 > dFluidDensity_dComponentConcentration() { return m_dFluidDens_dCompConc; }
  arrayView3d< real64 const > dFluidDensity_dComponentConcentration() const { return m_dFluidDens_dCompConc; }

  arrayView3d< real64 > componentDensity() { return m_componentDensity; }
  arrayView3d< real64 const > componentDensity() const { return m_componentDensity; }

  arrayView3d< real64 > dComponentDensity_dPressure() { return m_dCompDens_dPres; }
  arrayView3d< real64 const > dComponentDensity_dPressure() const { return m_dCompDens_dPres; }

  arrayView4d< real64 > dComponentDensity_dComponentConcentration() { return m_dCompDens_dCompConc; }
  arrayView4d< real64 const > dComponentDensity_dComponentConcentration() const { return m_dCompDens_dCompConc; }

  arrayView2d< real64 > fluidViscosity(){ return m_fluidViscosity; }
  arrayView2d< real64 const > fluidViscosity() const { return m_fluidViscosity; }

  arrayView2d< real64 > dFluidViscosity_dPressure() { return m_dFluidVisc_dPres; }
  arrayView2d< real64 const > dFluidViscosity_dPressure() const { return m_dFluidVisc_dPres; }

  arrayView3d< real64 > dFluidViscosity_dComponentConcentration() { return m_dFluidVisc_dCompConc; }
  arrayView3d< real64 const > dFluidViscosity_dComponentConcentration() const { return m_dFluidVisc_dCompConc; }

  arrayView2d< real64 > dViscosity_dProppantConcentration() { return m_dViscosity_dProppantConc; }
  arrayView2d< real64 const > dViscosity_dProppantConcentration() const { return m_dViscosity_dProppantConc; }

  arrayView3d< real64 > dViscosity_dComponentConcentration() { return m_dViscosity_dCompConc; }
  arrayView3d< real64 const > dViscosity_dComponentConcentration() const { return m_dViscosity_dCompConc; }

  bool isNewtonianFluid() const { return m_isNewtonianFluid; }

protected:

  virtual void postInputInitialization() override;

  string_array m_componentNames;

  array1d< real64 > m_defaultComponentDensity;
  array1d< real64 > m_defaultComponentCompressibility;
  array1d< real64 > m_defaultComponentViscosity;

  array1d< real64 > m_nIndices;
  array1d< real64 > m_Ks;

  array2d< real64 > m_dDensity_dProppantConc;
  array3d< real64 > m_dDensity_dCompConc;

  array3d< real64 > m_componentDensity;
  array3d< real64 > m_dCompDens_dPres;
  array4d< real64 > m_dCompDens_dCompConc;

  array2d< real64 > m_fluidDensity;
  array2d< real64 > m_dFluidDens_dPres;
  array3d< real64 > m_dFluidDens_dCompConc;

  array2d< real64 > m_fluidViscosity;
  array2d< real64 > m_dFluidVisc_dPres;
  array3d< real64 > m_dFluidVisc_dCompConc;

  array2d< real64 > m_dViscosity_dProppantConc;
  array3d< real64 > m_dViscosity_dCompConc;

  bool m_isNewtonianFluid;

private:

  // *** Data repository keys
  struct viewKeyStruct
  {
    static constexpr char const * componentNamesString() { return "componentNames"; }

    static constexpr char const * defaultComponentDensityString() { return "defaultComponentDensity"; }
    static constexpr char const * defaultCompressibilityString() { return "defaultCompressibility"; }
    static constexpr char const * defaultComponentViscosityString() { return "defaultComponentViscosity"; }

    static constexpr char const * flowBehaviorIndexString() { return "flowBehaviorIndex"; }
    static constexpr char const * flowConsistencyIndexString() { return "flowConsistencyIndex"; }
  };
};

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SLURRYFLUIDBASE_HPP_
