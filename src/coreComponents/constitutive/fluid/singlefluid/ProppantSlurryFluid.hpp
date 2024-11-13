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
 * @file ProppantSlurryFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PROPPANTSLURRYFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PROPPANTSLURRYFLUID_HPP_

#include "constitutive/fluid/singlefluid/SlurryFluidBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Update class for the proppant slurry model suitable for lambda capture.
 */
class ProppantSlurryFluidUpdate final : public SlurryFluidBaseUpdate
{
public:

  /**
   * @brief
   * @param compressibility
   * @param referenceProppantDensity
   * @param referencePressure
   * @param referenceDensity
   * @param referenceViscosity
   * @param maxProppantConcentration
   * @param defaultComponentDensity
   * @param defaultCompressibility
   * @param defaultComponentViscosity
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
  ProppantSlurryFluidUpdate( real64 const compressibility,
                             real64 const referenceProppantDensity,
                             real64 const referencePressure,
                             real64 const referenceDensity,
                             real64 const referenceViscosity,
                             real64 const maxProppantConcentration,
                             arrayView1d< real64 const > const & defaultComponentDensity,
                             arrayView1d< real64 const > const & defaultCompressibility,
                             arrayView1d< real64 const > const & defaultComponentViscosity,
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
    : SlurryFluidBaseUpdate( defaultComponentDensity,
                             defaultCompressibility,
                             defaultComponentViscosity,
                             nIndices,
                             Ks,
                             isNewtonianFluid,
                             density,
                             dDens_dPres,
                             dDens_dProppantConc,
                             dDens_dCompConc,
                             componentDensity,
                             dCompDens_dPres,
                             dCompDens_dCompConc,
                             fluidDensity,
                             dFluidDens_dPres,
                             dFluidDens_dCompConc,
                             fluidViscosity,
                             dFluidVisc_dPres,
                             dFluidVisc_dCompConc,
                             viscosity,
                             dVisc_dPres,
                             dVisc_dProppantConc,
                             dVisc_dCompConc ),
    m_compressibility( compressibility ),
    m_referenceProppantDensity( referenceProppantDensity ),
    m_referencePressure( referencePressure ),
    m_referenceDensity( referenceDensity ),
    m_referenceViscosity( referenceViscosity ),
    m_maxProppantConcentration( maxProppantConcentration )
  {}

  /// Default copy constructor
  ProppantSlurryFluidUpdate( ProppantSlurryFluidUpdate const & ) = default;

  /// Default move constructor
  ProppantSlurryFluidUpdate( ProppantSlurryFluidUpdate && ) = default;

  /// Deleted copy assignment operator
  ProppantSlurryFluidUpdate & operator=( ProppantSlurryFluidUpdate const & ) = delete;

  /// Deleted move assignment operator
  ProppantSlurryFluidUpdate & operator=( ProppantSlurryFluidUpdate && ) = delete;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const proppantConcentration,
                       arraySlice1d< real64 const > const & componentConcentration,
                       real64 const shearRate,
                       integer const isProppantBoundary ) const override
  {
    updateFluidProperty( k,
                         q,
                         pressure,
                         componentConcentration,
                         shearRate );

    compute( proppantConcentration,
             m_fluidDensity[k][q],
             m_dFluidDens_dPres[k][q],
             m_dFluidDens_dCompConc[k][q],
             m_fluidViscosity[k][q],
             m_dFluidVisc_dPres[k][q],
             m_dFluidVisc_dCompConc[k][q],
             isProppantBoundary,
             m_density[k][q],
             m_dDensity_dPressure[k][q],
             m_dDensity_dProppantConc[k][q],
             m_dDensity_dCompConc[k][q],
             m_viscosity[k][q],
             m_dViscosity_dPressure[k][q],
             m_dViscosity_dProppantConc[k][q],
             m_dViscosity_dCompConc[k][q] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void updateFluidProperty( localIndex const k,
                                    localIndex const q,
                                    real64 const pressure,
                                    arraySlice1d< real64 const > const & componentConcentration,
                                    real64 const shearRate ) const override
  {
    GEOS_UNUSED_VAR( shearRate );
    computeFluidDensity( pressure,
                         componentConcentration,
                         m_componentDensity[k][q],
                         m_dCompDens_dPres[k][q],
                         m_dCompDens_dCompConc[k][q],
                         m_fluidDensity[k][q],
                         m_dFluidDens_dPres[k][q],
                         m_dFluidDens_dCompConc[k][q] );

    computeFluidViscosity( m_componentDensity[k][q],
                           m_dCompDens_dPres[k][q],
                           m_dCompDens_dCompConc[k][q],
                           m_fluidDensity[k][q],
                           m_dFluidDens_dPres[k][q],
                           m_dFluidDens_dCompConc[k][q],
                           m_fluidViscosity[k][q],
                           m_dFluidVisc_dPres[k][q],
                           m_dFluidVisc_dCompConc[k][q] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void updateComponentDensity( localIndex const k,
                                       localIndex const q,
                                       real64 const pressure,
                                       arraySlice1d< real64 const > const & componentConcentration ) const override
  {
    computeComponentDensity( pressure,
                             componentConcentration,
                             m_componentDensity[k][q],
                             m_dCompDens_dPres[k][q],
                             m_dCompDens_dCompConc[k][q] );
  }

private:

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeFluidDensity( real64 const & pressure,
                            arraySlice1d< real64 const > const & componentConcentration,
                            arraySlice1d< real64 > const & componentDensity,
                            arraySlice1d< real64 > const & dComponentDensity_dPressure,
                            arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration,
                            real64 & fluidDensity,
                            real64 & dFluidDensity_dPressure,
                            arraySlice1d< real64 > const & dFluidDensity_dComponentConcentration ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeComponentDensity( real64 const pressure,
                                arraySlice1d< real64 const > const & componentConcentration,
                                arraySlice1d< real64 > const & componentDensity,
                                arraySlice1d< real64 > const & dComponentDensity_dPressure,
                                arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeFluidViscosity( arraySlice1d< real64 const > const & componentDensity,
                              arraySlice1d< real64 const > const & dComponentDensity_dPressure,
                              arraySlice2d< real64 const > const & dComponentDensity_dComponentConcentration,
                              real64 const fluidDensity,
                              real64 const dFluidDensity_dPressure,
                              arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                              real64 & fluidViscosity,
                              real64 & dFluidViscosity_dPressure,
                              arraySlice1d< real64 > const & dFluidViscosity_dComponentConcentration ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void compute( real64 const & proppantConcentration,
                real64 const & fluidDensity,
                real64 const & dFluidDensity_dPressure,
                arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                real64 const & fluidViscosity,
                real64 const & dFluidViscosity_dPressure,
                arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration,
                integer const & isProppantBoundary,
                real64 & density,
                real64 & dDensity_dPressure,
                real64 & dDensity_dProppantConcentration,
                arraySlice1d< real64 > const & dDensity_dComponentConcentration,
                real64 & viscosity,
                real64 & dViscosity_dPressure,
                real64 & dViscosity_dProppantConcentration,
                arraySlice1d< real64 > const & dViscosity_dComponentConcentration ) const;

  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantConcentration;

};

class ProppantSlurryFluid : public SlurryFluidBase
{
public:

  ProppantSlurryFluid( string const & name, Group * const parent );

  virtual ~ProppantSlurryFluid() override;

  // *** ConstitutiveBase interface
  static string catalogName() { return "ProppantSlurryFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  using KernelWrapper = ProppantSlurryFluidUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  // *** Data repository keys

  struct viewKeyStruct
  {
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referenceDensityString() { return "referenceDensity"; }
    static constexpr char const * referenceProppantDensityString() { return "referenceProppantDensity"; }
    static constexpr char const * maxProppantConcentrationString() { return "maxProppantConcentration"; }
    static constexpr char const * referenceViscosityString() { return "referenceViscosity"; }
  };

  /**
   * @brief get the default value for the density. For the proppant, we employ the reference value
   * as default.
   *
   * @return the default density value;
   */
  real64 defaultDensity() const override final {return m_referenceDensity; };
  /**
   * @brief get the default value for the density. For the proppant, we employ the reference value
   * as default.
   *
   * @return the default viscosity value;
   */
  real64 defaultViscosity() const override final {return m_referenceViscosity; };

protected:

  virtual void postInputInitialization() override;

private:

  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantConcentration;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
ProppantSlurryFluidUpdate::
  computeFluidDensity( real64 const & pressure,
                       arraySlice1d< real64 const > const & componentConcentration,
                       arraySlice1d< real64 > const & componentDensity,
                       arraySlice1d< real64 > const & dComponentDensity_dPressure,
                       arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration,
                       real64 & fluidDensity,
                       real64 & dFluidDensity_dPressure,
                       arraySlice1d< real64 > const & dFluidDensity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  real64 const baseFluidDensity = m_referenceDensity * exp( m_compressibility * (pressure - m_referencePressure));
  real64 const dBaseFluidDensity_dPressure = m_compressibility * baseFluidDensity;

  fluidDensity = baseFluidDensity;
  dFluidDensity_dPressure = dBaseFluidDensity_dPressure;
  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidDensity_dComponentConcentration[c] = 0.0;
  }

  for( localIndex c = 0; c < NC; ++c )
  {
    real64 const density = m_defaultComponentDensity[c] * exp( m_defaultComponentCompressibility[c] * (pressure - m_referencePressure));

    componentDensity[c] = componentConcentration[c] * density;
    dComponentDensity_dPressure[c] = m_defaultComponentCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
    {
      dComponentDensity_dComponentConcentration[c][i] = 0.0;
    }
    dComponentDensity_dComponentConcentration[c][c] = density;

    fluidDensity += componentDensity[c] - componentConcentration[c] * baseFluidDensity;
    dFluidDensity_dPressure += dComponentDensity_dPressure[c] - componentConcentration[c] * dBaseFluidDensity_dPressure;
    dFluidDensity_dComponentConcentration[c] += dComponentDensity_dComponentConcentration[c][c] - baseFluidDensity;
  }

  // TODO: why? - Sergey
  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidDensity_dComponentConcentration[c] = 0.0;
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
ProppantSlurryFluidUpdate::
  computeComponentDensity( real64 const pressure,
                           arraySlice1d< real64 const > const & componentConcentration,
                           arraySlice1d< real64 > const & componentDensity,
                           arraySlice1d< real64 > const & dComponentDensity_dPressure,
                           arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  for( localIndex c = 0; c < NC; ++c )
  {
    real64 const density = m_defaultComponentDensity[c] * exp( m_defaultComponentCompressibility[c] * (pressure - m_referencePressure));

    componentDensity[c] = componentConcentration[c] * density;
    dComponentDensity_dPressure[c] = m_defaultComponentCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
    {
      dComponentDensity_dComponentConcentration[c][i] = 0.0;
    }
    dComponentDensity_dComponentConcentration[c][c] = density;
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
ProppantSlurryFluidUpdate::
  computeFluidViscosity( arraySlice1d< real64 const > const & componentDensity,
                         arraySlice1d< real64 const > const & dComponentDensity_dPressure,
                         arraySlice2d< real64 const > const & GEOS_UNUSED_PARAM( dComponentDensity_dComponentConcentration ),
                         real64 const fluidDensity,
                         real64 const dFluidDensity_dPressure,
                         arraySlice1d< real64 const > const & GEOS_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
                         real64 & fluidViscosity,
                         real64 & dFluidViscosity_dPressure,
                         arraySlice1d< real64 > const & dFluidViscosity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  fluidViscosity = m_referenceViscosity;
  dFluidViscosity_dPressure = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidViscosity_dComponentConcentration[c] = 0.0;
  }

  for( localIndex c1 = 0; c1 < NC; ++c1 )
  {
    fluidViscosity += componentDensity[c1] / fluidDensity * (m_defaultComponentViscosity[c1] - m_referenceViscosity);

    dFluidViscosity_dPressure +=
      (dComponentDensity_dPressure[c1] / fluidDensity - componentDensity[c1] / fluidDensity / fluidDensity * dFluidDensity_dPressure) *
      (m_defaultComponentViscosity[c1] - m_referenceViscosity);
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
ProppantSlurryFluidUpdate::
  compute( real64 const & proppantConcentration,
           real64 const & fluidDensity,
           real64 const & dFluidDensity_dPressure,
           arraySlice1d< real64 const > const & GEOS_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
           real64 const & fluidViscosity,
           real64 const & dFluidViscosity_dPressure,
           arraySlice1d< real64 const > const & GEOS_UNUSED_PARAM( dFluidViscosity_dComponentConcentration ),
           integer const & isProppantBoundary,
           real64 & density,
           real64 & dDensity_dPressure,
           real64 & dDensity_dProppantConcentration,
           arraySlice1d< real64 > const & dDensity_dComponentConcentration,
           real64 & viscosity,
           real64 & dViscosity_dPressure,
           real64 & dViscosity_dProppantConcentration,
           arraySlice1d< real64 > const & dViscosity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  real64 effectiveConcentration = proppantConcentration;
  if( effectiveConcentration > 0.95 * m_maxProppantConcentration || isProppantBoundary == 1 )
  {
    effectiveConcentration = 0.0;
  }

  density = (1.0 - effectiveConcentration) * fluidDensity + effectiveConcentration * m_referenceProppantDensity;
  dDensity_dPressure = (1.0 - effectiveConcentration) * dFluidDensity_dPressure;

  dDensity_dProppantConcentration = 0.0;
  for( localIndex c = 0; c < NC; ++c )
  {
    dDensity_dComponentConcentration[c] = 0.0;
  }

  real64 const coef = pow( 1.0 + 1.25 *  effectiveConcentration / (1.0 - effectiveConcentration / m_maxProppantConcentration), 2.0 );
  viscosity = fluidViscosity * coef;
  dViscosity_dPressure = dFluidViscosity_dPressure * coef;

  dViscosity_dProppantConcentration = 0.0;
  for( localIndex c = 0; c < NC; ++c )
  {
    dViscosity_dComponentConcentration[c] = 0.0;
  }
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PROPPANTSLURRYFLUID_HPP_ */
