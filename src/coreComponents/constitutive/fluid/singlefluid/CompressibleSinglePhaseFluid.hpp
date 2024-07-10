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
 * @file CompressibleSinglePhaseFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Update class for the model suitable for lambda capture.
 * @tparam DENS_EAT type of density exponent approximation
 * @tparam VISC_EAT type of viscosity exponent approximation
 */
template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT >
class CompressibleSinglePhaseUpdate : public SingleFluidBaseUpdate
{
public:

  using DensRelationType  = ExponentialRelation< real64, DENS_EAT >;
  using ViscRelationType  = ExponentialRelation< real64, VISC_EAT >;

  CompressibleSinglePhaseUpdate( DensRelationType const & densRelation,
                                 ViscRelationType const & viscRelation,
                                 arrayView2d< real64 > const & density,
                                 arrayView2d< real64 > const & dDens_dPres,
                                 arrayView2d< real64 > const & viscosity,
                                 arrayView2d< real64 > const & dVisc_dPres )
    : SingleFluidBaseUpdate( density,
                             dDens_dPres,
                             viscosity,
                             dVisc_dPres ),
    m_densRelation( densRelation ),
    m_viscRelation( viscRelation )
  {}

  /// Default copy constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate && ) = default;

  /// Deleted copy assignment operator
  CompressibleSinglePhaseUpdate & operator=( CompressibleSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  CompressibleSinglePhaseUpdate & operator=( CompressibleSinglePhaseUpdate && ) = delete;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const override
  {
    m_densRelation.compute( pressure, density, dDensity_dPressure );
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 const GEOS_UNUSED_PARAM( temperature ),
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & GEOS_UNUSED_PARAM( dDensity_dTemperature ),
                        real64 & viscosity,
                        real64 & dViscosity_dPressure,
                        real64 & GEOS_UNUSED_PARAM( dViscosity_dTemperature ),
                        real64 & GEOS_UNUSED_PARAM( internalEnergy ),
                        real64 & GEOS_UNUSED_PARAM( dInternalEnergy_dPressure ),
                        real64 & GEOS_UNUSED_PARAM( dInternalEnergy_dTemperature ),
                        real64 & GEOS_UNUSED_PARAM( enthalpy ),
                        real64 & GEOS_UNUSED_PARAM( dEnthalpy_dPressure ),
                        real64 & GEOS_UNUSED_PARAM( dEnthalpy_dTemperature ) ) const override
  {
    m_densRelation.compute( pressure, density, dDensity_dPressure );
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
  }


  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure ) const override
  {
    compute( pressure,
             m_density[k][q],
             m_dDens_dPres[k][q],
             m_viscosity[k][q],
             m_dVisc_dPres[k][q] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const GEOS_UNUSED_PARAM( temperature ) ) const override
  {
    compute( pressure,
             m_density[k][q],
             m_dDens_dPres[k][q],
             m_viscosity[k][q],
             m_dVisc_dPres[k][q] );
  }

private:

  /// Relationship between the fluid density and pressure
  DensRelationType m_densRelation;

  /// Relationship between the fluid viscosity and pressure
  ViscRelationType m_viscRelation;

};

class CompressibleSinglePhaseFluid : public SingleFluidBase
{
public:

  CompressibleSinglePhaseFluid( string const & name, Group * const parent );

  virtual ~CompressibleSinglePhaseFluid() override;

  static string catalogName() { return "CompressibleSinglePhaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using KernelWrapper = CompressibleSinglePhaseUpdate< ExponentApproximationType::Linear, ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct
  {
    static constexpr char const * defaultDensityString() { return "defaultDensity"; }
    static constexpr char const * defaultViscosityString() { return "defaultViscosity"; }
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * viscosibilityString() { return "viscosibility"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referenceDensityString() { return "referenceDensity"; }
    static constexpr char const * referenceViscosityString() { return "referenceViscosity"; }
    static constexpr char const * densityModelTypeString() { return "densityModelType"; }
    static constexpr char const * viscosityModelTypeString() { return "viscosityModelType"; }
  };

  real64 defaultDensity() const override final { return m_defaultDensity; }
  real64 defaultViscosity() const override final { return m_defaultViscosity; }

protected:

  virtual void postInputInitialization() override;

  /// default density value
  real64 m_defaultDensity;

  /// default viscosity value
  real64 m_defaultViscosity;

  /// scalar fluid bulk modulus parameter
  real64 m_compressibility;

  /// scalar fluid viscosity exponential coefficient
  real64 m_viscosibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// reference density parameter
  real64 m_referenceDensity;

  /// reference viscosity parameter
  real64 m_referenceViscosity;

  /// type of density model in terms of pressure
  ExponentApproximationType m_densityModelType;

  /// type of viscosity model
  ExponentApproximationType m_viscosityModelType;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
