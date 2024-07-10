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
 * @file AquiferBoundaryCondition.hpp
 */


#ifndef GEOS_FIELDSPECIFICATION_AQUIFERBOUNDARYCONDITION_HPP
#define GEOS_FIELDSPECIFICATION_AQUIFERBOUNDARYCONDITION_HPP

#include "FieldSpecificationBase.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

/**
 * @class AquiferBoundaryCondition
 * Holds data and methods to apply a traction boundary condition
 */
class AquiferBoundaryCondition : public FieldSpecificationBase
{
public:

  /**
   * @class KernelWrapper
   *
   * A nested class encapsulating the kernel function doing the computing the average influx rate
   */
  class KernelWrapper
  {
public:

    /**
     * @brief Constructor of the kernel wrapper
     * @param[in] initialPressure the initial pressure in the aquifer
     * @param[in] density the water density in the aquifer
     * @param[in] gravCoef the elevation * gravVector in the aquifer
     * @param[in] timeConstant the time constant of the aquifer
     * @param[in] influxConstant the influx constant of the aquifer
     * @param[in] cumulativeFlux the cumulative flux of the aquifer
     * @param[in] pressureInfluenceFunction the pressure influence function of the aquifer
     */
    KernelWrapper( real64 initialPressure,
                   real64 density,
                   real64 gravCoef,
                   real64 timeConstant,
                   real64 influxConstant,
                   real64 cumulativeFlux,
                   TableFunction::KernelWrapper pressureInfluenceFunction )
      : m_initialPressure( initialPressure ),
      m_density( density ),
      m_gravCoef( gravCoef ),
      m_timeConstant( timeConstant ),
      m_influxConstant( influxConstant ),
      m_cumulativeFlux( cumulativeFlux ),
      m_pressureInfluenceFunction( pressureInfluenceFunction )
    {}

    /**
     * @brief Compute the aquifer-reservoir volumetric flux
     * @param[in] timeAtBeginningOfStep the time at the beginning of the step
     * @param[in] dt the time step size
     * @param[in] reservoirPressure the reservoir pressure
     * @param[in] reservoirPressure_n the reservoir pressure at the beginning of the time step
     * @param[in] reservoirGravCoef the elevation * gravVector in the aquifer
     * @param[in] areaFraction the area fraction for the face
     * @param[out] dAquiferVolFlux_dPres the derivative of the aquifer-reservoir volumetric flux
     * @return the aquifer-reservoir volumetric flux
     */
    GEOS_HOST_DEVICE
    inline real64
    compute( real64 const & timeAtBeginningOfStep,
             real64 const & dt,
             real64 const & reservoirPressure,
             real64 const & reservoirPressure_n,
             real64 const & reservoirGravCoef,
             real64 const & areaFraction,
             real64 & dAquiferVolFlux_dPres ) const;

private:

    // Physical parameters

    /// Aquifer initial pressure
    real64 m_initialPressure;

    /// Aquifer water density
    real64 m_density;

    /// Aquifer gravity coefficient
    real64 m_gravCoef;

    /// Aquifer time constant
    real64 m_timeConstant;

    /// Aquifer influx constant
    real64 m_influxConstant;

    /// Aquifer cumulative influx
    real64 m_cumulativeFlux;

    /// Pressure influence function
    TableFunction::KernelWrapper m_pressureInfluenceFunction;

  };

  /// @copydoc FieldSpecificationBase(string const &, Group *)
  AquiferBoundaryCondition( string const & name, Group * parent );

  /// deleted default constructor
  AquiferBoundaryCondition() = delete;

  /// default destructor
  virtual ~AquiferBoundaryCondition() = default;

  /// deleted copy constructor
  AquiferBoundaryCondition( AquiferBoundaryCondition const & ) = delete;

  /// defaulted move constructor
  AquiferBoundaryCondition( AquiferBoundaryCondition && ) = default;

  /// deleted copy assignment operator
  AquiferBoundaryCondition & operator=( AquiferBoundaryCondition const & ) = delete;

  /// deleted move assignment operator
  AquiferBoundaryCondition & operator=( AquiferBoundaryCondition && ) = delete;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "Aquifer"; }

  /**
   * @brief Create the wrapper performing in-kernel aquifer flow rate computation
   * @return the kernel wrapper
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @brief Setter for the R1Tensor storing the gravity vector
   * @param[in] gravityVector the gravity vector
   */
  void setGravityVector( R1Tensor const & gravityVector );

  /**
   * @brief Increment the cumulative flux for this aquifer
   * @param[in] fluxIncrement the new fluxes multiplied by dt
   */
  void saveConvergedState( real64 const fluxIncrement ) { m_cumulativeFlux += fluxIncrement; }

  /**
   * @brief Setter for the water phase index
   * @param[in] waterPhaseIndex the value of the water phase index
   */
  void setWaterPhaseIndex( integer const waterPhaseIndex ) { m_waterPhaseIndex = waterPhaseIndex; }

  /**
   * @brief Getter for the water phase index
   * @return the value of the water phase index
   */
  integer getWaterPhaseIndex() const { return m_waterPhaseIndex; }

  /**
   * @brief Getter for the aquifer water phase density
   * @return the value of the water phase density
   */
  real64 const & getWaterPhaseDensity() const { return m_density; }

  /**
   * @brief Getter for the aquifer water phase composition
   * @return an array storing the water phase component fractions
   */
  arrayView1d< real64 const > getWaterPhaseComponentFraction() const { return m_phaseComponentFraction.toViewConst(); }

  /**
   * @brief Getter for the aquifer water phase component names
   * @return an array storing the water phase component names
   */
  arrayView1d< string const > getWaterPhaseComponentNames() const { return m_phaseComponentNames.toViewConst(); }

  /**
   * @brief Flag to allow all phases to flow into the aquifer
   * @return true if we allow all phases to flow into the aquifer, false otherwise
   */
  bool allowAllPhasesIntoAquifer() const { return m_allowAllPhasesIntoAquifer; }

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationBase::viewKeyStruct
  {

    // aquifer geological properties

    /// @return The key for porosity
    constexpr static char const * aquiferPorosityString() { return "aquiferPorosity"; }

    /// @return The key for permeability
    constexpr static char const * aquiferPermeabilityString() { return "aquiferPermeability"; }

    // aquifer fluid properties

    /// @return The key for initial pressure
    constexpr static char const * aquiferInitialPressureString() { return "aquiferInitialPressure"; }

    /// @return The key for viscosity
    constexpr static char const * aquiferWaterViscosityString() { return "aquiferWaterViscosity"; }

    /// @return The key for density
    constexpr static char const * aquiferWaterDensityString() { return "aquiferWaterDensity"; }

    /// @return The key for phase component fraction
    constexpr static char const * aquiferWaterPhaseComponentFractionString() { return "aquiferWaterPhaseComponentFraction"; }

    /// @return The key for phase component names
    constexpr static char const * aquiferWaterPhaseComponentNamesString() { return "aquiferWaterPhaseComponentNames"; }

    /// @return The key for total compressibility
    constexpr static char const * aquiferTotalCompressibilityString() { return "aquiferTotalCompressibility"; }

    /// @return The key for the flag deciding whether we allow all phases into aquifer or not
    constexpr static char const * allowAllPhasesIntoAquiferString() { return "allowAllPhasesIntoAquifer"; }

    // aquifer geometry

    /// @return The key for elevation
    constexpr static char const * aquiferElevationString() { return "aquiferElevation"; }

    /// @return The key for thickness
    constexpr static char const * aquiferThicknessString() { return "aquiferThickness"; }

    /// @return The key for inner radius
    constexpr static char const * aquiferInnerRadiusString() { return "aquiferInnerRadius"; }

    /// @return The key for angle
    constexpr static char const * aquiferAngleString() { return "aquiferAngle"; }

    // table influence function

    /// @return The key for the pressure influence function
    constexpr static char const * pressureInfluenceFunctionNameString() { return "pressureInfluenceFunctionName"; }

    // cumulative flux

    /// @return The key for the cumulative aquifer flux
    constexpr static char const * cumulativeFluxString() { return "cumulativeFlux"; }

  };


protected:

  virtual void postInputInitialization() override final;

private:

  /**
   * @brief Sets up the default pressure influence function from the Carter-Tracy model
   */
  void setupDefaultPressureInfluenceFunction();

  /**
   * @brief Compute the aquifer time constant as a function of water properties, and aquifer geology / geometry
   */
  void computeTimeConstant();

  /**
   * @brief Compute the aquifer influx constant as a function of aquifer geology and geometry
   */
  void computeInfluxConstant();


  // Physical parameters

  /// Gravity vector
  R1Tensor m_gravityVector;

  /// Porosity of the aquifer
  real64 m_porosity;

  /// Permeability of the aquifer
  real64 m_permeability;

  /// Initial pressure
  real64 m_initialPressure;

  /// Flag to allow all phases to flow into the aquifer
  integer m_allowAllPhasesIntoAquifer;

  /// Water phase index
  integer m_waterPhaseIndex;

  /// Water viscosity
  real64 m_viscosity;

  /// Water density
  real64 m_density;

  /// Water phase component fraction
  array1d< real64 > m_phaseComponentFraction;

  /// Water phase component names
  array1d< string > m_phaseComponentNames;

  /// Total compressibility (rock + water)
  real64 m_totalCompressibility;

  /// Aquifer elevation
  real64 m_elevation;

  /// Aquifer thickness
  real64 m_thickness;

  /// Aquifer inner radius
  real64 m_innerRadius;

  /// Aquifer angle
  real64 m_angle;

  /// Aquifer time constant
  real64 m_timeConstant;

  /// Aquifer influx constant
  real64 m_influxConstant;

  /// Cumulative aquifer flux
  real64 m_cumulativeFlux;

  /// Name of the pressure influence table
  string m_pressureInfluenceFunctionName;

};

GEOS_HOST_DEVICE
real64
AquiferBoundaryCondition::KernelWrapper::
  compute( real64 const & timeAtBeginningOfStep,
           real64 const & dt,
           real64 const & reservoirPressure,
           real64 const & reservoirPressure_n,
           real64 const & reservoirGravCoef,
           real64 const & areaFraction,
           real64 & dAquiferVolFlux_dPres ) const
{
  // compute the dimensionless time (equation 5.5 of the Eclipse TD)
  real64 const dimensionlessTimeAtBeginningOfStep = timeAtBeginningOfStep / m_timeConstant;
  real64 const dimensionlessTimeAtEndOfStep = ( timeAtBeginningOfStep + dt ) / m_timeConstant;

  // compute the pressure influence and its derivative wrt to dimensionless time
  real64 dPresInfluence_dTime = 0;
  real64 const presInfluence = m_pressureInfluenceFunction.compute( &dimensionlessTimeAtEndOfStep, &dPresInfluence_dTime );

  // compute the potential difference between the reservoir (old pressure) and the aquifer
  real64 const potDiff = m_initialPressure - reservoirPressure_n - m_density * ( m_gravCoef - reservoirGravCoef );

  // compute the a (equation 5.8 of the Eclipse TD)
  real64 const timeConstantInv = 1.0 / m_timeConstant;
  real64 const denom = presInfluence - dimensionlessTimeAtBeginningOfStep * dPresInfluence_dTime;
  real64 const a = timeConstantInv * ( m_influxConstant * potDiff - m_cumulativeFlux * dPresInfluence_dTime ) / denom;

  // compute the b (equation 5.9 of the Eclipse TD)
  real64 const b = timeConstantInv * m_influxConstant / denom;

  // compute the average inflow rate Q (equation 5.7 of the Eclipse TD)
  real64 const aquiferVolFlux =  areaFraction * ( a - b * ( reservoirPressure - reservoirPressure_n ) );
  dAquiferVolFlux_dPres = -areaFraction * b;

  return aquiferVolFlux;
}


} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_AQUIFERBOUNDARYCONDITION_HPP */
