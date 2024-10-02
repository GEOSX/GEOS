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
 * @file LohrenzBrayClarkViscosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_LOHRENZBRAYCLARKVISCOSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_LOHRENZBRAYCLARKVISCOSITY_HPP_

#include "FunctionBase.hpp"

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ModelParameters;

class LohrenzBrayClarkViscosityUpdate final : public FunctionBaseUpdate
{
public:
/**
 * @brief Mixing types for phase viscosity
 */
  enum class MixingType : integer
  {
    HERNING_ZIPPERER,
    WILKE,
    BROKAW
  };

public:
  LohrenzBrayClarkViscosityUpdate( MixingType const mixing_type,
                                   arrayView1d< real64 const > const & componentCriticalVolume );

  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 const & density,
                arraySlice1d< real64 const, USD2 > const & dDensity,
                real64 & viscosity,
                arraySlice1d< real64, USD2 > const & dViscosity,
                bool useMass ) const;

  GEOS_HOST_DEVICE
  void setMixingType( MixingType const mixing_type )
  {
    m_mixing_type = mixing_type;
  }

private:
  /**
   * @brief Estimate pure component properties at dilute-gas conditions
   * @details This estimates pure component properties at dilute-gas conditions (pressure near atmospheric) using
   *          Stiel and Thodos [1961] correlation: https://doi.org/10.1002/aic.690070416
   *          Dilute viscosity is solely temperature dependent
   *          Units are converted so componentViscosity is in centipoise to match original reference
   * @param[in] numComponents The number of components
   * @param[in] componentProperties Physical properties of the components
   * @param[in] temperature The temperature
   * @param[out] componentDiluteViscosity The compoment dilute viscosity
   * @param[out] dComponentDiluteViscosity_dTemperature The derivatives of compoment dilute viscosity w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  void computeComponentDiluteViscosity_StielThodos( integer const numComponents,
                                                    ComponentProperties::KernelWrapper const & componentProperties,
                                                    real64 const temperature,
                                                    arraySlice1d< real64 > const componentDiluteViscosity,
                                                    arraySlice1d< real64 > const dComponentDiluteViscosity_dTemperature ) const
  {
    for( integer ic = 0; ic < numComponents; ++ic )
    {
      real64 const criticalPressure = componentProperties.m_componentCriticalPressure[ic];
      real64 const criticalTemperature = componentProperties.m_componentCriticalTemperature[ic];
      real64 const molarWeight = componentProperties.m_componentMolarWeight[ic];

      real64 reducedTemperature = temperature / criticalTemperature;
      real64 inverseComponentChi = 0.0;
      real64 discardDerivative = 0.0;
      inverseChiParameter( criticalPressure, criticalTemperature, molarWeight, inverseComponentChi,
                           discardDerivative, discardDerivative, discardDerivative );

      if( molarWeight < 2.1e-3 ) // hydrogen correlation, Stiel & Thodos, 1961, Eq. 12
      {
        componentDiluteViscosity[ic] = 90.71e-5 * pow( 0.1375*temperature - 1.67, 0.625 );
        dComponentDiluteViscosity_dTemperature[ic] = 90.71e-5 * 0.625 * 0.1375 * pow( 0.1375*temperature - 1.67, -0.375 );
      }
      else if( reducedTemperature <= 1.5 ) // nonpolar gas correlation at low temp, Eq. 9
      {
        componentDiluteViscosity[ic] = 34e-5 * pow( reducedTemperature, 0.94 ) * inverseComponentChi;
        dComponentDiluteViscosity_dTemperature[ic] = 34e-5 * 0.94 * pow( reducedTemperature, -0.06 ) * inverseComponentChi / criticalTemperature;
      }
      else // nonpolar gas correlation at high temp, Eq. 10
      {
        componentDiluteViscosity[ic] = 17.78e-5 * pow( 4.58*reducedTemperature-1.67, 0.625 ) * inverseComponentChi;
        dComponentDiluteViscosity_dTemperature[ic] = 17.78e-5 * 4.58 * 0.625 * pow( 4.58*reducedTemperature-1.67, -0.375 ) * inverseComponentChi / criticalTemperature;
      }
    }
  }

  /**
   * @brief Estimate phase viscosity at dilute-gas conditions using Herning and Zipperer [1936]
   * @details This estimates the phase viscosity properties at dilute-gas conditions (pressure near atmospheric) using
   *          Herning and Zipperer [1936] mixing rule.
   *          Herning, F. and Zipperer, L,: “Calculation of the Viscosity of Technical Gas Mixtures from the
   *          Viscosity of Individual Gases, german”, Gas u. Wasserfach (1936) 79, No. 49, 69.
   * @param[in] numComponents The number of components
   * @param[in] componentProperties Physical properties of the components
   * @param[in] temperature The temperature
   * @param[in] phaseComposition The composition of the phase
   * @param[in] componentDiluteViscosity The compoment dilute viscosity
   * @param[in] dComponentDiluteViscosity_dTemperature The derivatives of compoment dilute viscosity w.r.t. temperature
   * @param[out] phaseViscosity The phase viscosity
   * @param[out] dPhaseViscosity The phase viscosity derivatives
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void computePhaseDiluteViscosity_HerningZipperer( integer const numComponents,
                                                    ComponentProperties::KernelWrapper const & componentProperties,
                                                    real64 const temperature,
                                                    arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                                    arraySlice1d< real64 const > const & componentDiluteViscosity,
                                                    arraySlice1d< real64 const > const & dComponentDiluteViscosity_dTemperature,
                                                    real64 & phaseViscosity,
                                                    arraySlice1d< real64, USD2 > const & dPhaseViscosity ) const;

  /**
   * @brief Estimate phase viscosity at dilute-gas conditions using Wilke [1950]
   * @details This estimates the phase viscosity properties at dilute-gas conditions (pressure near atmospheric) using
   *          Wilke [1950] mixing rule. https://doi.org/10.1063/1.1747673
   * @param[in] numComponents The number of components
   * @param[in] componentProperties Physical properties of the components
   * @param[in] temperature The temperature
   * @param[in] phaseComposition The composition of the phase
   * @param[in] componentDiluteViscosity The compoment dilute viscosity
   * @param[in] dComponentDiluteViscosity_dTemperature The derivatives of compoment dilute viscosity w.r.t. temperature
   * @param[out] phaseViscosity The phase viscosity
   * @param[out] dPhaseViscosity The phase viscosity derivatives
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void computePhaseDiluteViscosity_Wilke( integer const numComponents,
                                          ComponentProperties::KernelWrapper const & componentProperties,
                                          real64 const temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          arraySlice1d< real64 const > const & componentDiluteViscosity,
                                          arraySlice1d< real64 const > const & dComponentDiluteViscosity_dTemperature,
                                          real64 & phaseViscosity,
                                          arraySlice1d< real64, USD2 > const & dPhaseViscosity ) const;

  /**
   * @brief Estimate phase viscosity at dilute-gas conditions using Brokaw[1968]
   * @details This estimates the phase viscosity properties at dilute-gas conditions (pressure near atmospheric) using
   *          Brokaw[1968] mixing rule.
   *          Brokaw, R. S. (1968). Viscosity of Gas Mixtures. United States: National Aeronautics and Space Administration.
   * @param[in] numComponents The number of components
   * @param[in] componentProperties Physical properties of the components
   * @param[in] temperature The temperature
   * @param[in] phaseComposition The composition of the phase
   * @param[in] componentDiluteViscosity The compoment dilute viscosity
   * @param[in] dComponentDiluteViscosity_dTemperature The derivatives of compoment dilute viscosity w.r.t. temperature
   * @param[out] phaseViscosity The phase viscosity
   * @param[out] dPhaseViscosity The phase viscosity derivatives
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void computePhaseDiluteViscosity_Brokaw( integer const numComponents,
                                           ComponentProperties::KernelWrapper const & componentProperties,
                                           real64 const temperature,
                                           arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                           arraySlice1d< real64 const > const & componentDiluteViscosity,
                                           arraySlice1d< real64 const > const & dComponentDiluteViscosity_dTemperature,
                                           real64 & phaseViscosity,
                                           arraySlice1d< real64, USD2 > const & dPhaseViscosity ) const;

  /**
   * @brief Estimates phase viscosity using Lohrenz, Bray & Clark [1964]
   * @details This estimates the phase viscosity at given (P,T) conditions using the Lohrenz, Bray & Clark [1964] correlation.
   *          This is an additional term added to the dilute gas estimate.
   *          https://doi.org/10.2118/915-PA
   * @param[in] numComponents The number of components
   * @param[in] componentProperties Physical properties of the components
   * @param[in] phaseComposition The composition of the phase
   * @param[in] phaseDensity The phase density
   * @param[in] dPhaseDensity The derivatives of the phase density
   * @param[out] phaseViscosity The phase viscosity
   * @param[out] dPhaseViscosity The phase viscosity derivatives
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void computePhaseViscosity_LohrenzBrayClark( integer const numComponents,
                                               ComponentProperties::KernelWrapper const & componentProperties,
                                               arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                               real64 const phaseDensity,
                                               arraySlice1d< real64 const, USD2 > const & dPhaseDensity,
                                               real64 & phaseViscosity,
                                               arraySlice1d< real64, USD2 > const & dPhaseViscosity ) const;

  /**
   * @brief Computes inverse chi parameter
   * @details Computes "1/chi" parameter (inverse of the viscosity-reducing parameter) from [ST 1961, LBC 1964].
   *          Using units of (K, atm, amu).
   * @param[in] criticalPressure The component critical pressure
   * @param[in] criticalTemperature The component critical temperature
   * @param[in] molarWeight The component molar weight
   * @param[out] value The inverse chi parameter
   * @param[out] derivP Derivative of the inverse chi parameter w.r.t. pressure
   * @param[out] derivT Derivative of the inverse chi parameter w.r.t. temperature
   * @param[out] derivM Derivative of the inverse chi parameter w.r.t. molar weight
   */
  GEOS_HOST_DEVICE
  void inverseChiParameter( real64 const criticalPressure,
                            real64 const criticalTemperature,
                            real64 const molarWeight,
                            real64 & value,
                            real64 & derivP,
                            real64 & derivT,
                            real64 & derivM ) const
  {
    real64 T  = pow( criticalTemperature, 1.0/6.0 );
    real64 dT = (1.0/6.0) * pow( criticalTemperature, -5.0/6.0 );

    real64 constexpr sqrt1000 = 31.6227766017; // note: kg/mol to atomic mass units
    real64 M  = sqrt1000 * sqrt( molarWeight );
    real64 dM = 0.5 * sqrt1000 / sqrt( molarWeight );

    real64 P  = pow( criticalPressure / PA_TO_ATM, 2.0/3.0 ); // note: pascal to atm conversion
    real64 dP = pow( PA_TO_ATM, -2.0/3.0 ) * pow( criticalPressure, -1.0/3.0 ) * 2.0/3.0;

    value  = M*P/T;
    derivP = M*dP/T;
    derivT = -M*P*dT/(T*T);
    derivM = dM*P/T;
  }

private:
  MixingType m_mixing_type;
  arrayView1d< real64 const > m_componentCriticalVolume;

private:
  // Conversion factor from cP to Pa.s
  static constexpr real64 CP_TO_PAS = 1.0e-3;
  // Conversion from Pa to atm
  static constexpr real64 PA_TO_ATM = 1.01325e+5;
};

class LohrenzBrayClarkViscosity : public FunctionBase
{
public:
  LohrenzBrayClarkViscosity( string const & name,
                             ComponentProperties const & componentProperties,
                             integer const phaseIndex,
                             ModelParameters const & modelParameters );

  static string catalogName() { return "LohrenzBrayClark"; }

  FunctionType functionType() const override
  {
    return FunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = LohrenzBrayClarkViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Parameters for the LBC viscosity model
  class Parameters : public ModelParameters
  {
public:
    Parameters( std::unique_ptr< ModelParameters > parameters );
    ~Parameters() override = default;

    string m_componentMixingType;

private:
    void registerParametersImpl( MultiFluidBase * fluid ) override;
    void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

    struct viewKeyStruct
    {
      static constexpr char const * componentMixingTypeString() { return "viscosityMixingRule"; }
    };
  };

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  ModelParameters const & m_parameters;
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LohrenzBrayClarkViscosityUpdate::MixingType,
              "HerningZipperer",
              "Wilke",
              "Brokaw" );

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

// Implementation
#include "LohrenzBrayClarkViscosityImpl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_LOHRENZBRAYCLARKVISCOSITY_HPP_
