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
 * @file BlackOilFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_BLACKOILFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_BLACKOILFLUID_HPP_

#include "constitutive/fluid/multifluid/blackOil/BlackOilFluidBase.hpp"
#include "constitutive/fluid/multifluid/blackOil/PVTOData.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "math/interpolation/Interpolation.hpp"

namespace geos
{

namespace constitutive
{

class BlackOilFluid : public BlackOilFluidBase
{
public:
  /// Number of components supported by the model
  static constexpr integer NC_BO = 3;
  /// Number of hydrocarbon components supported by the model
  static constexpr integer HNC_BO = NC_BO - 1;
  /// Number of phases supported by the model
  static constexpr integer NP_BO = 3;

  /**
   * @brief Constructor for the Black-Oil class
   * @param[in] name the name of the class
   * @param[in] parent the parent group registering the Black-Oil fluid
   */
  BlackOilFluid( string const & name, Group * const parent );

  static string catalogName() { return "BlackOilFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Kernel wrapper class for BlackOilFluid
   *        This kernel can be called on the GPU
   */
  class KernelWrapper final : public BlackOilFluidBase::KernelWrapper
  {
public:
    GEOS_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseProp::SliceType const phaseDensity,
                          PhaseProp::SliceType const phaseMassDensity,
                          PhaseProp::SliceType const phaseViscosity,
                          PhaseProp::SliceType const phaseEnthalpy,
                          PhaseProp::SliceType const phaseInternalEnergy,
                          PhaseComp::SliceType const phaseCompFraction,
                          FluidProp::SliceType const totalDensity ) const override;

    GEOS_HOST_DEVICE
    virtual void update( localIndex k,
                         localIndex q,
                         real64 pressure,
                         real64 temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class BlackOilFluid;

    /**
     * @brief Constructor for the class doing in-kernel Black-Oil fluid updates
     * @param[in] PVTO structure containing oil phase formation volume factor and viscosity information
     * @param[in] phaseTypes type of phases
     * @param[in] phaseOrder order of phases
     * @param[in] hydrocarbonPhaseOrder order of the hydrocarbon phases in the model
     * @param[in] surfacePhaseMassDensity surface phase mass densities provided by the user
     * @param[in] formationVolFractionTables gas formation volume table
     * @param[in] viscosityTables gas viscosity table
     * @param[in] waterParams parameters (Bo, visc) for the water phase
     * @param[in] componentMolarWeight component molecular weights
     * @param[in] useMass flag to decide whether we return mass or molar densities
     * @param[in] phaseFraction phase fractions (+ derivatives) in the cell
     * @param[in] phaseDensity phase mass/molar densities (+ derivatives) in the cell
     * @param[in] phaseMassDensity phase mass densities (+ derivatives) in the cell
     * @param[in] phaseViscosity phase viscosities (+ derivatives) in the cell
     * @param[in] phaseCompFraction phase component fractions (+ derivatives) in the cell
     * @param[in] totalDensity total density in the cell
     */
    KernelWrapper( PVTOData const & PVTO,
                   arrayView1d< integer const > phaseTypes,
                   arrayView1d< integer const > phaseOrder,
                   arrayView1d< integer const > hydrocarbonPhaseOrder,
                   arrayView1d< real64 const > surfacePhaseMassDensity,
                   arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                   arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                   BlackOilFluidBase::WaterParams const waterParams,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    /**
     * @brief Utility function computing mass/molar densities and viscosity (keeping derivatives)
     * @param[in] needDerivs flag to decide whether derivatives are computed or not
     * @param[in] pressure pressure in the cell
     * @param[in] composition component fractions in the cell
     * @param[in] phaseFrac phase fractions in the cell
     * @param[out] phaseDensity phase mass/molar densities (+ derivatives) in the cell
     * @param[out] phaseMassDensity phase mass densities (+ derivatives) in the cell
     * @param[out] phaseViscosity phase viscosities (+ derivatives) in the cell
     * @param[out] phaseMolecularWeight phase molecular weights in the cell
     * @param[out] dPhaseMolecularWeight derivative of phase molecular weights wrt pressure, temperature and comp fractions
     */
    GEOS_HOST_DEVICE
    void computeDensitiesViscosities( bool const needDerivs,
                                      real64 const pressure,
                                      real64 const composition[NC_BO],
                                      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
                                      PhaseProp::SliceType const & phaseDens,
                                      PhaseProp::SliceType const & phaseMassDens,
                                      PhaseProp::SliceType const & phaseVisc,
                                      real64 phaseMolecularWeight[NP_BO],
                                      real64 dPhaseMolecularWeight[NP_BO][NC_BO+2] ) const;

    /**
     * @brief Utility function to compute phase fractions and phase component fractions (keeping derivatives)
     * @param[in] needDerivs flag to decide whether derivatives are computed or not
     * @param[in] pressure pressure in the cell
     * @param[in] composition component fractions in the cell
     * @param[out] phaseFraction phase fractions (+ derivatives) in the cell
     * @param[out] phaseCompFraction phase component fractions (+ derivatives) in the cell
     */
    GEOS_HOST_DEVICE
    void computeEquilibrium( bool const needDerivs,
                             real64 const pressure,
                             real64 const composition[NC_BO],
                             PhaseProp::SliceType const & phaseFraction,
                             PhaseComp::SliceType const & phaseCompFraction ) const;

    /**
     * @brief Utility function to compute Rs as a function of the bubble-point pressure
     * @param[in] presBub bubble point pressure
     * @param[out] Rs ratio of volume of gas to the volume of oil at standard conditions
     * @param[out] dRs_dPres derivative of the ratio of volume of gas to the volume of oil at standard conditions wrt pressure
     */
    GEOS_HOST_DEVICE
    void computeRs( real64 const presBub,
                    real64 & Rs,
                    real64 & dRs_dPres ) const;

    /**
     * @brief Utility function to compute Bo and Visc (and derivatives) as a function of Rs in the saturated case
     * @param[in] Rs ratio of volume of gas to the volume of oil at standard conditions
     * @param[in] dRs_dPres derivative of the ratio of volume of gas to the volume of oil at standard conditions wrt pressure
     * @param[out] Bo oil formation volume factor
     * @param[out] dBo_dPres derivative of oil formation volume factor wrt pressure
     * @param[out] visc oil viscosity
     * @param[out] dVisc_dPres derivative of oil viscosity wrt pressure
     */
    GEOS_HOST_DEVICE
    void computeSaturatedBoViscosity( real64 const Rs,
                                      real64 const dRs_dPres,
                                      real64 & Bo,
                                      real64 & dBo_dPres,
                                      real64 & visc,
                                      real64 & dVisc_dPres )  const;

    /**
     * @brief Utility function to compute Bo and Visc (and derivatives) as a function of Rs in the undersaturated case
     * @param[in] needDerivs flag to decide whether derivatives are computed or not
     * @param[in] pres pressure in the cell
     * @param[in] Rs ratio of volume of gas to the volume of oil at standard conditions
     * @param[in] dRs_dComp derivatives of the ratio of volume of gas to the volume of oil at standard conditions wrt comp fractions
     * @param[out] Bo oil formation volume factor
     * @param[out] dBo_dPres derivative of oil formation volume factor wrt pressure
     * @param[out] dBo_dComp derivatives of oil formation volume factor wrt comp fractions
     * @param[out] visc oil viscosity
     * @param[out] dVisc_dPres derivative of oil viscosity wrt pressure
     * @param[out] dVisc_dComp derivatives of oil viscosity wrt comp fractions
     */
    GEOS_HOST_DEVICE
    void computeUndersaturatedBoViscosity( bool const needDerivs,
                                           real64 const pres,
                                           real64 const Rs,
                                           real64 const dRs_dComp[],
                                           real64 & Bo,
                                           real64 & dBo_dPres,
                                           real64 dBo_dComp[],
                                           real64 & visc,
                                           real64 & dVisc_dPres,
                                           real64 dVisc_dComp[] ) const;

    /**
     * @brief Utility function to compute Bo and Visc (no derivatives) as a function of Rs in the undersaturated case
     * @param[in] Rs ratio of volume of gas to the volume of oil at standard conditions
     * @param[in] pres pressure in the cell
     * @param[out] Bo oil formation volume factor
     * @param[out] visc oil viscosity
     */
    GEOS_HOST_DEVICE
    void computeUndersaturatedBoViscosity( real64 const Rs,
                                           real64 const pres,
                                           real64 & Bo,
                                           real64 & visc ) const;

    /**
     * @brief Utility function to compute the mass and molar densities as a function of Rs and Bo
     * @param[in] needDerivs flag to decide whether derivatives are computed or not
     * @param[in] useMass flag to decide whether we return mass or molar densities
     * @param[in] Rs ratio of volume of gas to the volume of oil at standard conditions
     * @param[in] dRs_dPres derivative of the ratio of volume of gas to the volume of oil at standard conditions wrt pressure
     * @param[in] dRs_dComp derivatives of the ratio of volume of gas to the volume of oil at standard conditions wrt comp fractions
     * @param[in] Bo oil formation volume factor
     * @param[in] dBo_dPres derivative of oil formation volume factor wrt pressure
     * @param[in] dBo_dComp derivatives of oil formation volume factor wrt comp fractions
     * @param[out] dens oil density
     * @param[out] dDens_dC derivatives of oil density wrt pressure, temperature, and comp fractions
     */
    GEOS_HOST_DEVICE
    void computeMassMoleDensity( bool const needDerivs,
                                 bool const useMass,
                                 real64 const Rs,
                                 real64 const dRs_dPres,
                                 real64 const dRs_dComp[HNC_BO],
                                 real64 const Bo,
                                 real64 const dBo_dPres,
                                 real64 const dBo_dComp[HNC_BO],
                                 real64 & dens,
                                 arraySlice1d< real64, multifluid::USD_PHASE_DC - 3 > const & dDens ) const;

    /// Data needed to update the oil phase properties
    PVTOData::KernelWrapper m_PVTOView;
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:

  virtual void postInputInitialization() override;

private:

  virtual void readInputDataFromTableFunctions() override;

  virtual void readInputDataFromPVTFiles() override;

  /**
   * @brief Read all the PVT table provided by the user in Eclipse format
   * @param[in] oilTable the oil table data read from file
   * @param[in] oilSurfaceMassDensity the oil phase surface mass density
   * @param[in] oilSurfaceMolecularWeight the oil phase surface molecular weight
   * @param[in] gasSurfaceMassDensity the oil phase surface mass density
   * @param[in] gasSurfaceMolecularWeight the oil phase surface molecular weight
   */
  void fillPVTOData( array1d< array1d< real64 > > const & oilTable,
                     real64 oilSurfaceMassDensity,
                     real64 oilSurfaceMolecularWeight,
                     real64 gasSurfaceMassDensity,
                     real64 gasSurfaceMolecularWeight );

  /**
   * @brief Form the tables:
   *    - m_undersaturatedPressure
   *    - m_undersaturatedBO
   *    - m_undersaturatedVisc
   *  using the input data
   */
  void createUndersaturatedProperties();

  /**
   * @brief Extrapolate the tables:
   *    - m_undersaturatedPressure
   *    - m_undersaturatedBO
   *    - m_undersaturatedVisc
   *  up to the value of m_maxRelativePressure (if the data is missing for a given branch)
   */
  void extendUndersaturatedProperties();

  /**
   * @brief Refine the undersaturated tables and copy the data into the final 2D arrays that will be used in the kernel
   */
  void refineUndersaturatedTables( integer numRefinedPresPoints );

  /**
   * @brief Check the monotonicity of the PVTO table values
   */
  void checkTableConsistency() const;

  /// The data needed to compute the oil phase properties
  PVTOData m_PVTO;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           PhaseProp::SliceType const phaseFraction,
           PhaseProp::SliceType const phaseDensity,
           PhaseProp::SliceType const phaseMassDensity,
           PhaseProp::SliceType const phaseViscosity,
           PhaseProp::SliceType const phaseEnthalpy,
           PhaseProp::SliceType const phaseInternalEnergy,
           PhaseComp::SliceType const phaseCompFraction,
           FluidProp::SliceType const totalDensity ) const
{
  GEOS_UNUSED_VAR( temperature, phaseEnthalpy, phaseInternalEnergy );

  real64 compMoleFrac[NC_BO]{};
  real64 dCompMoleFrac_dCompMassFrac[NC_BO][NC_BO]{};

  real64 phaseMolecularWeight[NP_BO]{};
  real64 dPhaseMolecularWeight[NP_BO][NC_BO+2]{};

  // 1. Convert to moles if necessary

  if( m_useMass )
  {
    convertToMoleFractions( composition,
                            compMoleFrac,
                            dCompMoleFrac_dCompMassFrac );
  }
  else
  {
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Compute phase fractions and phase component fractions

  computeEquilibrium( true, // we want derivatives
                      pressure,
                      compMoleFrac,
                      phaseFraction,
                      phaseCompFraction );

  // 3. Compute phase densities and viscosities

  computeDensitiesViscosities( true, // we want derivatives
                               pressure,
                               compMoleFrac,
                               phaseFraction.value,
                               phaseDensity,
                               phaseMassDensity,
                               phaseViscosity,
                               phaseMolecularWeight,
                               dPhaseMolecularWeight );

  // 4. If mass variables used instead of molar, perform the conversion

  if( m_useMass )
  {
    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight,
                            phaseFraction,
                            phaseCompFraction,
                            phaseDensity.derivs,
                            phaseViscosity.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );
  }

  // 5. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );

}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  computeEquilibrium( bool const needDerivs,
                      real64 const pressure,
                      real64 const composition[NC_BO],
                      PhaseProp::SliceType const & phaseFraction,
                      PhaseComp::SliceType const & phaseCompFraction ) const
{
  using Deriv = multifluid::DerivativeOffset;
  using PT = BlackOilFluid::PhaseType;

  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];
  integer const ipWater = m_phaseOrder[PT::WATER];

  integer const icOil   = ipOil;
  integer const icGas   = ipGas;
  integer const icWater = ipWater;

  real64 const zo = composition[icOil];
  real64 const zg = composition[icGas];
  real64 const zw = composition[icWater];

  // 1. Make everything zero first

  auto setZero = []( real64 & val ){ val = 0.0; };
  LvArray::forValuesInSlice( phaseFraction.value, setZero );
  LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.value, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );

  // 2. Check feed first, and if only water is present (e.g., water inj), then skip

  if( zw >= 1.0 - MultiFluidConstants::minForSpeciesPresence )
  {
    phaseFraction.value[ipWater] = zw;
    if( needDerivs )
    {
      phaseFraction.derivs[ipWater][Deriv::dC+icWater] = 1.0;
    }
    phaseCompFraction.value[ipWater][icWater] = 1.0;
    return;
  }

  // 3. Compute Rs (saturated) as a function of pressure

  // oil
  real64 const & oilSurfaceMoleDensity = m_PVTOView.m_surfaceMoleDensity[PT::OIL];
  real64 RsSat = 0.0;
  real64 dRsSat_dP = 0.0;
  computeRs( pressure, RsSat, dRsSat_dP );
  if( RsSat < MultiFluidConstants::minForSpeciesPresence )
  {
    RsSat = MultiFluidConstants::minForSpeciesPresence;
  }

  // gas
  real64 const & gasSurfaceMoleDensity = m_PVTOView.m_surfaceMoleDensity[PT::GAS];

  // 4. Use the saturated Rs and density to compute the gas phase fraction
  //    Note : we assume the oil component cannot enter the gas phase

  real64 const Kg = ( oilSurfaceMoleDensity + gasSurfaceMoleDensity * RsSat ) / ( RsSat * gasSurfaceMoleDensity );
  real64 const dKg_dP = -oilSurfaceMoleDensity / gasSurfaceMoleDensity * dRsSat_dP / ( RsSat * RsSat );
  real64 const gasPhaseFraction = zo / ( 1.0 - Kg ) + zg;
  real64 const dGasPhaseFraction_dP = zo / ( ( 1.0 - Kg ) * ( 1.0 - Kg ) ) *  dKg_dP;
  real64 const dGasPhaseFraction_dzo = 1.0 / ( 1.0 - Kg );
  real64 const dGasPhaseFraction_dzg = 1.0;

  // 4. Update phase fraction and phase component fractions

  // 4.1 The gas phase is present
  if( ( gasPhaseFraction > 0.0 ) && ( gasPhaseFraction < 1.0 ) )
  {

    // phase fractions
    phaseFraction.value[ipOil] = 1.0 - gasPhaseFraction - zw;
    phaseFraction.value[ipGas] = gasPhaseFraction;
    phaseFraction.value[ipWater] =  zw;

    if( needDerivs )
    {
      phaseFraction.derivs[ipOil][Deriv::dP] = -dGasPhaseFraction_dP;
      phaseFraction.derivs[ipGas][Deriv::dP] = dGasPhaseFraction_dP;
      phaseFraction.derivs[ipOil][Deriv::dC+icOil] = -dGasPhaseFraction_dzo;
      phaseFraction.derivs[ipOil][Deriv::dC+icGas] = -dGasPhaseFraction_dzg;
      phaseFraction.derivs[ipOil][Deriv::dC+icWater] = -1.0;
      phaseFraction.derivs[ipGas][Deriv::dC+icOil] = dGasPhaseFraction_dzo;
      phaseFraction.derivs[ipGas][Deriv::dC+icGas] = dGasPhaseFraction_dzg;
      phaseFraction.derivs[ipWater][Deriv::dC+icWater] = 1.0;
    }

    // oil
    real64 const tmp = ( oilSurfaceMoleDensity + gasSurfaceMoleDensity * RsSat );
    real64 const tmpOil = oilSurfaceMoleDensity / tmp;
    real64 const dTmpOil_dP = -oilSurfaceMoleDensity * gasSurfaceMoleDensity * dRsSat_dP / ( tmp * tmp );
    phaseCompFraction.value[ipOil][icOil] = tmpOil;
    phaseCompFraction.value[ipOil][icGas] = 1.0 - tmpOil;
    phaseCompFraction.value[ipOil][icWater] = 0.0;

    if( needDerivs )
    {
      phaseCompFraction.derivs[ipOil][icOil][Deriv::dP] = dTmpOil_dP;
      phaseCompFraction.derivs[ipOil][icGas][Deriv::dP] = -dTmpOil_dP;
    }

    // gas
    real64 const tmpGas = gasSurfaceMoleDensity / ( gasSurfaceMoleDensity );
    phaseCompFraction.value[ipGas][icOil] = 1.0 - tmpGas;
    phaseCompFraction.value[ipGas][icGas] = tmpGas;
    phaseCompFraction.value[ipGas][icWater] = 0.0;

    // water
    phaseCompFraction.value[ipWater][icOil] = 0.0;
    phaseCompFraction.value[ipWater][icGas] = 0.0;
    phaseCompFraction.value[ipWater][icWater] = 1.0;

  }
  // 4.2 The gas phase is absent
  else
  {

    // phase fractions
    phaseFraction.value[ipOil] = 1.0 - zw;
    phaseFraction.value[ipGas] = 0.0;
    phaseFraction.value[ipWater] = zw;

    // oil
    phaseCompFraction.value[ipOil][icOil] = zo / ( 1 - zw );
    phaseCompFraction.value[ipOil][icGas] = zg / ( 1 - zw );
    phaseCompFraction.value[ipOil][icWater] = 0.0;

    // gas
    phaseCompFraction.value[ipWater][icOil] = 0.0;
    phaseCompFraction.value[ipWater][icGas] = 0.0;
    phaseCompFraction.value[ipWater][icWater] = 1.0;

    if( needDerivs )
    {
      phaseFraction.derivs[ipOil][Deriv::dC+icWater] = -1.0;
      phaseFraction.derivs[ipWater][Deriv::dC+icWater] = 1.0;
      phaseCompFraction.derivs[ipOil][icOil][Deriv::dC+icOil] = 1 / ( 1 - zw );
      phaseCompFraction.derivs[ipOil][icOil][Deriv::dC+icWater] = zo / (( 1 - zw )*( 1 - zw ));
      phaseCompFraction.derivs[ipOil][icGas][Deriv::dC+icGas] = 1 / ( 1 - zw );
      phaseCompFraction.derivs[ipOil][icGas][Deriv::dC+icWater] = zg / (( 1 - zw )*( 1 - zw ));
    }
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  computeDensitiesViscosities( bool const needDerivs,
                               real64 const pressure,
                               real64 const composition[NC_BO],
                               arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
                               PhaseProp::SliceType const & phaseDens,
                               PhaseProp::SliceType const & phaseMassDens,
                               PhaseProp::SliceType const & phaseVisc,
                               real64 phaseMolecularWeight[NP_BO],
                               real64 dPhaseMolecularWeight[NP_BO][NC_BO+2] ) const
{
  using Deriv = multifluid::DerivativeOffset;
  using PT = BlackOilFluid::PhaseType;

  integer const ipOil = m_phaseOrder[PT::OIL];
  integer const ipGas = m_phaseOrder[PT::GAS];
  integer const ipWater = m_phaseOrder[PT::WATER];

  integer const icOil = ipOil;
  integer const icGas = ipGas;

  bool const isWater = (ipWater >= 0 && phaseFrac[ipWater] > 0);
  bool const isGas = (ipGas >= 0 && phaseFrac[ipGas] > 0);
  bool const isOil = (ipOil >= 0 && phaseFrac[ipOil] > 0);

  auto setZero = []( real64 & val ){ val = 0.0; };
  LvArray::forValuesInSlice( phaseMassDens.value, setZero );
  LvArray::forValuesInSlice( phaseMassDens.derivs, setZero );
  LvArray::forValuesInSlice( phaseDens.value, setZero );
  LvArray::forValuesInSlice( phaseDens.derivs, setZero );
  LvArray::forValuesInSlice( phaseVisc.value, setZero );
  LvArray::forValuesInSlice( phaseVisc.derivs, setZero );

  // 1. Gas phase: look up in the formation vol factor tables

  if( isGas )
  {
    // interpolate in the table to get the phase formation vol factor and its derivative wrt pressure
    real64 derivative = 0.0;
    real64 const fvf = m_formationVolFactorTables[0].compute( &pressure, &derivative );

    // we are ready to update the densities
    real64 const fvfInv = 1.0 / fvf;

    phaseMassDens.value[ipGas] = m_surfacePhaseMassDensity[ipGas] * fvfInv;
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ipGas];
    phaseDens.value[ipGas] = phaseMassDens.value[ipGas] * mult;
    phaseMolecularWeight[ipGas] = m_componentMolarWeight[ipGas];

    if( needDerivs )
    {
      phaseMassDens.derivs[ipGas][Deriv::dP] = -derivative * phaseMassDens.value[ipGas] * fvfInv;
      phaseDens.derivs[ipGas][Deriv::dP] = phaseMassDens.derivs[ipGas][Deriv::dP] * mult;
      dPhaseMolecularWeight[ipGas][Deriv::dP] = 0.0;
      for( integer ic = 0; ic < NC_BO; ic++ )
      {
        phaseMassDens.derivs[ipGas][Deriv::dC+ic]  = 0.0;
        phaseDens.derivs[ipGas][Deriv::dC+ic]      = 0.0;
        dPhaseMolecularWeight[ipGas][Deriv::dC+ic] = 0.0;
      }
    }

    phaseVisc.value[ipGas] = m_viscosityTables[0].compute( &pressure, &(phaseVisc.derivs)[ipGas][Deriv::dP] );
  }

  // 2. Water phase: use the constant formation volume factor and compressibility provided by the user

  if( isWater )
  {
    // if water is present
    real64 const expCompDeltaPres = std::exp( -m_waterParams.compressibility * ( pressure - m_waterParams.referencePressure ) );
    real64 const dExpCompDeltaPres_dPres = -m_waterParams.compressibility * expCompDeltaPres;
    real64 const denom = m_waterParams.formationVolFactor * expCompDeltaPres;
    real64 const dDenom_dPres = m_waterParams.formationVolFactor * dExpCompDeltaPres_dPres;
    real64 const denomInv = 1.0 / denom;
    phaseMassDens.value[ipWater] = m_surfacePhaseMassDensity[ipWater] * denomInv;
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ipWater];
    phaseDens.value[ipWater] = phaseMassDens.value[ipWater] * mult;
    phaseVisc.value[ipWater] = m_waterParams.viscosity;
    phaseMolecularWeight[ipWater] = m_componentMolarWeight[ipWater];

    if( needDerivs )
    {
      phaseMassDens.derivs[ipWater][Deriv::dP] = -dDenom_dPres * phaseMassDens.value[ipWater] * denomInv;
      phaseDens.derivs[ipWater][Deriv::dP] = phaseMassDens.derivs[ipWater][Deriv::dP] * mult;
      dPhaseMolecularWeight[ipWater][Deriv::dP] = 0.0;
      for( integer ic = 0; ic < NC_BO; ic++ )
      {
        phaseMassDens.derivs[ipWater][Deriv::dC+ic]  = 0.0;
        phaseDens.derivs[ipWater][Deriv::dC+ic]      = 0.0;
        dPhaseMolecularWeight[ipWater][Deriv::dC+ic] = 0.0;
      }
    }

  }

  // 3. Oil phase: make the distinction between saturated and unsaturated conditions

  if( isOil )
  {

    real64 Rs = 0.0;
    real64 dRs_dC[HNC_BO]{};
    real64 dRs_dP = 0.0;
    real64 Bo = 0.0;
    real64 dBo_dP = 0.0;
    real64 dBo_dC[HNC_BO]{};
    real64 visc = 0.0;
    real64 dVisc_dP = 0.0;
    real64 dVisc_dC[HNC_BO]{};

    // saturated conditions
    if( isGas )
    {

      // compute Rs as a function of pressure
      computeRs( pressure, Rs, dRs_dP );

      // compute saturated properties (Bo, viscosity) as a function of Rs
      computeSaturatedBoViscosity( Rs, dRs_dP, Bo, dBo_dP, visc, dVisc_dP );

    }
    // unsaturated conditions
    else
    {

      // compute Rs as a function of composition
      real64 const densRatio = m_PVTOView.m_surfaceMoleDensity[PT::OIL] / m_PVTOView.m_surfaceMoleDensity[PT::GAS];
      Rs = densRatio * composition[icGas] / composition[icOil];
      dRs_dC[PT::OIL] = -densRatio * composition[icGas] / (composition[icOil] * composition[icOil]);
      dRs_dC[PT::GAS] =  densRatio  / composition[icOil];

      // compute undersaturated properties (Bo, viscosity) by two-step interpolation in undersaturated tables
      // this part returns numerical derivatives
      computeUndersaturatedBoViscosity( needDerivs, pressure, Rs, dRs_dC, Bo, dBo_dP, dBo_dC,
                                        visc, dVisc_dP, dVisc_dC );

    }

    // compute densities
    computeMassMoleDensity( needDerivs, true, Rs, dRs_dP, dRs_dC, Bo, dBo_dP, dBo_dC,
                            phaseMassDens.value[ipOil], phaseMassDens.derivs[ipOil] );
    computeMassMoleDensity( needDerivs, false, Rs, dRs_dP, dRs_dC, Bo, dBo_dP, dBo_dC,
                            phaseDens.value[ipOil], phaseDens.derivs[ipOil] );

    phaseMolecularWeight[ipOil] = phaseMassDens.value[ipOil] / phaseDens.value[ipOil];
    real64 const tmp = 1. / ( phaseDens.value[ipOil] * phaseDens.value[ipOil] );
    dPhaseMolecularWeight[ipOil][Deriv::dP] =
      tmp * ( phaseMassDens.derivs[ipOil][Deriv::dP] * phaseDens.value[ipOil] - phaseDens.derivs[ipOil][Deriv::dP] * phaseMassDens.value[ipOil] );
    for( integer ic = 0; ic < NC_BO; ++ic )
    {
      dPhaseMolecularWeight[ipOil][Deriv::dC+ic] =
        tmp * ( phaseMassDens.derivs[ipOil][Deriv::dC+ic] * phaseDens.value[ipOil] - phaseDens.derivs[ipOil][Deriv::dC+ic] * phaseMassDens.value[ipOil] );
    }

    if( m_useMass )
    {
      phaseDens.value[ipOil] = phaseMassDens.value[ipOil];
      if( needDerivs )
      {
        phaseDens.derivs[ipOil][Deriv::dP] = phaseMassDens.derivs[ipOil][Deriv::dP];
        for( integer ic = 0; ic < NC_BO; ++ic )
        {
          phaseDens.derivs[ipOil][Deriv::dC+ic] = phaseMassDens.derivs[ipOil][Deriv::dC+ic];
        }
      }
    }

    // copy viscosity into the final array
    // TODO: skip this step
    phaseVisc.value[ipOil] = visc;
    if( needDerivs )
    {
      phaseVisc.derivs[ipOil][Deriv::dP] = dVisc_dP;
      phaseVisc.derivs[ipOil][Deriv::dC+icOil] = dVisc_dC[PT::OIL];
      phaseVisc.derivs[ipOil][Deriv::dC+icGas] = dVisc_dC[PT::GAS];
    }
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  computeRs( real64 const presBub,
             real64 & Rs,
             real64 & dRs_dPres ) const
{

  integer const idx = LvArray::sortedArrayManipulation::find( m_PVTOView.m_bubblePressure.begin(),
                                                              m_PVTOView.m_bubblePressure.size(),
                                                              presBub );
  integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( m_PVTOView.m_bubblePressure.size()-1 ) );
  integer const iLow = iUp-1;
  interpolation::linearInterpolation( presBub - m_PVTOView.m_bubblePressure[iLow], m_PVTOView.m_bubblePressure[iUp] - presBub,
                                      m_PVTOView.m_Rs[iLow], m_PVTOView.m_Rs[iUp],
                                      Rs, dRs_dPres );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  computeSaturatedBoViscosity( real64 const Rs,
                               real64 const dRs_dPres,
                               real64 & Bo,
                               real64 & dBo_dPres,
                               real64 & visc,
                               real64 & dVisc_dPres ) const
{
  arrayView1d< real64 const > const & RsVec = m_PVTOView.m_Rs;
  arrayView1d< real64 const > const & BoVec = m_PVTOView.m_saturatedBo;
  arrayView1d< real64 const > const & viscVec = m_PVTOView.m_saturatedViscosity;

  integer const idx = LvArray::sortedArrayManipulation::find( RsVec.begin(),
                                                              RsVec.size(),
                                                              Rs );
  integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( RsVec.size()-1 ) );
  integer const iLow = iUp-1;

  interpolation::linearInterpolation( Rs - RsVec[iLow], RsVec[iUp] - Rs,
                                      BoVec[iLow], BoVec[iUp],
                                      Bo, dBo_dPres );
  interpolation::linearInterpolation( Rs - RsVec[iLow], RsVec[iUp] - Rs,
                                      viscVec[iLow], viscVec[iUp],
                                      visc, dVisc_dPres );

  // chain rule
  dBo_dPres *= dRs_dPres;
  dVisc_dPres *= dRs_dPres;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  computeUndersaturatedBoViscosity( bool const needDerivs,
                                    real64 const P,
                                    real64 const Rs,
                                    real64 const dRs_dComp[],
                                    real64 & Bo,
                                    real64 & dBo_dPres,
                                    real64 dBo_dComp[],
                                    real64 & visc,
                                    real64 & dVisc_dPres,
                                    real64 dVisc_dComp[] ) const
{
  computeUndersaturatedBoViscosity( Rs, P, Bo, visc );

  if( needDerivs )
  {
    // numerical derivatives

    // 1. dPres
    real64 const eps = 1e-6;
    real64 const inv_eps = 1.0 / eps;
    real64 const P_eps = P + eps;
    real64 Bo_eps   = 0.0;
    real64 visc_eps = 0.0;
    computeUndersaturatedBoViscosity( Rs, P_eps, Bo_eps, visc_eps );
    dBo_dPres = (Bo_eps - Bo) * inv_eps;
    dVisc_dPres = (visc_eps - visc) * inv_eps;

    // 2. dRs
    real64 const Rs_eps = Rs + eps;
    computeUndersaturatedBoViscosity( Rs_eps, P, Bo_eps, visc_eps );
    real64 const dBo_dRs =   (Bo_eps - Bo) * inv_eps;
    real64 const dVisc_dRs = (visc_eps - visc) * inv_eps;

    // 3. chainrule to dComp
    for( integer i = 0; i < HNC_BO; ++i )
    {
      dBo_dComp[i] = dBo_dRs * dRs_dComp[i];
      dVisc_dComp[i] = dVisc_dRs * dRs_dComp[i];
    }
  }

}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  computeUndersaturatedBoViscosity( real64 const Rs,
                                    real64 const pres,
                                    real64 & Bo,
                                    real64 & visc ) const
{
  // Step 1: interpolate for presBub
  arrayView1d< real64 const > const RsVec = m_PVTOView.m_Rs;
  integer idx = LvArray::sortedArrayManipulation::find( RsVec.begin(),
                                                        RsVec.size(),
                                                        Rs );
  integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( RsVec.size()-1 ) );
  integer const iLow = iUp-1;

  real64 const presBub = interpolation::linearInterpolation( Rs - m_PVTOView.m_Rs[iLow], m_PVTOView.m_Rs[iUp] - Rs,
                                                             m_PVTOView.m_bubblePressure[iLow], m_PVTOView.m_bubblePressure[iUp] );
  real64 const deltaPres = pres - presBub;

  // Step 2: get indices in undersatured pressure table
  real64 const deltaRsUp = LvArray::math::abs( m_PVTOView.m_Rs[iUp] - Rs );
  real64 const deltaRsLow = LvArray::math::abs( Rs - m_PVTOView.m_Rs[iLow] );
  arraySlice1d< real64 const > const & presUp = m_PVTOView.m_undersaturatedPressure2d[iUp];
  arraySlice1d< real64 const > const & presLow = m_PVTOView.m_undersaturatedPressure2d[iLow];

  idx = LvArray::sortedArrayManipulation::find( presUp.begin(),
                                                presUp.size(),
                                                deltaPres );
  integer const iUpP  = LvArray::math::min( LvArray::math::max( idx, 1 ), LvArray::integerConversion< integer >( presUp.size()-1 ) );
  integer const iLowP = iUpP-1;

  // Step 3: interpolate for Bo
  arraySlice1d< real64 const > const & BoUp = m_PVTOView.m_undersaturatedBo2d[iUp];
  arraySlice1d< real64 const > const & BoLow = m_PVTOView.m_undersaturatedBo2d[iLow];

  real64 const BoInterpLow = interpolation::linearInterpolation( deltaPres-presLow[iLowP], presLow[iUpP]-deltaPres,
                                                                 BoLow[iLowP], BoLow[iUpP] );
  real64 const BoInterpUp = interpolation::linearInterpolation( deltaPres-presUp[iLowP], presUp[iUpP]-deltaPres,
                                                                BoUp[iLowP], BoUp[iUpP] );
  Bo = interpolation::linearInterpolation( deltaRsLow, deltaRsUp, BoInterpLow, BoInterpUp );

  // Step 4: interpolate for viscosity
  arraySlice1d< real64 const > const & viscUp = m_PVTOView.m_undersaturatedViscosity2d[iUp];
  arraySlice1d< real64 const > const & viscLow = m_PVTOView.m_undersaturatedViscosity2d[iLow];

  real64 const viscInterpLow = interpolation::linearInterpolation( deltaPres-presLow[iLowP], presLow[iUpP]-deltaPres,
                                                                   viscLow[iLowP], viscLow[iUpP] );
  real64 const viscInterpUp = interpolation::linearInterpolation( deltaPres-presUp[iLowP], presUp[iUpP]-deltaPres,
                                                                  viscUp[iLowP], viscUp[iUpP] );
  visc = interpolation::linearInterpolation( deltaRsLow, deltaRsUp, viscInterpLow, viscInterpUp );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  computeMassMoleDensity( bool const needDerivs,
                          bool const useMass,
                          real64 const Rs,
                          real64 const dRs_dPres,
                          real64 const dRs_dComp[HNC_BO],
                          real64 const Bo,
                          real64 const dBo_dPres,
                          real64 const dBo_dComp[HNC_BO],
                          real64 & dens,
                          arraySlice1d< real64, multifluid::USD_PHASE_DC - 3 > const & dDens ) const
{
  using Deriv = multifluid::DerivativeOffset;
  using PT = BlackOilFluid::PhaseType;

  real64 const oilDens = (useMass)? m_PVTOView.m_surfaceMassDensity[PT::OIL]:
                         m_PVTOView.m_surfaceMoleDensity[PT::OIL];
  real64 const gasDens = (useMass)? m_PVTOView.m_surfaceMassDensity[PT::GAS]:
                         m_PVTOView.m_surfaceMoleDensity[PT::GAS];

  real64 const Binv = 1. / Bo;
  real64 const tmp = ( oilDens + gasDens * Rs );
  dens =  Binv * tmp;
  if( needDerivs )
  {
    integer const icOil = m_phaseOrder[PT::OIL];
    integer const icGas = m_phaseOrder[PT::GAS];

    dDens[Deriv::dP] = Binv * Binv * (Bo * gasDens * dRs_dPres - tmp * dBo_dPres);
    dDens[Deriv::dC+icOil] = Binv * Binv * (Bo * gasDens * dRs_dComp[PT::OIL] - tmp * dBo_dComp[PT::OIL]);
    dDens[Deriv::dC+icGas] = Binv * Binv * (Bo * gasDens * dRs_dComp[PT::GAS] - tmp * dBo_dComp[PT::GAS]);
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BlackOilFluid::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
