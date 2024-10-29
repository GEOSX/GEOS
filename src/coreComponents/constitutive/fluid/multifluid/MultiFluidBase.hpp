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
 * @file MultiFluidBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{
namespace constitutive
{

class MultiFluidBase : public ConstitutiveBase
{
public:

  MultiFluidBase( string const & name,
                  Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** MultiFluid-specific interface

  /**
   * @brief Maximum supported number of fluid components (species)
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_COMPONENTS = MultiFluidConstants::MAX_NUM_COMPONENTS;

  /**
   * @brief Maximum supported number of fluid phases
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_PHASES = MultiFluidConstants::MAX_NUM_PHASES;

  /**
   * @return number of fluid components (species) in the model
   */
  integer numFluidComponents() const { return LvArray::integerConversion< integer >( m_componentNames.size() ); }

  /**
   * @brief Getter for the fluid component names
   * @return an array storing the component names
   */
  arrayView1d< string const > componentNames() const { return m_componentNames; }

  /**
   * @brief Getter for the fluid component molar weights
   * @return an arrayView1d storing the component molar weights
   */
  arrayView1d< real64 const > componentMolarWeights() const { return m_componentMolarWeight; }

  /**
   * @return number of fluid phases in the model
   */
  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  /**
   * @brief Getter for the fluid phase names
   * @return an array storing the phase names
   */
  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  /**
   * @brief Getter for the water phase index
   * @return the water phase index
   */
  virtual integer getWaterPhaseIndex() const = 0;

  /**
   * @brief Get the thermal flag.
   * @return boolean value indicating whether the model can be used to assemble the energy balance equation or not
   * @detail if isThermal is true, the constitutive model compute the enthalpy and internal energy of the phase.
   *         This can be used to check the compatibility of the constitutive model with the solver
   */
  virtual bool isThermal() const { return false; }

  /**
   * @brief Get the mass flag.
   * @return boolean value indicating whether the model is using mass-based quantities (as opposed to mole-based)
   */
  bool getMassFlag() const { return m_useMass; }

  /**
   * @brief Set the mass flag.
   * @param flag boolean value indicating whether the model should use mass-based quantities (as opposed to mole-based)
   *
   * @note This affects both input (compositions) and output quantities. The flag should be set prior to calling
   * any compute or state update methods.
   */
  void setMassFlag( bool const flag ) { m_useMass = flag; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseFraction() const
  { return m_phaseFraction.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseFraction() const
  { return m_phaseFraction.derivs; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity() const
  { return m_phaseDensity.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity_n() const
  { return m_phaseDensity_n; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseDensity() const
  { return m_phaseDensity.derivs; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseMassDensity() const
  { return m_phaseMassDensity.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseMassDensity() const
  { return m_phaseMassDensity.derivs; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseViscosity() const
  { return m_phaseViscosity.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseViscosity() const
  { return m_phaseViscosity.derivs; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > phaseCompFraction() const
  { return m_phaseCompFraction.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > phaseCompFraction_n() const
  { return m_phaseCompFraction_n; }

  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > dPhaseCompFraction() const
  { return m_phaseCompFraction.derivs; }

  arrayView2d< real64 const, multifluid::USD_FLUID > totalDensity() const
  { return m_totalDensity.value; }

  arrayView2d< real64 const, multifluid::USD_FLUID > totalDensity_n() const
  { return m_totalDensity_n; }

  arrayView3d< real64 const, multifluid::USD_FLUID_DC > dTotalDensity() const
  { return m_totalDensity.derivs; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseEnthalpy() const
  { return m_phaseEnthalpy.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseEnthalpy_n() const
  { return m_phaseEnthalpy_n; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseEnthalpy() const
  { return m_phaseEnthalpy.derivs; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseInternalEnergy() const
  { return m_phaseInternalEnergy.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseInternalEnergy_n() const
  { return m_phaseInternalEnergy_n; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseInternalEnergy() const
  { return m_phaseInternalEnergy.derivs; }

  /**
   * @brief Initialize the model
   * @param[in] phaseVolFraction an array containing the initial phase volume fractions
   */
  virtual void initializeState() const;

  /**
   * @brief Save the phase densities, component fractions, enthalpies and internal energies (for accumulation)
   */
  virtual void saveConvergedState() const override;

  /**
   * @brief If m_checkPVTTablesRanges, Check if the input values are in the expected PVT tables bounds
   * @param pressure input pressure to check
   * @param temperature input temperature to check (in K)
   * @throw a SimulationError if one of the input values is out of bound.
   */
  virtual void checkTablesParameters( real64 pressure, real64 temperature ) const = 0;

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * componentNamesString() { return "componentNames"; }
    static constexpr char const * componentMolarWeightString() { return "componentMolarWeight"; }
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * useMassString() { return "useMass"; }
    static constexpr char const * checkPVTTablesRangesString() { return "checkPVTTablesRanges"; }
  };


public:

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
  using FluidProp = MultiFluidVar< real64, 2, multifluid::LAYOUT_FLUID, multifluid::LAYOUT_FLUID_DC >;

public:
  class KernelWrapper
  {
public:

    /// @cond DO_NOT_DOCUMENT
    /// We need these SMFs to avoid host-device errors with CUDA.
    KernelWrapper() = default;
    KernelWrapper( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper && ) = default;
    /// @endcond

    /**
     * @brief Get number of elements in this wrapper.
     * @return number of elements
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    localIndex numElems() const { return m_phaseFraction.value.size( 0 ); }

    /**
     * @brief Get number of gauss points per element.
     * @return number of gauss points per element
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    localIndex numGauss() const { return m_phaseFraction.value.size( 1 ); }

    /**
     * @brief Get number of fluid components.
     * @return number of components
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    integer numComponents() const { return LvArray::integerConversion< integer >( m_componentMolarWeight.size() ); }

    /**
     * @brief Get number of fluid phases.
     * @return number of phases
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    integer numPhases() const { return LvArray::integerConversion< integer >( m_phaseFraction.value.size( 2 ) ); }

    GEOS_HOST_DEVICE arrayView3d< real64 const, multifluid::USD_PHASE > phaseFraction() const
    { return m_phaseFraction.value; }

    GEOS_HOST_DEVICE arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity() const
    { return m_phaseDensity.value; }

    GEOS_HOST_DEVICE arrayView3d< real64 const, multifluid::USD_PHASE > phaseMassDensity() const
    { return m_phaseMassDensity.value; }

    GEOS_HOST_DEVICE arrayView3d< real64 const, multifluid::USD_PHASE > phaseViscosity() const
    { return m_phaseViscosity.value; }

    GEOS_HOST_DEVICE arrayView4d< real64 const, multifluid::USD_PHASE_COMP > phaseCompFraction() const
    { return m_phaseCompFraction.value; }

    GEOS_HOST_DEVICE arrayView2d< real64 const, multifluid::USD_FLUID > totalDensity() const
    { return m_totalDensity.value; }

    GEOS_HOST_DEVICE arrayView3d< real64 const, multifluid::USD_PHASE > phaseEnthalpy() const
    { return m_phaseEnthalpy.value; }

    GEOS_HOST_DEVICE arrayView3d< real64 const, multifluid::USD_PHASE > phaseInternalEnergy() const
    { return m_phaseInternalEnergy.value; }

    /**
     * @brief Compute function to update properties in a cell without returning derivatives.
     * @details This delegates the call to the fluid wrapper using the value and derivative function.
     *          This is used for initialisation and boundary conditions.
     * @param[in] fluidWrapper the actual fluid kernel
     * @param[in] pressure pressure in the cell
     * @param[in] temperature temperature in the cell
     * @param[in] composition mass/molar component fractions in the cell
     * @param[out] phaseFraction phase fractions in the cell
     * @param[out] phaseDensity phase mass/molar density in the cell
     * @param[out] phaseMassDensity phase mass density in the cell
     * @param[out] phaseViscosity phase viscosity in the cell
     * @param[out] phaseEnthalpy phase enthalpy in the cell
     * @param[out] phaseInternalEnergy phase internal energy in the cell
     * @param[out] phaseCompFraction phase component fraction in the cell
     * @param[out] totalDensity total mass/molar density in the cell
     */
    template< typename FLUIDWRAPPER >
    GEOS_HOST_DEVICE
    static void computeValues( FLUIDWRAPPER const fluidWrapper,
                               real64 const pressure,
                               real64 const temperature,
                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                               arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                               arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                               arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                               arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                               arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                               arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                               arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                               real64 & totalDensity );

protected:

    /**
     * @brief Constructor for the kernel wrapper
     * @param[in] componentMolarWeight the array of component molar weight
     * @param[in] useMass the flag to decide whether the solver works in units of mass or moles
     * @param[out] phaseFraction the array of phase fractions (+ derivatives)
     * @param[out] phaseDensity the array of phase densities (+ derivatives)
     * @param[out] phaseMassDensity the array of phase mass densities (+derivatives)
     * @param[out] phaseViscosity the array of phase viscosities (+derivatives)
     * @param[out] phaseEnthalpy the array of phase enthalpy (+derivatives)
     * @param[out] phaseInternalEnergy the array of phase internal energy (+derivatives)
     * @param[out] phaseCompFraction the array of phase component fractions (+derivatives)
     * @param[out] totalDensity the total density (+derivatives)
     */
    KernelWrapper( arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity )
      : m_componentMolarWeight( std::move( componentMolarWeight ) ),
      m_useMass( useMass ),
      m_phaseFraction( std::move( phaseFraction ) ),
      m_phaseDensity( std::move( phaseDensity ) ),
      m_phaseMassDensity( std::move( phaseMassDensity ) ),
      m_phaseViscosity( std::move( phaseViscosity ) ),
      m_phaseEnthalpy( std::move( phaseEnthalpy ) ),
      m_phaseInternalEnergy( std::move( phaseInternalEnergy ) ),
      m_phaseCompFraction( std::move( phaseCompFraction ) ),
      m_totalDensity( std::move( totalDensity ) )
    { }

    /**
     * @brief Utility function to convert mass fractions to mole fractions
     * @tparam maxNumComp the max number of components
     * @tparam OUT_ARRAY the type of array storing the component mole fractions
     * @param[in] composition the component mass fractions
     * @param[out] compMoleFrac the newly converted component mole fractions
     * @detail The template is needed because PVTPackage expects a std::vector
     */
    template< integer maxNumComp, typename OUT_ARRAY >
    GEOS_HOST_DEVICE
    void convertToMoleFractions( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition,
                                 OUT_ARRAY && compMoleFrac ) const;

    /**
     * @brief Utility function to convert mass fractions to mole fractions and keep derivatives
     * @tparam maxNumComp the max number of components
     * @tparam OUT_ARRAY the type of array storing the component mole fractions
     * @param[in] composition the component mass fractions
     * @param[in] componentMolarWeight the component molar weight
     * @param[out] compMoleFrac the newly converted component mole fractions
     * @param[out] dCompMoleFrac_dCompMassFrac the derivatives of the newly converted component mole fractions
     * @detail The template is needed because PVTPackage expects a std::vector
     */
    template< integer maxNumComp, typename OUT_ARRAY >
    GEOS_HOST_DEVICE
    void convertToMoleFractions( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition,
                                 OUT_ARRAY && compMoleFrac,
                                 real64 ( &dCompMoleFrac_dCompMassFrac )[maxNumComp][maxNumComp] ) const;

    /**
     * @brief Utility function to convert mole fractions to mass fractions
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] phaseMolecularWeight the phase molecular weight computed by the constitutive model
     * @param[inout] phaseFrac the phase fractions in moles that will be converted to mass
     * @param[inout] phaseCompFrac the phase component fractions in moles that will be converted to mass
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOS_HOST_DEVICE
    void convertToMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                 arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                                 arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const phaseCompFrac ) const;

    /**
     * @brief Utility function to convert mole fractions to mass fractions and keep derivatives
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] dCompMoleFrac_dCompMassFrac the derivatives of mole fractions wrt mass fractions
     * @param[in] phaseMolecularWeight the phase molecular weight computed by the constitutive model
     * @param[in] dPhaseMolecularWeight the derivatives of phase molecular weights wrt pressure, temperature, and comp fractions
     * @param[inout] phaseFrac the phase fractions in moles that will be converted to mass
     * @param[inout] phaseCompFrac the phase component fractions in moles that will be converted to mass
     * @param[inout] dPhaseDens the derivatives of phase densities wrt pressure, temperature, and comp fractions
     * @param[inout] dPhaseVisc the derivatives of phase viscosities wrt pressure, temperature, and comp fractions
     * @param[inout] dPhaseEnthalpy the derivatives of phase enthalpy wrt pressure, temperature, and comp fractions
     * @param[inout] dPhaseInternalEnergy the derivatives of phase internal energy wrt pressure, temperature, and comp fractions
     * @detail This function performs three conversions
     *    1) Conversion of phase mass fractions into phase mole fractions
     *    2) Conversion of phase component mass fractions into phase component mole fractions
     *    3) Conversion of derivatives wrt mass fractions into derivatives wrt mole fractions
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOS_HOST_DEVICE
    void convertToMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                                 real64 const (&phaseMolecularWeight)[maxNumPhase],
                                 real64 const (&dPhaseMolecularWeight)[maxNumPhase][maxNumComp+2],
                                 PhaseProp::SliceType const phaseFrac,
                                 PhaseComp::SliceType const phaseCompFrac,
                                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens,
                                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc,
                                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseEnthalpy,
                                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseInternalEnergy ) const;

    /**
     * @brief Utility function to compute the internal energy from pressure, enthalpy, and density
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] pressure the pressure in the cell
     * @param[in] phaseFrac the phase fractions
     * @param[in] phaseMassDens the phase mass densities
     * @param[out] phaseEnthalpy the phase enthalpies
     * @param[out] phaseInternalEnergy the phase internal energy
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOS_HOST_DEVICE
    void computeInternalEnergy( real64 const & pressure,
                                arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                                arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseMassDens,
                                arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseEnthalpy,
                                arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseInternalEnergy ) const;

    /**
     * @brief Utility function to compute the internal energy from pressure, enthalpy, and density and keep derivatives
     * @param[in] pressure the pressure in the cellx
     * @param[in] phaseFrac the phase fractions (+ derivatives)
     * @param[in] phaseMassDens the phase mass densities (+ derivatives)
     * @param[out] phaseEnthalpy the phase enthalpies (+ derivatives)
     * @param[out] phaseInternalEnergy the phase internal energy (+ derivatives)
     */
    GEOS_HOST_DEVICE
    void computeInternalEnergy( real64 const & pressure,
                                PhaseProp::SliceType const phaseFrac,
                                PhaseProp::SliceType const phaseMassDens,
                                PhaseProp::SliceType const phaseEnthalpy,
                                PhaseProp::SliceType const phaseInternalEnergy ) const;

    /**
     * @brief Utility function to convert mole fractions to mass fractions
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] phaseFrac the phase fractions properly converted
     * @param[in] phaseFrac the phase densities in mass or moles
     * @param[out] totalDens the total density
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOS_HOST_DEVICE
    void computeTotalDensity( arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseDens,
                              real64 & totalDens ) const;

    /**
     * @brief Utility function to convert mole fractions to mass fractions and keep derivatives
     * @param[in] phaseFrac the phase fractions properly converted (+ derivatives)
     * @param[in] phaseDens the phase densities in mass or moles (+ derivatives)
     * @param[out] totalDens the total density (+ derivatives)
     */
    GEOS_HOST_DEVICE
    void computeTotalDensity( PhaseProp::SliceType const phaseFrac,
                              PhaseProp::SliceType const phaseDens,
                              FluidProp::SliceType const totalDens ) const;


    /// View on the component molar weights
    arrayView1d< real64 const > m_componentMolarWeight;

    /// Flag to decide whether the solver writes mole or mass balance
    bool m_useMass;

    /// Views on the phase properties
    PhaseProp::ViewType m_phaseFraction;
    PhaseProp::ViewType m_phaseDensity;
    PhaseProp::ViewType m_phaseMassDensity;
    PhaseProp::ViewType m_phaseViscosity;
    PhaseProp::ViewType m_phaseEnthalpy;
    PhaseProp::ViewType m_phaseInternalEnergy;
    PhaseComp::ViewType m_phaseCompFraction;
    FluidProp::ViewType m_totalDensity;

public:
    /**
     * @brief Calculate the total fluid compressibility
     * @param i Element index
     * @param q Quadrature node index
     * @return The total fluid compressibility
     */
    GEOS_HOST_DEVICE
    real64 totalCompressibility( integer const i, integer const q ) const
    {
      real64 const totalFluidDensity = totalDensity()( i, q );
      real64 const dTotalFluidDensity_dP = m_totalDensity.derivs( i, q, multifluid::DerivativeOffset::dP );
      return 0.0 < totalFluidDensity ? dTotalFluidDensity_dP / totalFluidDensity : 0.0;
    }

    /**
     * @brief Extract the phase mole fractions for a phase
     * @param i Element index
     * @param q Quadrature node index
     * @param p Phase index
     * @param[out] moleFractions The calculated mole fractions
     */
    template< typename OUT_ARRAY >
    GEOS_HOST_DEVICE
    void phaseCompMoleFraction( integer const i,
                                integer const q,
                                integer const p,
                                OUT_ARRAY && moleFractions ) const
    {
      integer const numComponents = m_componentMolarWeight.size( 0 );
      for( integer ic = 0; ic < numComponents; ++ic )
      {
        moleFractions[ic] = m_phaseCompFraction.value( i, q, p, ic );
      }
      detail::convertToMoleFractions( numComponents,
                                      m_componentMolarWeight,
                                      moleFractions,
                                      moleFractions );
    }

private:

    /**
     * @brief Utility function to convert phase mole fractions to phase mass fractions and keep derivatives
     * @tparam maxNumDof the max number of dofs
     * @tparam maxNumPhase the max number of phases
     * @param[in] phaseMolecularWeight the phase molecular weight computed by the constitutive model
     * @param[in] dPhaseMolecularWeight the derivatives of phase molecular weights wrt pressure, temperature, and comp fractions
     * @param[inout] phaseFrac the phase fractions in moles that will be converted to mass
     */
    template< integer maxNumDof, integer maxNumPhase >
    GEOS_HOST_DEVICE
    void convertToPhaseMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                      real64 const (&dPhaseMolecularWeight)[maxNumPhase][maxNumDof],
                                      PhaseProp::SliceType const phaseFrac ) const;

    /**
     * @brief Utility function to convert phase component mole fractions to phase component mass fractions and keep derivatives
     * @tparam maxNumDof the max number of dofs
     * @tparam maxNumPhase the max number of phases
     * @param[in] phaseMolecularWeight the phase molecular weight computed by the constitutive model
     * @param[in] dPhaseMolecularWeight the derivatives of phase molecular weights wrt pressure, temperature, and comp fractions
     * @param[in] phaseFrac the phase fractions in the cell
     * @param[inout] phaseCompFrac the phase component fractions in moles that will be converted to mass
     */
    template< integer maxNumDof, integer maxNumPhase >
    GEOS_HOST_DEVICE
    void convertToPhaseCompMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                          real64 const (&dPhaseMolecularWeight)[maxNumPhase][maxNumDof],
                                          PhaseProp::SliceType const phaseFrac,
                                          PhaseComp::SliceType const phaseCompFrac ) const;

    /**
     * @brief Utility function to convert derivatives wrt mole fractions into derivatives wrt mass fractions
     * @tparam maxNumComp the max number of components
     * @param[in] dCompMoleFrac_dCompMassFrac the derivatives of mole fractions wrt mass fractions
     * @param[inout] phaseFrac the phase fractions in the cell
     * @param[inout] phaseCompFrac the phase component fractions in the cell
     * @param[inout] dPhaseDens the derivatives of phase densities wrt pressure, temperature, and comp fractions
     * @param[inout] dPhaseVisc the derivatives of phase viscosities wrt pressure, temperature, and comp fractions
     * @param[inout] dPhaseEnthalpy the derivatives of phase enthalpy wrt pressure, temperature, and comp fractions
     * @param[inout] dPhaseInternalEnergy the derivatives of phase internal energy wrt pressure, temperature, and comp fractions
     */
    template< integer maxNumComp >
    GEOS_HOST_DEVICE
    void computeDerivativesWrtMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                                             PhaseProp::SliceType const phaseFrac,
                                             PhaseComp::SliceType const phaseCompFrac,
                                             arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens,
                                             arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc,
                                             arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseEnthalpy,
                                             arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseInternalEnergy ) const;

    /**
     * @brief Main compute function to update properties in a cell with derivatives (used in Newton iterations)
     * @param[in] pressure pressure in the cell
     * @param[in] temperature temperature in the cell
     * @param[in] composition mass/molar component fractions in the cell
     * @param[out] phaseFraction phase fractions in the cell  (+ derivatives)
     * @param[out] phaseDensity phase mass/molar density in the cell (+ derivatives)
     * @param[out] phaseMassDensity phase mass density in the cell (+ derivatives)
     * @param[out] phaseViscosity phase viscosity in the cell (+ derivatives)
     * @param[out] phaseEnthalpy phase enthalpy in the cell (+ derivatives)
     * @param[out] phaseInternalEnergy phase internal energy in the cell (+ derivatives)
     * @param[out] phaseCompFraction phase component fraction in the cell (+ derivatives)
     * @param[out] totalDensity total mass/molar density in the cell (+ derivatives)
     */
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
                          FluidProp::SliceType const totalDensity ) const = 0;

    /**
     * @brief Update function seen by the solver and calling the compute function with appropriate parameters
     * @param[in] k index of the cell
     * @param[in] q index of the quadrature point
     * @param[in] pressure pressure in the cell
     * @param[in] temperature temperature in the cell
     * @param[in] composition mass/molar component fractions in the cell
     */
    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const = 0;
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

  // flag indicating whether input/output component fractions are treated as mass fractions
  int m_useMass;

  /// Enable an error when checkTableParameters() is called and the input pressure or temperature of the PVT tables is out of range
  integer m_checkPVTTablesRanges;

  // general fluid composition information

  array1d< string > m_componentNames;
  array1d< real64 > m_componentMolarWeight;
  array1d< string > m_phaseNames;

  // constitutive data

  PhaseProp m_phaseFraction;
  PhaseProp m_phaseDensity;
  PhaseProp m_phaseMassDensity;
  PhaseProp m_phaseViscosity;
  PhaseProp m_phaseEnthalpy;
  PhaseProp m_phaseInternalEnergy;
  PhaseComp m_phaseCompFraction;
  FluidProp m_totalDensity;

  // backup data

  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseDensity_n;
  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseEnthalpy_n;
  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseInternalEnergy_n;
  array4d< real64, multifluid::LAYOUT_PHASE_COMP > m_phaseCompFraction_n;
  array2d< real64, multifluid::LAYOUT_FLUID > m_totalDensity_n;

};

template< typename FLUIDWRAPPER >
GEOS_HOST_DEVICE
void
MultiFluidBase::KernelWrapper::computeValues( FLUIDWRAPPER const fluidWrapper,
                                              real64 const pressure,
                                              real64 const temperature,
                                              arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                                              arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                                              real64 & totalDensity )
{
  integer constexpr maxNumPhase = MAX_NUM_PHASES;
  integer constexpr maxNumComp = MAX_NUM_COMPONENTS;
  integer constexpr maxNumDof = MAX_NUM_COMPONENTS + 2;
  integer const numPhase = fluidWrapper.numPhases();
  integer const numComp = fluidWrapper.numComponents();

  // Allocate data for derivatives. All properties of the same type will use the same memory
  // space for the derivatives. The derivatives returned will clearly be garbage but the values
  // should be correct.
  StackArray< real64, 4, maxNumDof * maxNumPhase, multifluid::LAYOUT_PHASE_DC > dPhaseProp( 1, 1, numPhase, numComp+2 );
  StackArray< real64, 5, maxNumDof * maxNumComp * maxNumPhase, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComp( 1, 1, numPhase, numComp, numComp+2 );
  StackArray< real64, 3, maxNumDof, multifluid::LAYOUT_FLUID_DC > dFluidProp( 1, 1, numComp+2 );

  // Wrap the output in multi variable objects
  PhaseProp::SliceType phaseFractionWrapper { phaseFraction, dPhaseProp[0][0] };
  PhaseProp::SliceType phaseDensityWrapper { phaseDensity, dPhaseProp[0][0] };
  PhaseProp::SliceType phaseMassDensityWrapper { phaseMassDensity, dPhaseProp[0][0] };
  PhaseProp::SliceType phaseViscosityWrapper { phaseViscosity, dPhaseProp[0][0] };
  PhaseProp::SliceType phaseEnthalpyWrapper { phaseEnthalpy, dPhaseProp[0][0] };
  PhaseProp::SliceType phaseInternalEnergyWrapper { phaseInternalEnergy, dPhaseProp[0][0] };
  PhaseComp::SliceType phaseCompFractionWrapper { phaseCompFraction, dPhaseComp[0][0] };
  FluidProp::SliceType totalDensityWrapper { totalDensity, dFluidProp[0][0] };

  // Pass on to fluid kernel
  fluidWrapper.compute( pressure,
                        temperature,
                        composition,
                        phaseFractionWrapper,
                        phaseDensityWrapper,
                        phaseMassDensityWrapper,
                        phaseViscosityWrapper,
                        phaseEnthalpyWrapper,
                        phaseInternalEnergyWrapper,
                        phaseCompFractionWrapper,
                        totalDensityWrapper );
}

template< integer maxNumComp, typename OUT_ARRAY >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  convertToMoleFractions( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition,
                          OUT_ARRAY && compMoleFrac ) const
{
  real64 dCompMoleFrac_dCompMassFrac[maxNumComp][maxNumComp]{};

  convertToMoleFractions( composition,
                          compMoleFrac,
                          dCompMoleFrac_dCompMassFrac );
}

template< integer maxNumComp, typename OUT_ARRAY >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  convertToMoleFractions( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition,
                          OUT_ARRAY && compMoleFrac,
                          real64 (& dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp] ) const
{
  integer const numComps = numComponents();

  real64 totalMolality = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
    compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
    dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
    totalMolality += compMoleFrac[ic];
  }

  real64 const totalMolalityInv = 1.0 / totalMolality;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    compMoleFrac[ic] *= totalMolalityInv;

    for( integer jc = 0; jc < numComps; ++jc )
    {
      dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
      dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
    }
  }
}

template< integer maxNumComp, integer maxNumPhase >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  convertToMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const phaseCompFrac ) const
{
  using namespace multifluid;

  integer constexpr maxNumDof = maxNumComp + 2;
  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  real64 dCompMoleFrac_dCompMassFrac[maxNumComp][maxNumComp]{};
  real64 dPhaseMolecularWeight[maxNumPhase][maxNumDof]{};

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseFrac( 1, 1, numPhase, numComp+2 );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac, dPhaseFrac[0][0] };

  StackArray< real64, 5, maxNumDof *maxNumComp *maxNumPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFrac( 1, 1, numPhase, numComp, numComp+2 );
  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 >
  phaseCompFracAndDeriv { phaseCompFrac, dPhaseCompFrac[0][0] };

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseDens( 1, 1, numPhase, numComp+2 );
  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseVisc( 1, 1, numPhase, numComp+2 );
  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseEnthalpy( 1, 1, numPhase, numComp+2 );
  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseInternalEnergy( 1, 1, numPhase, numComp+2 );

  convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                          phaseMolecularWeight,
                          dPhaseMolecularWeight,
                          phaseFracAndDeriv,
                          phaseCompFracAndDeriv,
                          dPhaseDens[0][0],
                          dPhaseVisc[0][0],
                          dPhaseEnthalpy[0][0],
                          dPhaseInternalEnergy[0][0] );
}

template< integer maxNumComp, integer maxNumPhase >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  convertToMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                          real64 const (&phaseMolecularWeight)[maxNumPhase],
                          real64 const (&dPhaseMolecularWeight)[maxNumPhase][maxNumComp+2],
                          PhaseProp::SliceType const phaseFrac,
                          PhaseComp::SliceType const phaseCompFrac,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseEnthalpy,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseInternalEnergy ) const
{
  convertToPhaseMassFractions( phaseMolecularWeight,
                               dPhaseMolecularWeight,
                               phaseFrac );
  convertToPhaseCompMassFractions( phaseMolecularWeight,
                                   dPhaseMolecularWeight,
                                   phaseFrac,
                                   phaseCompFrac );
  computeDerivativesWrtMassFractions( dCompMoleFrac_dCompMassFrac,
                                      phaseFrac,
                                      phaseCompFrac,
                                      dPhaseDens,
                                      dPhaseVisc,
                                      dPhaseEnthalpy,
                                      dPhaseInternalEnergy );
}

template< integer maxNumDof, integer maxNumPhase >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  convertToPhaseMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                               real64 const (&dPhaseMolecularWeight)[maxNumPhase][maxNumDof],
                               PhaseProp::SliceType const phaseFrac ) const
{
  using Deriv = multifluid::DerivativeOffset;

  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  real64 totalMass{};
  real64 dTotalMass[maxNumDof]{};

  // 1. Compute mass of each phase and total mass (on a 1-mole basis)
  for( integer ip = 0; ip < numPhase; ++ip )
  {

    real64 const nu = phaseFrac.value[ip];

    phaseFrac.value[ip] *= phaseMolecularWeight[ip];
    phaseFrac.derivs[ip][Deriv::dP] = phaseFrac.derivs[ip][Deriv::dP] * phaseMolecularWeight[ip] + nu * dPhaseMolecularWeight[ip][Deriv::dP];
    phaseFrac.derivs[ip][Deriv::dT] = phaseFrac.derivs[ip][Deriv::dT] * phaseMolecularWeight[ip] + nu * dPhaseMolecularWeight[ip][Deriv::dT];

    totalMass += phaseFrac.value[ip];
    dTotalMass[Deriv::dP] += phaseFrac.derivs[ip][Deriv::dP];
    dTotalMass[Deriv::dT] += phaseFrac.derivs[ip][Deriv::dT];

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFrac.derivs[ip][Deriv::dC+jc] =
        phaseFrac.derivs[ip][Deriv::dC+jc] * phaseMolecularWeight[ip] + nu * dPhaseMolecularWeight[ip][Deriv::dC+jc];
      dTotalMass[Deriv::dC+jc] += phaseFrac.derivs[ip][Deriv::dC+jc];
    }
  }

  // 2. Normalize to get mass fractions
  real64 const totalMassInv = 1.0 / totalMass;
  for( integer ip = 0; ip < numPhase; ++ip )
  {
    phaseFrac.value[ip] *= totalMassInv;
    phaseFrac.derivs[ip][Deriv::dP] = ( phaseFrac.derivs[ip][Deriv::dP] - phaseFrac.value[ip] * dTotalMass[Deriv::dP] ) * totalMassInv;
    phaseFrac.derivs[ip][Deriv::dT] = ( phaseFrac.derivs[ip][Deriv::dT] - phaseFrac.value[ip] * dTotalMass[Deriv::dT] ) * totalMassInv;

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFrac.derivs[ip][Deriv::dC+jc] = ( phaseFrac.derivs[ip][Deriv::dC+jc] - phaseFrac.value[ip] * dTotalMass[Deriv::dC+jc] ) * totalMassInv;
    }
  }
}

template< integer maxNumDof, integer maxNumPhase >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  convertToPhaseCompMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                   real64 const (&dPhaseMolecularWeight)[maxNumPhase][maxNumDof],
                                   PhaseProp::SliceType const phaseFrac,
                                   PhaseComp::SliceType const phaseCompFrac ) const
{
  using Deriv = multifluid::DerivativeOffset;

  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  for( integer ip = 0; ip < numPhase; ++ip )
  {

    // Note: for Black-Oil, phaseMolecularWeight can be zero for absent gas
    // TODO: implement forExistingPhase lambda
    bool const phaseExists = (phaseFrac.value[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const phaseMolecularWeightInv = 1.0 / phaseMolecularWeight[ip];

    for( integer ic = 0; ic < numComp; ++ic )
    {
      real64 const compMW = m_componentMolarWeight[ic];

      phaseCompFrac.value[ip][ic] = phaseCompFrac.value[ip][ic] * compMW * phaseMolecularWeightInv;
      phaseCompFrac.derivs[ip][ic][Deriv::dP] =
        ( phaseCompFrac.derivs[ip][ic][Deriv::dP] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMolecularWeight[ip][Deriv::dP] ) * phaseMolecularWeightInv;
      phaseCompFrac.derivs[ip][ic][Deriv::dT] =
        ( phaseCompFrac.derivs[ip][ic][Deriv::dT] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMolecularWeight[ip][Deriv::dT] ) * phaseMolecularWeightInv;

      for( integer jc = 0; jc < numComp; ++jc )
      {
        phaseCompFrac.derivs[ip][ic][Deriv::dC+jc] =
          phaseMolecularWeightInv *
          ( phaseCompFrac.derivs[ip][ic][Deriv::dC+jc] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMolecularWeight[ip][Deriv::dC+jc] );
      }
    }
  }
}

template< integer maxNumComp >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  computeDerivativesWrtMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                                      PhaseProp::SliceType const phaseFrac,
                                      PhaseComp::SliceType const phaseCompFrac,
                                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens,
                                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc,
                                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseEnthalpy,
                                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseInternalEnergy ) const
{
  using Deriv = multifluid::DerivativeOffset;

  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  real64 work[maxNumComp]{};
  for( integer ip = 0; ip < numPhase; ++ip )
  {
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, phaseFrac.derivs[ip], work, Deriv::dC );
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, dPhaseDens[ip], work, Deriv::dC );
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, dPhaseVisc[ip], work, Deriv::dC );
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, dPhaseEnthalpy[ip], work, Deriv::dC );
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, dPhaseInternalEnergy[ip], work, Deriv::dC );

    for( integer ic = 0; ic < numComp; ++ic )
    {
      applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, phaseCompFrac.derivs[ip][ic], work, Deriv::dC );
    }
  }
}

template< integer maxNumComp, integer maxNumPhase >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  computeInternalEnergy( real64 const & pressure,
                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseMassDens,
                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseEnthalpy,
                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseInternalEnergy ) const
{
  using namespace multifluid;

  integer constexpr maxNumDof = maxNumComp + 2;
  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseFrac( 1, 1, numPhase, numComp+2 );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac, dPhaseFrac[0][0] };

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseMassDens( 1, 1, numPhase, numComp+2 );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseMassDensAndDeriv { phaseMassDens, dPhaseMassDens[0][0] };

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseEnthalpy( 1, 1, numPhase, numComp+2 );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseEnthalpyAndDeriv { phaseEnthalpy, dPhaseEnthalpy[0][0] };

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseInternalEnergy( 1, 1, numPhase, numComp+2 );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseInternalEnergyAndDeriv { phaseInternalEnergy, dPhaseInternalEnergy[0][0] };

  computeInternalEnergy( pressure,
                         phaseFracAndDeriv,
                         phaseMassDensAndDeriv,
                         phaseEnthalpyAndDeriv,
                         phaseInternalEnergyAndDeriv );
}

GEOS_HOST_DEVICE
inline void
MultiFluidBase::KernelWrapper::
  computeInternalEnergy( real64 const & pressure,
                         PhaseProp::SliceType const phaseFrac,
                         PhaseProp::SliceType const phaseMassDens,
                         PhaseProp::SliceType const phaseEnthalpy,
                         PhaseProp::SliceType const phaseInternalEnergy ) const
{
  integer const numPhase = numPhases();
  integer const numComp = numComponents();
  integer const numDOF = numComp + 2;
  for( integer ip = 0; ip < numPhase; ++ip )
  {

    bool const phaseExists = (phaseFrac.value[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const densInv = 1.0 / phaseMassDens.value[ip];
    real64 const densInvSquared = densInv * densInv;
    phaseInternalEnergy.value[ip] = phaseEnthalpy.value[ip] - pressure * densInv;
    for( integer idof = 0; idof < numDOF; ++idof )
    {
      phaseInternalEnergy.derivs[ip][idof] = phaseEnthalpy.derivs[ip][idof] + pressure * phaseMassDens.derivs[ip][idof] * densInvSquared;
    }
    phaseInternalEnergy.derivs[ip][multifluid::DerivativeOffset::dP] -= densInv;
  }
}

template< integer maxNumComp, integer maxNumPhase >
GEOS_HOST_DEVICE
inline void
MultiFluidBase::KernelWrapper::
  computeTotalDensity( arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                       arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseDens,
                       real64 & totalDens ) const
{
  using namespace multifluid;

  integer constexpr maxNumDof = maxNumComp + 2;
  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  StackArray< real64, 3, maxNumDof, LAYOUT_FLUID_DC > dTotalDens( 1, 1, numComp+2 );
  MultiFluidVarSlice< real64, 0, USD_FLUID - 2, USD_FLUID_DC - 2 > totalDensAndDeriv { totalDens, dTotalDens[0][0] };

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseFrac( 1, 1, numPhase, numComp+2 );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseFracAndDeriv { phaseFrac, dPhaseFrac[0][0] };

  StackArray< real64, 4, maxNumDof *maxNumPhase, LAYOUT_PHASE_DC > dPhaseDens( 1, 1, numPhase, numComp+2 );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseDensAndDeriv { phaseDens, dPhaseDens[0][0] };

  computeTotalDensity( phaseFracAndDeriv,
                       phaseDensAndDeriv,
                       totalDensAndDeriv );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
MultiFluidBase::KernelWrapper::
  computeTotalDensity( PhaseProp::SliceType const phaseFraction,
                       PhaseProp::SliceType const phaseDensity,
                       FluidProp::SliceType const totalDensity ) const
{
  using Deriv = multifluid::DerivativeOffset;

  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  totalDensity.value = 0.0;
  LvArray::forValuesInSlice( totalDensity.derivs, []( real64 & val ){ val = 0.0; } );

  // 1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
  for( integer ip = 0; ip < numPhase; ++ip )
  {
    bool const phaseExists = (phaseFraction.value[ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const densInv = 1.0 / phaseDensity.value[ip];
    real64 const value = phaseFraction.value[ip] * densInv;

    totalDensity.value += value;
    totalDensity.derivs[Deriv::dP] += ( phaseFraction.derivs[ip][Deriv::dP] - value * phaseDensity.derivs[ip][Deriv::dP] ) * densInv;
    totalDensity.derivs[Deriv::dT] += ( phaseFraction.derivs[ip][Deriv::dT] - value * phaseDensity.derivs[ip][Deriv::dT] ) * densInv;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity.derivs[Deriv::dC+ic] += ( phaseFraction.derivs[ip][Deriv::dC+ic] - value * phaseDensity.derivs[ip][Deriv::dC+ic] ) * densInv;
    }
  }

  // 2. Invert the previous quantity to get actual density
  totalDensity.value = 1.0 / totalDensity.value;
  real64 const minusDens2 = -totalDensity.value * totalDensity.value;
  totalDensity.derivs[Deriv::dP] *= minusDens2;
  totalDensity.derivs[Deriv::dT] *= minusDens2;
  for( integer ic = 0; ic < numComp; ++ic )
  {
    totalDensity.derivs[Deriv::dC+ic] *= minusDens2;
  }
}

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
