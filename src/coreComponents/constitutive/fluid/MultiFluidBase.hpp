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
 * @file MultiFluidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"

namespace geosx
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
  static constexpr integer MAX_NUM_COMPONENTS = 16;

  /**
   * @brief Maximum supported number of fluid phases
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_PHASES = 4;

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

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseFraction_dPressure() const
  { return m_phaseFraction.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseFraction_dTemperature() const
  { return m_phaseFraction.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseFraction_dGlobalCompFraction() const
  { return m_phaseFraction.dComp; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity() const
  { return m_phaseDensity.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseDensity_dPressure() const
  { return m_phaseDensity.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseDensity_dTemperature() const
  { return m_phaseDensity.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseDensity_dGlobalCompFraction() const
  { return m_phaseDensity.dComp; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseMassDensity() const
  { return m_phaseMassDensity.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseMassDensity_dPressure() const
  { return m_phaseMassDensity.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseMassDensity_dTemperature() const
  { return m_phaseMassDensity.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseMassDensity_dGlobalCompFraction() const
  { return m_phaseMassDensity.dComp; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseViscosity() const
  { return m_phaseViscosity.value; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseViscosity_dPressure() const
  { return m_phaseViscosity.dPres; }

  arrayView3d< real64 const, multifluid::USD_PHASE > dPhaseViscosity_dTemperature() const
  { return m_phaseViscosity.dTemp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseViscosity_dGlobalCompFraction() const
  { return m_phaseViscosity.dComp; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > phaseCompFraction() const
  { return m_phaseCompFraction.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > dPhaseCompFraction_dPressure() const
  { return m_phaseCompFraction.dPres; }

  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > dPhaseCompFraction_dTemperature() const
  { return m_phaseCompFraction.dTemp; }

  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > dPhaseCompFraction_dGlobalCompFraction() const
  { return m_phaseCompFraction.dComp; }

  arrayView2d< real64 const, multifluid::USD_FLUID > totalDensity() const
  { return m_totalDensity.value; }

  arrayView2d< real64 const, multifluid::USD_FLUID > dTotalDensity_dPressure() const
  { return m_totalDensity.dPres; }

  arrayView2d< real64 const, multifluid::USD_FLUID > dTotalDensity_dTemperature() const
  { return m_totalDensity.dTemp; }

  arrayView3d< real64 const, multifluid::USD_FLUID_DC > dTotalDensity_dGlobalCompFraction() const
  { return m_totalDensity.dComp; }

  arrayView2d< real64 const, multifluid::USD_FLUID > initialTotalMassDensity() const
  { return m_initialTotalMassDensity.toViewConst(); }


  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * componentNamesString() { return "componentNames"; }
    static constexpr char const * componentMolarWeightString() { return "componentMolarWeight"; }
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * useMassString() { return "useMass"; }
  };

protected:

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
  using FluidProp = MultiFluidVar< real64, 2, multifluid::LAYOUT_FLUID, multifluid::LAYOUT_FLUID_DC >;

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
    GEOSX_HOST_DEVICE
    localIndex numElems() const { return m_phaseFraction.value.size( 0 ); }

    /**
     * @brief Get number of gauss points per element.
     * @return number of gauss points per element
     */
    GEOSX_HOST_DEVICE
    localIndex numGauss() const { return m_phaseFraction.value.size( 1 ); }

    /**
     * @brief Get number of fluid components.
     * @return number of components
     */
    GEOSX_HOST_DEVICE
    integer numComponents() const { return LvArray::integerConversion< integer >( m_componentMolarWeight.size() ); }

    /**
     * @brief Get number of fluid phases.
     * @return number of phases
     */
    GEOSX_HOST_DEVICE
    integer numPhases() const { return LvArray::integerConversion< integer >( m_phaseFraction.value.size( 2 ) ); }

protected:

    /**
     * @brief Constructor for the kernel wrapper
     * @param[in] componentMolarWeight the array of component molar weight
     * @param[in] useMass the flag to decide whether the solver works in units of mass or moles
     * @param[out] phaseFraction the array of phase fractions (+ derivatives)
     * @param[out] phaseDensity the array of phase densities (+ derivatives)
     * @param[out] phaseMassDensity the array of phase mass densities (+derivatives)
     * @param[out] phaseViscosity the array of phase viscosities (+derivatives)
     * @param[out] phaseCompFraction the array of phase component fractions (+derivatives)
     * @param[out] totalDensity the total density (+derivatives)
     */
    KernelWrapper( arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity )
      : m_componentMolarWeight( std::move( componentMolarWeight ) ),
      m_useMass( useMass ),
      m_phaseFraction( std::move( phaseFraction ) ),
      m_phaseDensity( std::move( phaseDensity ) ),
      m_phaseMassDensity( std::move( phaseMassDensity ) ),
      m_phaseViscosity( std::move( phaseViscosity ) ),
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
    GEOSX_HOST_DEVICE
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
    GEOSX_HOST_DEVICE
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
    GEOSX_HOST_DEVICE
    void convertToMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                 arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                                 arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const phaseCompFrac ) const;

    /**
     * @brief Utility function to convert mole fractions to mass fractions and keep derivatives
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] dCompMoleFrac_dCompMassFrac the derivatives of mole fractions wrt mass fractions
     * @param[in] phaseMolecularWeight the phase molecular weight computed by the constitutive model
     * @param[in] dPhaseMolecularWeight_dPres the derivatives of phase molecular weights wrt pressure
     * @param[in] dPhaseMolecularWeight_dTemp the derivatives of phase molecular weights wrt temperature
     * @param[in] dPhaseMolecularWeight_dGlobalCompFrac the derivatives of phase molecular weights wrt comp fractions
     * @param[inout] phaseFrac the phase fractions in moles that will be converted to mass
     * @param[inout] phaseCompFrac the phase component fractions in moles that will be converted to mass
     * @param[inout] dPhaseDens_dGlobalCompFrac the derivatives of phase densities wrt comp fractions
     * @param[inout] dPhaseVisc_dGlobalCompFrac the derivatives of phase viscosities wrt comp fractions
     * @detail This function performs three conversions
     *    1) Conversion of phase mass fractions into phase mole fractions
     *    2) Conversion of phase component mass fractions into phase component mole fractions
     *    3) Conversion of derivatives wrt mass fractions into derivatives wrt mole fractions
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOSX_HOST_DEVICE
    void convertToMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                                 real64 const (&phaseMolecularWeight)[maxNumPhase],
                                 real64 const (&dPhaseMolecularWeight_dPres)[maxNumPhase],
                                 real64 const (&dPhaseMolecularWeight_dTemp)[maxNumPhase],
                                 real64 const (&dPhaseMolecularWeight_dGlobalCompFrac)[maxNumPhase][maxNumComp],
                                 PhaseProp::SliceType const phaseFrac,
                                 PhaseComp::SliceType const phaseCompFrac,
                                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens_dGlobalCompFrac,
                                 arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc_dGlobalCompFrac ) const;

    /**
     * @brief Utility function to convert mole fractions to mass fractions
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] phaseFrac the phase fractions properly converted
     * @param[in] phaseFrac the phase densities in mass or moles
     * @param[out] totalDens the total density
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOSX_HOST_DEVICE
    void computeTotalDensity( arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                              arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseDens,
                              real64 & totalDens ) const;

    /**
     * @brief Utility function to convert mole fractions to mass fractions and keep derivatives
     * @param[in] phaseFrac the phase fractions properly converted (+ derivatives)
     * @param[in] phaseFrac the phase densities in mass or moles (+ derivatives)
     * @param[out] totalDens the total density (+ derivatives)
     */
    GEOSX_HOST_DEVICE
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
    PhaseComp::ViewType m_phaseCompFraction;
    FluidProp::ViewType m_totalDensity;

private:

    /**
     * @brief Utility function to convert phase mole fractions to phase mass fractions and keep derivatives
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] phaseMolecularWeight the phase molecular weight computed by the constitutive model
     * @param[in] dPhaseMolecularWeight_dPres the derivatives of phase molecular weights wrt pressure
     * @param[in] dPhaseMolecularWeight_dTemp the derivatives of phase molecular weights wrt temperature
     * @param[in] dPhaseMolecularWeight_dGlobalCompFrac the derivatives of phase molecular weights wrt comp fractions
     * @param[inout] phaseFrac the phase fractions in moles that will be converted to mass
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOSX_HOST_DEVICE
    void convertToPhaseMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                      real64 const (&dPhaseMolecularWeight_dPres)[maxNumPhase],
                                      real64 const (&dPhaseMolecularWeight_dTemp)[maxNumPhase],
                                      real64 const (&dPhaseMolecularWeight_dGlobalCompFrac)[maxNumPhase][maxNumComp],
                                      PhaseProp::SliceType const phaseFrac ) const;

    /**
     * @brief Utility function to convert phase component mole fractions to phase component mass fractions and keep derivatives
     * @tparam maxNumComp the max number of components
     * @tparam maxNumPhase the max number of phases
     * @param[in] phaseMolecularWeight the phase molecular weight computed by the constitutive model
     * @param[in] dPhaseMolecularWeight_dPres the derivatives of phase molecular weights wrt pressure
     * @param[in] dPhaseMolecularWeight_dTemp the derivatives of phase molecular weights wrt temperature
     * @param[in] dPhaseMolecularWeight_dGlobalCompFrac the derivatives of phase molecular weights wrt comp fractions
     * @param[inout] phaseFrac the phase fractions
     * @param[inout] phaseCompFrac the phase component fractions in moles that will be converted to mass
     */
    template< integer maxNumComp, integer maxNumPhase >
    GEOSX_HOST_DEVICE
    void convertToPhaseCompMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                          real64 const (&dPhaseMolecularWeight_dPres)[maxNumPhase],
                                          real64 const (&dPhaseMolecularWeight_dTemp)[maxNumPhase],
                                          real64 const (&dPhaseMolecularWeight_dGlobalCompFrac)[maxNumPhase][maxNumComp],
                                          PhaseProp::SliceType const phaseFrac,
                                          PhaseComp::SliceType const phaseCompFrac ) const;

    /**
     * @brief Utility function to convert derivatives wrt mole fractions into derivatives wrt mass fractions
     * @tparam maxNumComp the max number of components
     * @param[in] dCompMoleFrac_dCompMassFrac the derivatives of mole fractions wrt mass fractions
     * @param[inout] phaseFrac the phase fractions
     * @param[inout] phaseCompFrac the phase component fractions
     * @param[inout] dPhaseDens_dGlobalCompFrac the derivatives of phase densities wrt comp fractions
     * @param[inout] dPhaseVisc_dGlobalCompFrac the derivatives of phase viscosities wrt comp fractions
     */
    template< integer maxNumComp >
    GEOSX_HOST_DEVICE
    void computeDerivativesWrtMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                                             PhaseProp::SliceType const phaseFrac,
                                             PhaseComp::SliceType const phaseCompFrac,
                                             arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens_dGlobalCompFrac,
                                             arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc_dGlobalCompFrac ) const;


    /**
     * @brief Main compute function to update properties in a cell without returning derivatives (used at initialization)
     * @param[in] pressure pressure in the cell
     * @param[in] temperature temperature in the cell
     * @param[in] composition mass/molar component fractions in the cell
     * @param[out] phaseFraction phase fractions in the cell
     * @param[out] phaseDensity phase mass/molar density in the cell
     * @param[out] phaseMassDensity phase mass density in the cell
     * @param[out] phaseViscosity phase viscosity in the cell
     * @param[out] phaseCompFraction phase component fraction in the cell
     * @param[out] totalDensity total mass/molar density in the cell
     */
    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          real64 & totalDensity ) const = 0;

    /**
     * @brief Main compute function to update properties in a cell with derivatives (used in Newton iterations)
     * @param[in] pressure pressure in the cell
     * @param[in] temperature temperature in the cell
     * @param[in] composition mass/molar component fractions in the cell
     * @param[out] phaseFraction phase fractions in the cell  (+ derivatives)
     * @param[out] phaseDensity phase mass/molar density in the cell (+ derivatives)
     * @param[out] phaseMassDensity phase mass density in the cell (+ derivatives)
     * @param[out] phaseViscosity phase viscosity in the cell (+ derivatives)
     * @param[out] phaseCompFraction phase component fraction in the cell (+ derivatives)
     * @param[out] totalDensity total mass/molar density in the cell (+ derivatives)
     */
    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseProp::SliceType const phaseDensity,
                          PhaseProp::SliceType const phaseMassDensity,
                          PhaseProp::SliceType const phaseViscosity,
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
    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const = 0;
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

  virtual void postProcessInput() override;

  // flag indicating whether input/output component fractions are treated as mass fractions
  int m_useMass;

  // general fluid composition information

  array1d< string > m_componentNames;
  array1d< real64 > m_componentMolarWeight;
  array1d< string > m_phaseNames;

  // constitutive data

  PhaseProp m_phaseFraction;
  PhaseProp m_phaseDensity;
  PhaseProp m_phaseMassDensity;
  PhaseProp m_phaseViscosity;
  PhaseComp m_phaseCompFraction;
  FluidProp m_totalDensity;

  // initial data (used to compute the body force in the poromechanics solver)

  array2d< real64, multifluid::LAYOUT_FLUID > m_initialTotalMassDensity;

};

template< integer maxNumComp, typename OUT_ARRAY >
GEOSX_HOST_DEVICE
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
GEOSX_HOST_DEVICE
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
GEOSX_HOST_DEVICE
void
MultiFluidBase::KernelWrapper::
  convertToMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const phaseCompFrac ) const
{
  using namespace multifluid;

  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  real64 dCompMoleFrac_dCompMassFrac[maxNumComp][maxNumComp]{};
  real64 dPhaseMolecularWeight_dPres[maxNumPhase]{};
  real64 dPhaseMolecularWeight_dTemp[maxNumPhase]{};
  real64 dPhaseMolecularWeight_dComp[maxNumPhase][maxNumComp]{};

  StackArray< real64, 3, maxNumPhase, LAYOUT_PHASE > dPhaseFrac_dPres( 1, 1, numPhase );
  StackArray< real64, 3, maxNumPhase, LAYOUT_PHASE > dPhaseFrac_dTemp( 1, 1, numPhase );
  StackArray< real64, 4, maxNumComp *maxNumPhase, LAYOUT_PHASE_DC > dPhaseFrac_dComp( 1, 1, numPhase, numComp );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac, dPhaseFrac_dPres[0][0], dPhaseFrac_dTemp[0][0], dPhaseFrac_dComp[0][0] };

  StackArray< real64, 4, maxNumComp *maxNumPhase, LAYOUT_PHASE_COMP > dPhaseCompFrac_dPres( 1, 1, numPhase, numComp );
  StackArray< real64, 4, maxNumComp *maxNumPhase, LAYOUT_PHASE_COMP > dPhaseCompFrac_dTemp( 1, 1, numPhase, numComp );
  StackArray< real64, 5, maxNumComp *maxNumComp *maxNumPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFrac_dComp( 1, 1, numPhase, numComp, numComp );
  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 >
  phaseCompFracAndDeriv { phaseCompFrac, dPhaseCompFrac_dPres[0][0], dPhaseCompFrac_dTemp[0][0], dPhaseCompFrac_dComp[0][0] };

  StackArray< real64, 4, maxNumComp *maxNumPhase, LAYOUT_PHASE_DC > dPhaseDens_dComp( 1, 1, numPhase, numComp );
  StackArray< real64, 4, maxNumComp *maxNumPhase, LAYOUT_PHASE_DC > dPhaseVisc_dComp( 1, 1, numPhase, numComp );

  convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                          phaseMolecularWeight,
                          dPhaseMolecularWeight_dPres,
                          dPhaseMolecularWeight_dTemp,
                          dPhaseMolecularWeight_dComp,
                          phaseFracAndDeriv,
                          phaseCompFracAndDeriv,
                          dPhaseDens_dComp[0][0],
                          dPhaseVisc_dComp[0][0] );
}

template< integer maxNumComp, integer maxNumPhase >
GEOSX_HOST_DEVICE
void
MultiFluidBase::KernelWrapper::
  convertToMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                          real64 const (&phaseMolecularWeight)[maxNumPhase],
                          real64 const (&dPhaseMolecularWeight_dPres)[maxNumPhase],
                          real64 const (&dPhaseMolecularWeight_dTemp)[maxNumPhase],
                          real64 const (&dPhaseMolecularWeight_dGlobalCompFrac)[maxNumPhase][maxNumComp],
                          PhaseProp::SliceType const phaseFrac,
                          PhaseComp::SliceType const phaseCompFrac,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens_dGlobalCompFrac,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc_dGlobalCompFrac ) const
{
  convertToPhaseMassFractions( phaseMolecularWeight,
                               dPhaseMolecularWeight_dPres,
                               dPhaseMolecularWeight_dTemp,
                               dPhaseMolecularWeight_dGlobalCompFrac,
                               phaseFrac );

  convertToPhaseCompMassFractions( phaseMolecularWeight,
                                   dPhaseMolecularWeight_dPres,
                                   dPhaseMolecularWeight_dTemp,
                                   dPhaseMolecularWeight_dGlobalCompFrac,
                                   phaseFrac,
                                   phaseCompFrac );

  computeDerivativesWrtMassFractions( dCompMoleFrac_dCompMassFrac,
                                      phaseFrac,
                                      phaseCompFrac,
                                      dPhaseDens_dGlobalCompFrac,
                                      dPhaseVisc_dGlobalCompFrac );
}

template< integer maxNumComp, integer maxNumPhase >
GEOSX_HOST_DEVICE
void
MultiFluidBase::KernelWrapper::
  convertToPhaseMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                               real64 const (&dPhaseMolecularWeight_dPres)[maxNumPhase],
                               real64 const (&dPhaseMolecularWeight_dTemp)[maxNumPhase],
                               real64 const (&dPhaseMolecularWeight_dGlobalCompFrac)[maxNumPhase][maxNumComp],
                               PhaseProp::SliceType const phaseFrac ) const
{
  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  real64 totalMass{};
  real64 dTotalMass_dP{};
  real64 dTotalMass_dT{};
  real64 dTotalMass_dC[maxNumComp]{};

  // 1. Compute mass of each phase and total mass (on a 1-mole basis)
  for( integer ip = 0; ip < numPhase; ++ip )
  {

    real64 const nu = phaseFrac.value[ip];

    phaseFrac.value[ip] *= phaseMolecularWeight[ip];
    phaseFrac.dPres[ip] = phaseFrac.dPres[ip] * phaseMolecularWeight[ip] + nu * dPhaseMolecularWeight_dPres[ip];
    phaseFrac.dTemp[ip] = phaseFrac.dTemp[ip] * phaseMolecularWeight[ip] + nu * dPhaseMolecularWeight_dTemp[ip];

    totalMass += phaseFrac.value[ip];
    dTotalMass_dP += phaseFrac.dPres[ip];
    dTotalMass_dT += phaseFrac.dTemp[ip];

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFrac.dComp[ip][jc] =
        phaseFrac.dComp[ip][jc] * phaseMolecularWeight[ip] + nu * dPhaseMolecularWeight_dGlobalCompFrac[ip][jc];
      dTotalMass_dC[jc] += phaseFrac.dComp[ip][jc];
    }
  }

  // 2. Normalize to get mass fractions
  real64 const totalMassInv = 1.0 / totalMass;
  for( integer ip = 0; ip < numPhase; ++ip )
  {
    phaseFrac.value[ip] *= totalMassInv;
    phaseFrac.dPres[ip] = ( phaseFrac.dPres[ip] - phaseFrac.value[ip] * dTotalMass_dP ) * totalMassInv;
    phaseFrac.dTemp[ip] = ( phaseFrac.dTemp[ip] - phaseFrac.value[ip] * dTotalMass_dT ) * totalMassInv;

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFrac.dComp[ip][jc] = ( phaseFrac.dComp[ip][jc] - phaseFrac.value[ip] * dTotalMass_dC[jc] ) * totalMassInv;
    }
  }
}

template< integer maxNumComp, integer maxNumPhase >
GEOSX_HOST_DEVICE
void
MultiFluidBase::KernelWrapper::
  convertToPhaseCompMassFractions( real64 const (&phaseMolecularWeight)[maxNumPhase],
                                   real64 const (&dPhaseMolecularWeight_dPres)[maxNumPhase],
                                   real64 const (&dPhaseMolecularWeight_dTemp)[maxNumPhase],
                                   real64 const (&dPhaseMolecularWeight_dGlobalCompFrac)[maxNumPhase][maxNumComp],
                                   PhaseProp::SliceType const phaseFrac,
                                   PhaseComp::SliceType const phaseCompFrac ) const
{
  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    real64 const phaseMolecularWeightInv = 1.0 / phaseMolecularWeight[ip];

    for( integer ic = 0; ic < numComp; ++ic )
    {
      real64 const compMW = m_componentMolarWeight[ic];

      phaseCompFrac.value[ip][ic] = phaseCompFrac.value[ip][ic] * compMW * phaseMolecularWeightInv;
      phaseCompFrac.dPres[ip][ic] =
        ( phaseCompFrac.dPres[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMolecularWeight_dPres[ip] ) * phaseMolecularWeightInv;
      phaseCompFrac.dTemp[ip][ic] =
        ( phaseCompFrac.dTemp[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMolecularWeight_dTemp[ip] ) * phaseMolecularWeightInv;

      for( integer jc = 0; jc < numComp; ++jc )
      {
        phaseCompFrac.dComp[ip][ic][jc] =
          phaseMolecularWeightInv *
          ( phaseCompFrac.dComp[ip][ic][jc] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMolecularWeight_dGlobalCompFrac[ip][jc] );
      }
    }
  }
}

template< integer maxNumComp >
GEOSX_HOST_DEVICE
void
MultiFluidBase::KernelWrapper::
  computeDerivativesWrtMassFractions( real64 const (&dCompMoleFrac_dCompMassFrac)[maxNumComp][maxNumComp],
                                      PhaseProp::SliceType const phaseFrac,
                                      PhaseComp::SliceType const phaseCompFrac,
                                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseDens_dGlobalCompFrac,
                                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc_dGlobalCompFrac ) const
{
  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  real64 work[maxNumComp]{};
  for( integer ip = 0; ip < numPhase; ++ip )
  {
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, phaseFrac.dComp[ip], work );
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, dPhaseDens_dGlobalCompFrac[ip], work );
    applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, dPhaseVisc_dGlobalCompFrac[ip], work );
    for( integer ic = 0; ic < numComp; ++ic )
    {
      applyChainRuleInPlace( numComp, dCompMoleFrac_dCompMassFrac, phaseCompFrac.dComp[ip][ic], work );
    }
  }
}

template< integer maxNumComp, integer maxNumPhase >
GEOSX_HOST_DEVICE
inline void
MultiFluidBase::KernelWrapper::
  computeTotalDensity( arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseFrac,
                       arraySlice1d< real64, multifluid::USD_PHASE - 2 > const phaseDens,
                       real64 & totalDens ) const
{
  using namespace multifluid;

  integer const numPhase = numPhases();
  integer const numComp = numComponents();

  StackArray< real64, 2, 1, LAYOUT_FLUID > dTotalDens_dPres( 1, 1 );
  StackArray< real64, 2, 1, LAYOUT_FLUID > dTotalDens_dTemp( 1, 1 );
  StackArray< real64, 3, maxNumComp, LAYOUT_FLUID_DC > dTotalDens_dComp( 1, 1, numComp );
  MultiFluidVarSlice< real64, 0, USD_FLUID - 2, USD_FLUID_DC - 2 >
  totalDensAndDeriv { totalDens, dTotalDens_dPres[0][0], dTotalDens_dTemp[0][0], dTotalDens_dComp[0][0] };

  StackArray< real64, 3, maxNumPhase, LAYOUT_PHASE > dPhaseFrac_dPres( 1, 1, numPhase );
  StackArray< real64, 3, maxNumPhase, LAYOUT_PHASE > dPhaseFrac_dTemp( 1, 1, numPhase );
  StackArray< real64, 4, maxNumComp *maxNumPhase, LAYOUT_PHASE_DC > dPhaseFrac_dComp( 1, 1, numPhase, numComp );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac, dPhaseFrac_dPres[0][0], dPhaseFrac_dTemp[0][0], dPhaseFrac_dComp[0][0] };

  StackArray< real64, 3, maxNumPhase, LAYOUT_PHASE > dPhaseDens_dPres( 1, 1, numPhase );
  StackArray< real64, 3, maxNumPhase, LAYOUT_PHASE > dPhaseDens_dTemp( 1, 1, numPhase );
  StackArray< real64, 4, maxNumComp *maxNumPhase, LAYOUT_PHASE_DC > dPhaseDens_dComp( 1, 1, numPhase, numComp );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseDensAndDeriv { phaseDens, dPhaseDens_dPres[0][0], dPhaseDens_dTemp[0][0], dPhaseDens_dComp[0][0] };

  computeTotalDensity( phaseFracAndDeriv,
                       phaseDensAndDeriv,
                       totalDensAndDeriv );
}

GEOSX_HOST_DEVICE
inline void
MultiFluidBase::KernelWrapper::
  computeTotalDensity( PhaseProp::SliceType const phaseFraction,
                       PhaseProp::SliceType const phaseDensity,
                       FluidProp::SliceType const totalDensity ) const
{
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  totalDensity.value = 0.0;
  totalDensity.dPres = 0.0;
  LvArray::forValuesInSlice( totalDensity.dComp, []( real64 & val ){ val = 0.0; } );

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
    totalDensity.dPres += ( phaseFraction.dPres[ip] - value * phaseDensity.dPres[ip] ) * densInv;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity.dComp[ic] += ( phaseFraction.dComp[ip][ic] - value * phaseDensity.dComp[ip][ic] ) * densInv;
    }
  }

  // 2. Invert the previous quantity to get actual density
  totalDensity.value = 1.0 / totalDensity.value;
  real64 const minusDens2 = -totalDensity.value * totalDensity.value;
  totalDensity.dPres *= minusDens2;
  for( integer ic = 0; ic < numComp; ++ic )
  {
    totalDensity.dComp[ic] *= minusDens2;
  }
}


} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
