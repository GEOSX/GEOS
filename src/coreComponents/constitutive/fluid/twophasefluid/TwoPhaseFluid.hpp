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
 * @file TwoPhaseFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_TWOPHASEFLUID_TWOPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_TWOPHASEFLUID_TWOPHASEFLUID_HPP_

#include "common/DataLayouts.hpp"
#include "functions/TableFunction.hpp"
#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

#include "constitutive/ConstitutivePassThruHandler.hpp"


namespace geos
{
namespace constitutive
{

class TwoPhaseFluid : public ConstitutiveBase
{
public:

  TwoPhaseFluid( string const & name,
                 Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return "TwoPhaseFluid"; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  ///@}



  /**
   * @brief Getter for the fluid phase names
   * @return an array storing the phase names
   */
  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * tableFilesString() { return "tableFiles"; }
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * densityTableNamesString() { return "densityTableNames"; }
    static constexpr char const * viscosityTableNamesString() { return "viscosityTableNames"; }
  };

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity_n() const
  { return m_phaseDensity_n; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDensity() const
  { return m_phaseDensity.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseDensity() const
  { return m_phaseDensity.derivs; }

  arrayView3d< real64 const, multifluid::USD_PHASE > phaseViscosity() const
  { return m_phaseViscosity.value; }

  arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseViscosity() const
  { return m_phaseViscosity.derivs; }

  using PhaseProp = MultiFluidVar< real64, 3, constitutive::multifluid::LAYOUT_PHASE, constitutive::multifluid::LAYOUT_PHASE_DC >;


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
    localIndex numElems() const { return m_phaseDensity.value.size( 0 ); }

    /**
     * @brief Get number of gauss points per element.
     * @return number of gauss points per element
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    localIndex numGauss() const { return m_phaseDensity.value.size( 1 ); }


    GEOS_HOST_DEVICE
    void compute( real64 const pressure,
                  PhaseProp::SliceType const phaseDensity,
                  PhaseProp::SliceType const phaseViscosity ) const;

    GEOS_HOST_DEVICE
    void update( localIndex const k,
                 localIndex const q,
                 real64 const pressure ) const;

private:

    friend class TwoPhaseFluid;

    /**
     * @brief Constructor for the class doing in-kernel two-phase fluid updates
     * @param[in] densityTables density tables
     * @param[in] viscosityTables viscosity tables
     * @param[in] phaseDensity phase densities (+ derivatives) in the cell
     * @param[in] phaseViscosity phase viscosities (+ derivatives) in the cell
     */
    KernelWrapper(
      arrayView1d< TableFunction::KernelWrapper const > densityTables,
      arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
      PhaseProp::ViewType phaseDensity,
      PhaseProp::ViewType phaseViscosity );


protected:

    KernelWrapper(
      arrayView1d< TableFunction::KernelWrapper const > densityTables,
      arrayView1d< TableFunction::KernelWrapper const > viscosityTables
      );

    /// Table kernel wrappers to interpolate in the two phase (\rho vs p) tables
    arrayView1d< TableFunction::KernelWrapper const > m_densityTables;

    /// Table kernel wrappers to interpolate in the two phase (\mu vs p) tables
    arrayView1d< TableFunction::KernelWrapper const > m_viscosityTables;

    /**
     * @brief Utility function to compute densities as a function of pressure (keeping derivatives)
     * @param[in] pressure pressure in the cell
     * @param[out] phaseDensity the phase density in the cell (+ derivatives)
     */
    GEOS_HOST_DEVICE
    void computeDensities( real64 const pressure,
                           PhaseProp::SliceType const & phaseDensity ) const;

    /**
     * @brief Utility function to compute viscosities as a function of pressure (keeping derivatives)
     * @param[in] pressure pressure in the cell
     * @param[out] phaseViscosity the phase viscosities in the cell (+ derivatives)
     */
    GEOS_HOST_DEVICE
    void computeViscosities( real64 const pressure,
                             PhaseProp::SliceType const & phaseViscosity ) const;

    /// Views on the phase properties
    PhaseProp::ViewType m_phaseDensity;
    PhaseProp::ViewType m_phaseViscosity;

  };  //class KernelWrapper


  array1d< string > m_phaseNames;


  path_array m_tableFiles;

  /// Names of the density tables (one per phase)
  array1d< string > m_densityTableNames;

  /// Names of the viscosity tables (one per phase)
  array1d< string > m_viscosityTableNames;

  PhaseProp m_phaseDensity;
  PhaseProp m_phaseViscosity;

  /// Backup data
  array3d< real64, multifluid::LAYOUT_PHASE > m_phaseDensity_n;


  virtual void resizeFields( localIndex const size, localIndex const numPts );

  virtual void postInputInitialization() override;

  virtual void initializePostSubGroups() override;

  /// Table kernel wrappers to interpolate (\rho vs p) tables
  array1d< TableFunction const * > m_densityTables;
  /// Table kernel wrappers of m_densityTables
  array1d< TableFunction::KernelWrapper > m_densityTableKernels;

  /// Table kernel wrappers to interpolate (\mu vs p) tables
  array1d< TableFunction const * > m_viscosityTables;
  /// Table kernel wrappers of m_viscosityTables
  array1d< TableFunction::KernelWrapper > m_viscosityTableKernels;

  KernelWrapper createKernelWrapper();


private:
  void readInputDataFromTableFunctions();
  void readInputDataFromFileTableFunctions();

  /**
   * @brief Fill the fluid data (pressure, density, viscosity)
   * @param[in] ip the index of the phase
   * @param[in] tableValues the values in the fluid table
   */
  void fillData( integer const ip,
                 array1d< array1d< real64 > > const & tableValues );

  /**
   * @brief Check the monotonicity of the PVT relationship
   */
  void checkTableConsistency() const;
};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void TwoPhaseFluid::KernelWrapper::
  computeDensities( real64 const pressure,
                    PhaseProp::SliceType const & phaseDensity ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  LvArray::forValuesInSlice( phaseDensity.derivs, []( real64 & val ) { val = 0.0; } );

  for( integer iph = 0; iph < 2; ++iph )
  {
    // interpolate in the table to get the phase density and derivatives
    real64 dPhaseDens_dPres = 0.0;

    phaseDensity.value[iph] = m_densityTables[iph].compute( &pressure, &dPhaseDens_dPres );
    phaseDensity.derivs[iph][Deriv::dP] = dPhaseDens_dPres;
  }
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void TwoPhaseFluid::KernelWrapper::
  computeViscosities( real64 const pressure,
                      PhaseProp::SliceType const & phaseViscosity ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  LvArray::forValuesInSlice( phaseViscosity.derivs, []( real64 & val ) { val = 0.0; } );

  for( integer iph = 0; iph < 2; ++iph )
  {
    // interpolate in the table to get the phase viscosity and derivatives
    real64 dPhaseVisc_dPres = 0.0;
    phaseViscosity.value[iph] = m_viscosityTables[iph].compute( &pressure, &dPhaseVisc_dPres );
    phaseViscosity.derivs[iph][Deriv::dP] = dPhaseVisc_dPres;
  }
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void TwoPhaseFluid::KernelWrapper::
  compute( real64 const pressure,
           PhaseProp::SliceType const phaseDensity,
           PhaseProp::SliceType const phaseViscosity ) const
{
  computeDensities( pressure,
                    phaseDensity );

  computeViscosities( pressure,
                      phaseViscosity );
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void TwoPhaseFluid::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure
          ) const
{
  compute( pressure,
           m_phaseDensity( k, q ),
           m_phaseViscosity( k, q ) );
}


template< typename LAMBDA >
void constitutiveUpdatePassThru( TwoPhaseFluid const & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< TwoPhaseFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}


template< typename LAMBDA >
void constitutiveUpdatePassThru( TwoPhaseFluid & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< TwoPhaseFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

}  // namespace constitutive
}  // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_TWOPHASEFLUID_TWOPHASEFLUID_HPP_
