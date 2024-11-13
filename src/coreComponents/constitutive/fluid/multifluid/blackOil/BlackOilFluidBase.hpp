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
 * @file BlackOilFluidBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_BLACKOILFLUIDBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_BLACKOILFLUIDBASE_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "functions/TableFunction.hpp"


namespace geos
{

namespace constitutive
{

class BlackOilFluidBase : public MultiFluidBase
{
public:

  using exec_policy = parallelDevicePolicy<>;

  static constexpr integer MAX_NUM_PHASES = 3;

  struct PhaseType
  {
    enum : integer
    {
      OIL = 0,
      GAS = 1,
      WATER = 2,
    };
  };

  BlackOilFluidBase( string const & name, Group * const parent );

  /**
   * @copydoc MultiFluidBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  virtual void checkTablesParameters( real64 pressure, real64 temperature ) const override final;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * surfacePhaseMassDensitiesString() { return "surfaceDensities"; }
    static constexpr char const * tableFilesString() { return "tableFiles"; }
    static constexpr char const * formationVolumeFactorTableNamesString() { return "hydrocarbonFormationVolFactorTableNames"; }
    static constexpr char const * viscosityTableNamesString() { return "hydrocarbonViscosityTableNames"; }
    static constexpr char const * waterRefPressureString() { return "waterReferencePressure"; }
    static constexpr char const * waterFormationVolumeFactorString() { return "waterFormationVolumeFactor"; }
    static constexpr char const * waterCompressibilityString() { return "waterCompressibility"; }
    static constexpr char const * waterViscosityString() { return "waterViscosity"; }
  };

protected:

  struct WaterParams
  {
    real64 referencePressure = 0.0;  ///< Water reference pressure
    real64 formationVolFactor = 0.0; ///< Water formation volume factor
    real64 compressibility = 0.0;    ///< Water compressibility
    real64 viscosity = 0.0;          ///< Water viscosity
  };

  class KernelWrapper : public MultiFluidBase::KernelWrapper
  {
public:

    /// @cond DO_NOT_DOCUMENT
    /// We need these SMFs to avoid host-device errors with CUDA.
    KernelWrapper() = default;
    KernelWrapper( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper && ) = default;
    /// @endcond

protected:

    KernelWrapper( arrayView1d< integer const > phaseTypes,
                   arrayView1d< integer const > phaseOrder,
                   arrayView1d< integer const > hydrocarbonPhaseOrder,
                   arrayView1d< real64 const > surfacePhaseMassDensity,
                   arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                   arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                   WaterParams const waterParams,
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

    /// Phase ordering info
    arrayView1d< integer const > m_phaseTypes;
    arrayView1d< integer const > m_phaseOrder;
    arrayView1d< integer const > m_hydrocarbonPhaseOrder;

    /// Surface mass density for each phase
    arrayView1d< real64 const > m_surfacePhaseMassDensity;

    /// Table kernel wrappers to interpolate in the oil and gas (B vs p) tables
    arrayView1d< TableFunction::KernelWrapper const > m_formationVolFactorTables;

    /// Table kernel wrappers to interpolate in the oil and gas (\mu vs p) tables
    arrayView1d< TableFunction::KernelWrapper const > m_viscosityTables;

    /// Water parameters
    WaterParams m_waterParams;
  };

  virtual integer getWaterPhaseIndex() const override final;

  virtual void postInputInitialization() override;

  virtual void initializePostSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers
   */
  virtual void createAllKernelWrappers();

  /**
   * @brief Use the TableFunctions provided by the user to get the PVT data
   */
  virtual void readInputDataFromTableFunctions() = 0;

  /**
   * @brief Read all the PVT table provided by the user in Eclipse format
   */
  virtual void readInputDataFromPVTFiles() = 0;

  /**
   * @brief Fill the water data (formation vol factor, compressibility, etc)
   * @param[in] tableValues the values in the water table
   */
  void fillWaterData( array1d< array1d< real64 > > const & tableValues );

  /**
   * @brief Fill the hydrocarbon data (pressure, formation vol factor, viscosity)
   * @param[in] ip the index of the phase
   * @param[in] tableValues the values in the oil or gas table
   */
  void fillHydrocarbonData( integer const ip,
                            array1d< array1d< real64 > > const & tableValues );

  /**
   * @brief Check that the table values make sense
   * @param[in] table the values in the oil or gas table
   * @param[in] warningIfDecreasing flag to issue a warning if values are decreasing (otherwise, issue a warning if increasing)
   * @detail This function throws an error if the table is invalid
   */
  virtual void validateTable( TableFunction const & table,
                              bool warningIfDecreasing ) const;


  /**
   * @brief Check water parameters for correctness.
   */
  void validateWaterParams() const;

  // Input data

  // Black-oil table filenames
  path_array m_tableFiles;

  // Fluid data

  /// Names of the formation volume factor tables
  array1d< string > m_formationVolFactorTableNames;

  /// Names of the viscosity tables
  array1d< string > m_viscosityTableNames;

  /// Surface densities
  array1d< real64 > m_surfacePhaseMassDensity;

  /// Water parameters
  WaterParams m_waterParams;

  /// Data after processing of input

  // Phase ordering info (all the constitutive arrays (density, viscosity, fractions, etc) use the phase ordering of the XML

  /// Map from the phase ordering defined in the XML (by the user) to the phase ordering in struct PhaseType
  array1d< integer > m_phaseTypes;
  /// Map from the phase ordering in struct PhaseType to the phase ordering defined in the XML (by the user)
  array1d< integer > m_phaseOrder;
  /// Map from the phase ordering in struct PhaseType to the hydrocarbon phase ordering deduced from the XML
  array1d< integer > m_hydrocarbonPhaseOrder;

  /// Table kernel wrappers to interpolate in the oil and gas (B vs p) tables
  array1d< TableFunction const * > m_formationVolFactorTables;

  /// Table kernel wrappers of m_formationVolFactorTables
  array1d< TableFunction::KernelWrapper > m_formationVolFactorTableKernels;

  /// Table kernel wrappers to interpolate in the oil and gas (\mu vs p) tables
  array1d< TableFunction const * > m_viscosityTables;

  /// Table kernel wrappers of m_viscosityTables
  array1d< TableFunction::KernelWrapper > m_viscosityTableKernels;

};

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_BLACKOILFLUIDBASE_HPP_
