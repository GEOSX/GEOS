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
 * @file HydraulicApertureTable.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTURETABLE_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTURETABLE_HPP_

#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "functions/TableFunction.hpp"


namespace geos
{

namespace constitutive
{

/**
 * @class HydraulicApertureTableUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class HydraulicApertureTableUpdates
{
public:

  HydraulicApertureTableUpdates( TableFunction const & apertureTable )
    : m_apertureTable( apertureTable.createKernelWrapper() )
  {}

  /// Default copy constructor
  HydraulicApertureTableUpdates( HydraulicApertureTableUpdates const & ) = default;

  /// Default move constructor
  HydraulicApertureTableUpdates( HydraulicApertureTableUpdates && ) = default;

  /// Deleted default constructor
  HydraulicApertureTableUpdates() = default;

  /// Deleted copy assignment operator
  HydraulicApertureTableUpdates & operator=( HydraulicApertureTableUpdates const & ) = delete;

  /// Deleted move assignment operator
  HydraulicApertureTableUpdates & operator=( HydraulicApertureTableUpdates && ) =  delete;

  /**
   * @brief Evaluate the effective aperture, and its derivative wrt aperture
   * @param[in] aperture the model aperture/gap
   * @param[out] dHydraulicAperture_dAperture the derivative of the effective aperture wrt aperture
   * @return The hydraulic aperture that is always > 0
   */
  GEOS_HOST_DEVICE
  real64 computeHydraulicAperture( real64 const aperture,
                                   real64 const normalTraction,
                                   real64 & dHydraulicAperture_daperture,
                                   real64 & dHydraulicAperture_dNormalStress ) const;

protected:

  /// The aperture table function wrapper
  TableFunction::KernelWrapper m_apertureTable;
};


/**
 * @class HydraulicApertureTable
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class HydraulicApertureTable : public HydraulicApertureBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  HydraulicApertureTable( string const & name,
                          Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~HydraulicApertureTable() override;

  static string catalogName() { return "HydraulicApertureTable"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override final;



  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = HydraulicApertureTableUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @struct Structure to hold scoped key names
   */
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    /// string/key for aperture tolerance
    static constexpr char const * apertureToleranceString() { return "apertureTolerance"; }

    /// string/key for aperture table name
    static constexpr char const * apertureTableNameString() { return "apertureTableName"; }
  };


protected:

  virtual void postInputInitialization() override;

  /**
   * @brief Validate the values provided in the aperture table
   * @param[in] apertureTable the effective aperture vs aperture table
   */
  void validateApertureTable( TableFunction const & apertureTable ) const;

  /// The aperture tolerance to avoid floating point errors in expressions involving aperture
  real64 m_apertureTolerance;

  /// The name of the aperture table, if any
  string m_apertureTableName;

  /// Pointer to the function that limits the model aperture to a physically admissible value.
  TableFunction const * m_apertureTable;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 HydraulicApertureTableUpdates::computeHydraulicAperture( real64 const aperture,
                                                                real64 const normalTraction,
                                                                real64 & dHydraulicAperture_dAperture,
                                                                real64 & dHydraulicAperture_dNormalStress ) const
{
  GEOS_UNUSED_VAR( normalTraction, dHydraulicAperture_dNormalStress );
  return m_apertureTable.compute( &aperture, &dHydraulicAperture_dAperture );
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTURETABLE_HPP_ */
