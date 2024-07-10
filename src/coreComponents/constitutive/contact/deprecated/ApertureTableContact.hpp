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
 * @file ApertureTableContact.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_APERTURETABLECONTACT_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_APERTURETABLECONTACT_HPP_

#include "constitutive/contact/ContactBase.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class ApertureTableContactUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class ApertureTableContactUpdates : public ContactBaseUpdates
{
public:

  ApertureTableContactUpdates( real64 const & penaltyStiffness,
                               TableFunction const & apertureTable )
    : ContactBaseUpdates( penaltyStiffness ),
    m_apertureTable( apertureTable.createKernelWrapper() )
  {}

  /// Default copy constructor
  ApertureTableContactUpdates( ApertureTableContactUpdates const & ) = default;

  /// Default move constructor
  ApertureTableContactUpdates( ApertureTableContactUpdates && ) = default;

  /// Deleted default constructor
  ApertureTableContactUpdates() = delete;

  /// Deleted copy assignment operator
  ApertureTableContactUpdates & operator=( ApertureTableContactUpdates const & ) = delete;

  /// Deleted move assignment operator
  ApertureTableContactUpdates & operator=( ApertureTableContactUpdates && ) =  delete;

  GEOS_HOST_DEVICE
  inline
  virtual real64 computeEffectiveAperture( real64 const aperture,
                                           real64 & dEffectiveAperture_dAperture ) const;

  GEOS_HOST_DEVICE
  inline
  virtual void computeTraction( arraySlice1d< real64 const > const & dispJump,
                                arraySlice1d< real64 > const & tractionVector,
                                arraySlice2d< real64 > const & dTractionVector_dJump ) const;

  GEOS_HOST_DEVICE
  inline
  void addPressureToTraction( real64 const & pressure,
                              bool const isOpen,
                              arraySlice1d< real64 >const & tractionVector,
                              real64 & dTraction_dPressure ) const;

private:

  /// The aperture table function wrapper
  TableFunction::KernelWrapper m_apertureTable;

};


/**
 * @class ApertureTableContact
 */
class ApertureTableContact : public ContactBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  ApertureTableContact( string const & name,
                        Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ApertureTableContact() override;

  /**
   * @brief Name that is used to register this a type of "ApertureTableContact" in the object catalog
   * @return See description
   */
  static string catalogName() { return "Contact"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief accessor for aperture tolerance
   * @return the aperture tolerance
   */
  real64 apertureTolerance() const { return m_apertureTolerance; }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ApertureTableContactUpdates;

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

  virtual void initializePreSubGroups() override;

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
real64 ApertureTableContactUpdates::computeEffectiveAperture( real64 const aperture,
                                                              real64 & dEffectiveAperture_dAperture ) const
{
  return m_apertureTable.compute( &aperture, &dEffectiveAperture_dAperture );
}


} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_APERTURETABLECONTACT_HPP_ */
