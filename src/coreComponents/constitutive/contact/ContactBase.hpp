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
 * @file ContactBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_CONTACTBASE_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_CONTACTBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"


namespace geos
{

namespace constitutive
{

/**
 * @class ContactBaseUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class ContactBaseUpdates
{
public:

  ContactBaseUpdates( real64 const & penaltyStiffness,
                      real64 const & shearStiffness,
                      real64 const & displacementJumpThreshold,
                      TableFunction const & apertureTable )
    : m_penaltyStiffness( penaltyStiffness ),
    m_shearStiffness( shearStiffness ),
    m_displacementJumpThreshold( displacementJumpThreshold ),
    m_apertureTable( apertureTable.createKernelWrapper() )
  {}

  /// Default copy constructor
  ContactBaseUpdates( ContactBaseUpdates const & ) = default;

  /// Default move constructor
  ContactBaseUpdates( ContactBaseUpdates && ) = default;

  /// Deleted default constructor
  ContactBaseUpdates() = default;

  /// Deleted copy assignment operator
  ContactBaseUpdates & operator=( ContactBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  ContactBaseUpdates & operator=( ContactBaseUpdates && ) =  delete;

  /**
   * @brief Evaluate the effective aperture, and its derivative wrt aperture
   * @param[in] aperture the model aperture/gap
   * @param[out] dHydraulicAperture_dAperture the derivative of the effective aperture wrt aperture
   * @return The hydraulic aperture that is always > 0
   */
  GEOS_HOST_DEVICE
  virtual real64 computeHydraulicAperture( real64 const aperture,
                                           real64 & dHydraulicAperture_dAperture ) const;


  /**
   * @brief Evaluate the traction vector and its derivatives wrt to pressure and jump
   * @param[in] dispJump the displacement jump
   * @param[in] fractureState the fracture state
   * @param[out] tractionVector the traction vector
   * @param[out] dTractionVector_dJump the derivative of the traction vector wrt displacement jump
   */
  GEOS_HOST_DEVICE
  inline
  virtual void computeTraction( localIndex const k,
                                arraySlice1d< real64 const > const & oldDispJump,
                                arraySlice1d< real64 const > const & dispJump,
                                integer const & fractureState,
                                arraySlice1d< real64 > const & tractionVector,
                                arraySlice2d< real64 > const & dTractionVector_dJump ) const
  {GEOS_UNUSED_VAR( k, oldDispJump, dispJump, tractionVector, dTractionVector_dJump, fractureState );}

  /**
   * @brief Evaluate the traction vector and its derivatives wrt to pressure and jump
   * @param[in] dispJump the displacement jump
   * @param[in] tractionVector the traction vector
   * @param[out] fractureState the fracture state
   */
  GEOS_HOST_DEVICE
  inline
  virtual void updateFractureState( localIndex const k,
                                    arraySlice1d< real64 const > const & dispJump,
                                    arraySlice1d< real64 const > const & tractionVector,
                                    integer & fractureState ) const
  { GEOS_UNUSED_VAR( k, dispJump, tractionVector, fractureState ); }

  /**
   * @brief Update the traction with the pressure term
   * @param[in] pressure the pressure term
   * @param[in] isOpen a flag specifying whether the fracture is open or closed
   * @param[inout] traction the current tractionVector
   * @param[out] dTraction_dPressure the derivative of the fist component of traction wrt pressure
   * @return the updated traction
   */
  GEOS_HOST_DEVICE
  virtual void addPressureToTraction( real64 const & pressure,
                                      arraySlice1d< real64 >const & tractionVector,
                                      real64 & dTraction_dPressure ) const;

  /**
   * @brief Evaluate the limit tangential traction norm and return the derivative wrt normal traction
   * @param[in] normalTraction the normal traction
   * @param[out] dLimitTangentialTractionNorm_dTraction the derivative of the limit tangential traction norm wrt normal traction
   * @return the limit tangential traction norm
   */
  GEOS_HOST_DEVICE
  inline
  virtual real64 computeLimitTangentialTractionNorm( real64 const & normalTraction,
                                                     real64 & dLimitTangentialTractionNorm_dTraction ) const
  { GEOS_UNUSED_VAR( normalTraction, dLimitTangentialTractionNorm_dTraction ); return 0; };

protected:

  /// The penalty stiffness
  real64 m_penaltyStiffness;

  /// The shear stiffness
  real64 m_shearStiffness;

  /// A threshold valued to determine whether a fracture is open or not.
  real64 m_displacementJumpThreshold;

  /// The aperture table function wrapper
  TableFunction::KernelWrapper m_apertureTable;
};


/**
 * @class ContactBase
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class ContactBase : public ConstitutiveBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  ContactBase( string const & name,
               Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ContactBase() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief accessor for penalty stiffness
   * @return the stiffness
   */
  real64 stiffness() const { return m_penaltyStiffness; }


  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ContactBaseUpdates;

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
    /// string/key for penalty stiffness
    static constexpr char const * penaltyStiffnessString() { return "penaltyStiffness"; }

    /// string/key for shear stiffness
    static constexpr char const * shearStiffnessString() { return "shearStiffness"; }

    /// string/key for the displacement jump threshold value
    static constexpr char const * displacementJumpThresholdString() { return "displacementJumpThreshold"; }

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


  /// The value of penalty to penetration
  real64 m_penaltyStiffness;

  /// The value of the shear stiffness
  real64 m_shearStiffness;

  /// The aperture tolerance to avoid floating point errors in expressions involving aperture
  real64 m_apertureTolerance;

  /// The name of the aperture table, if any
  string m_apertureTableName;

  /// A threshold valued to determine whether a fracture is open or not.
  real64 m_displacementJumpThreshold;

  /// Pointer to the function that limits the model aperture to a physically admissible value.
  TableFunction const * m_apertureTable;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 ContactBaseUpdates::computeHydraulicAperture( real64 const aperture,
                                                     real64 & dHydraulicAperture_dAperture ) const
{
  return m_apertureTable.compute( &aperture, &dHydraulicAperture_dAperture );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ContactBaseUpdates::addPressureToTraction( real64 const & pressure,
                                                arraySlice1d< real64 > const & tractionVector,
                                                real64 & dTraction_dPressure ) const
{
  tractionVector[0] -= pressure;
  dTraction_dPressure = -1.0;
}


} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_CONTACTBASE_HPP_ */
