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
 * @file FrictionlessContact.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_FRICTIONLESSCONTACT_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_FRICTIONLESSCONTACT_HPP_

#include "constitutive/contact/ContactBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class FrictionlessContactUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class FrictionlessContactUpdates : public ContactBaseUpdates
{
public:

  FrictionlessContactUpdates( real64 const & penaltyStiffness,
                              real64 const & shearStiffness,
                              real64 const & displacementJumpThreshold,
                              TableFunction const & apertureTable )
    : ContactBaseUpdates( penaltyStiffness, shearStiffness, displacementJumpThreshold, apertureTable )
  {}

  /// Default copy constructor
  FrictionlessContactUpdates( FrictionlessContactUpdates const & ) = default;

  /// Default move constructor
  FrictionlessContactUpdates( FrictionlessContactUpdates && ) = default;

  /// Deleted default constructor
  FrictionlessContactUpdates() = default;

  /// Deleted copy assignment operator
  FrictionlessContactUpdates & operator=( FrictionlessContactUpdates const & ) = delete;

  /// Deleted move assignment operator
  FrictionlessContactUpdates & operator=( FrictionlessContactUpdates && ) =  delete;

  GEOS_HOST_DEVICE
  inline
  virtual void computeTraction( localIndex const k,
                                arraySlice1d< real64 const > const & oldDispJump,
                                arraySlice1d< real64 const > const & dispJump,
                                integer const & fractureState,
                                arraySlice1d< real64 > const & tractionVector,
                                arraySlice2d< real64 > const & dTractionVector_dJump ) const override final;


  GEOS_HOST_DEVICE
  inline
  virtual void updateFractureState( localIndex const k,
                                    arraySlice1d< real64 const > const & dispJump,
                                    arraySlice1d< real64 const > const & tractionVector,
                                    integer & fractureState ) const override final;


  /**
   * @brief Evaluate the limit tangential traction norm and return the derivative wrt normal traction
   * @param[in] normalTraction the normal traction
   * @param[out] dLimitTangentialTractionNorm_dTraction the derivative of the limit tangential traction norm wrt normal traction
   * @return the limit tangential traction norm
   */
  GEOS_HOST_DEVICE
  inline
  virtual real64 computeLimitTangentialTractionNorm( real64 const & normalTraction,
                                                     real64 & dLimitTangentialTractionNorm_dTraction ) const override final
  { GEOS_UNUSED_VAR( normalTraction, dLimitTangentialTractionNorm_dTraction ); return 0.0; }

private:
};


/**
 * @class FrictionlessContact
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class FrictionlessContact : public ContactBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  FrictionlessContact( string const & name,
                       Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~FrictionlessContact() override;

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return "FrictionlessContact"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FrictionlessContactUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @struct Structure to hold scoped key names
   */
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {};

protected:

};

GEOS_HOST_DEVICE
inline void FrictionlessContactUpdates::computeTraction( localIndex const k,
                                                         arraySlice1d< real64 const > const & oldDispJump,
                                                         arraySlice1d< real64 const > const & dispJump,
                                                         integer const & fractureState,
                                                         arraySlice1d< real64 > const & tractionVector,
                                                         arraySlice2d< real64 > const & dTractionVector_dJump ) const
{
  GEOS_UNUSED_VAR( k, oldDispJump );

  bool const isOpen = fractureState == fields::contact::FractureState::Open;

  tractionVector[0] = isOpen ? 0.0 : m_penaltyStiffness * dispJump[0];
  tractionVector[1] = 0.0;
  tractionVector[2] = 0.0;

  LvArray::forValuesInSlice( dTractionVector_dJump, []( real64 & val ){ val = 0.0; } );
  dTractionVector_dJump( 0, 0 ) = isOpen ? 0.0 : m_penaltyStiffness;
}

GEOS_HOST_DEVICE
inline void FrictionlessContactUpdates::updateFractureState( localIndex const k,
                                                             arraySlice1d< real64 const > const & dispJump,
                                                             arraySlice1d< real64 const > const & tractionVector,
                                                             integer & fractureState ) const
{
  GEOS_UNUSED_VAR( k, tractionVector );
  using namespace fields::contact;
  fractureState = dispJump[0] > m_displacementJumpThreshold ? FractureState::Open : FractureState::Stick;
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_FRICTIONLESSCONTACT_HPP_ */
