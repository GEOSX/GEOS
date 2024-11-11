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
 * @file FrictionBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_FRICTIONBASE_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_FRICTIONBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"


namespace geos
{

namespace constitutive
{

/**
 * @class FrictionBaseUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class FrictionBaseUpdates
{
public:

  FrictionBaseUpdates( real64 const & displacementJumpThreshold )
    : m_displacementJumpThreshold( displacementJumpThreshold )
  {}

  /// Default copy constructor
  FrictionBaseUpdates( FrictionBaseUpdates const & ) = default;

  /// Default move constructor
  FrictionBaseUpdates( FrictionBaseUpdates && ) = default;

  /// Deleted default constructor
  FrictionBaseUpdates() = default;

  /// Deleted copy assignment operator
  FrictionBaseUpdates & operator=( FrictionBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  FrictionBaseUpdates & operator=( FrictionBaseUpdates && ) =  delete;

  /**
   * @brief Evaluate the traction vector and its derivatives wrt to pressure and jump
   * @param[in] dispJump the displacement jump
   * @param[in] fractureState the fracture state
   * @param[out] tractionVector the traction vector
   * @param[out] dTractionVector_dJump the derivative of the traction vector wrt displacement jump
   */
  GEOS_HOST_DEVICE
  inline
  virtual void computeShearTraction( localIndex const k,
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
                                    arraySlice1d< real64 const > const & oldDispJump,
                                    arraySlice1d< real64 const > const & tractionVector,
                                    real64 const pressure,
                                    integer & fractureState ) const
  { GEOS_UNUSED_VAR( k, dispJump, oldDispJump, tractionVector, pressure, fractureState ); }

  /**
   * @brief Update the trial traction vector ( return mapping )
   * @param[in] oldDispJump the displacement jump of the previous time step
   * @param[in] dispJump the displacement jump of the current time step
   * @param[in] penalty  the penalty coefficients
   * @param[in] traction the traction vector
   * @param[in] symmetric flag to compute only the symmetric part of dTraction_dDispJump
   * @param[in] fixedLimitTau flag to keep fixed the tangential stress
   * @param[in] normalTractionTolerance normal traction tolerance (if tn > tol => fracture is open)
   * @param[in] tangentialTractionTolerance tangential traction tolerance (if tau < tol => tau = 0)
   * @param[out] dTraction_dDispJump matrix containing the derivatives of traction over the displacement jump
   * @param[out] tractionNew the new traction vector
   * @param[out] fractureState the fracture state
   */
  GEOS_HOST_DEVICE
  inline
  virtual void updateTraction( arraySlice1d< real64 const > const & oldDispJump,
                               arraySlice1d< real64 const > const & dispJump,
                               arraySlice1d< real64 const > const & penalty,
                               arraySlice1d< real64 const > const & traction,
                               bool const symmetric,
                               bool const fixedLimitTau,
                               real64 const normalTractionTolerance,
                               real64 const tangentialTractionTolerance,
                               real64 ( & dTraction_dDispJump )[3][3],
                               real64 ( & tractionNew )[3],
                               integer & fractureState ) const
  {
    GEOS_UNUSED_VAR( oldDispJump, dispJump, penalty, traction, symmetric, fixedLimitTau,
                     normalTractionTolerance, tangentialTractionTolerance,
                     dTraction_dDispJump, tractionNew, fractureState );
  }

  /**
   * @brief Update the traction vector only ( return mapping )
   * @param[in] dispJump the displacement jump of the current time step
   * @param[in] deltaDispJump the delta displacement jump of the current time step
   * @param[in] penalty  the penalty coefficients
   * @param[in] traction the traction vector
   * @param[out] tractionNew the new traction vector
   */
  GEOS_HOST_DEVICE
  inline
  virtual void updateTractionOnly( arraySlice1d< real64 const > const & dispJump,
                                   arraySlice1d< real64 const > const & deltaDispJump,
                                   arraySlice1d< real64 const > const & penalty,
                                   arraySlice1d< real64 const > const & traction,
                                   arraySlice1d< real64 > const & tractionNew ) const
  { GEOS_UNUSED_VAR( dispJump, deltaDispJump, penalty, traction, tractionNew ); }

  /**
   * @brief Check for the constraint satisfaction
   * @param[in] dispJump the displacement jump of the current time step
   * @param[in] deltaDispJump the delta displacement jump of the current time step
   * @param[in] tractionVector the traction vector
   * @param[in] fractureState the fracture check
   * @param[in] tolerance the tolerance
   * @param[out] condConv flag indicating the result of the check
   */
  GEOS_HOST_DEVICE
  inline
  virtual void constraintCheck( arraySlice1d< real64 const > const & dispJump,
                                arraySlice1d< real64 const > const & deltaDispJump,
                                arraySlice1d< real64 > const & tractionVector,
                                integer const fractureState,
                                real64 const normalTractionTolerance,
                                real64 const normalDisplacementTolerance,
                                real64 const slidingTolerance,
                                real64 const slidingCheckTolerance,
                                integer & condConv ) const
  {
    GEOS_UNUSED_VAR( dispJump, deltaDispJump, tractionVector, fractureState,
                     normalTractionTolerance, normalDisplacementTolerance, slidingTolerance,
                     slidingCheckTolerance, condConv );
  }

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

  /// A threshold valued to determine whether a fracture is open or not.
  real64 m_displacementJumpThreshold;

};


/**
 * @class FrictionBase
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class FrictionBase : public ConstitutiveBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  FrictionBase( string const & name,
                Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~FrictionBase() override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FrictionBaseUpdates;

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
    /// string/key for the displacement jump threshold value
    static constexpr char const * displacementJumpThresholdString() { return "displacementJumpThreshold"; }
  };

protected:

  /// A threshold valued to determine whether a fracture is open or not.
  real64 m_displacementJumpThreshold;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_FRICTIONBASE_HPP_ */
