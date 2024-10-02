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
 *  @file CoulombFriction.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_COULOMBFRICTION_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_COULOMBFRICTION_HPP_

#include "FrictionBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class CoulombFrictionUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class CoulombFrictionUpdates : public FrictionBaseUpdates
{
public:
  CoulombFrictionUpdates( real64 const & displacementJumpThreshold,
                          real64 const & shearStiffness,
                          real64 const & cohesion,
                          real64 const & frictionCoefficient,
                          arrayView2d< real64 > const & elasticSlip )
    : FrictionBaseUpdates( displacementJumpThreshold ),
    m_shearStiffness( shearStiffness ),
    m_cohesion( cohesion ),
    m_frictionCoefficient( frictionCoefficient ),
    m_elasticSlip( elasticSlip )
  {}

  /// Default copy constructor
  CoulombFrictionUpdates( CoulombFrictionUpdates const & ) = default;

  /// Default move constructor
  CoulombFrictionUpdates( CoulombFrictionUpdates && ) = default;

  /// Deleted default constructor
  CoulombFrictionUpdates() = delete;

  /// Deleted copy assignment operator
  CoulombFrictionUpdates & operator=( CoulombFrictionUpdates const & ) = delete;

  /// Deleted move assignment operator
  CoulombFrictionUpdates & operator=( CoulombFrictionUpdates && ) =  delete;

  /**
   * @brief Evaluate the limit tangential traction norm and return the derivative wrt normal traction
   * @param[in] normalTraction the normal traction
   * @param[out] dLimitTangentialTractionNorm_dTraction the derivative of the limit tangential traction norm wrt normal traction
   * @return the limit tangential traction norm
   */
  GEOS_HOST_DEVICE
  inline
  virtual real64 computeLimitTangentialTractionNorm( real64 const & normalTraction,
                                                     real64 & dLimitTangentialTractionNorm_dTraction ) const override final;

  GEOS_HOST_DEVICE
  inline
  virtual void computeShearTraction( localIndex const k,
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
                               real64 ( &dTraction_dDispJump )[3][3],
                               real64 ( &tractionNew )[3],
                               integer & fractureState ) const override final;


  GEOS_HOST_DEVICE
  inline
  virtual void updateTractionOnly( arraySlice1d< real64 const > const & dispJump,
                                   arraySlice1d< real64 const > const & deltaDispJump,
                                   arraySlice1d< real64 const > const & penalty,
                                   arraySlice1d< real64 const > const & traction,
                                   arraySlice1d< real64 > const & tractionNew ) const override final;

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
                                integer & condConv ) const override final;

private:
  /// The shear stiffness
  real64 m_shearStiffness;

  /// The cohesion for each upper level dimension (i.e. cell) of *this
  real64 m_cohesion;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  real64 m_frictionCoefficient;

  arrayView2d< real64 > m_elasticSlip;
};


/**
 * @class CoulombFriction
 *
 * Class to provide a CoulombFriction friction model.
 */
class CoulombFriction : public FrictionBase
{
public:

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  CoulombFriction( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~CoulombFriction() override;

  static string catalogName() { return "Coulomb"; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override final;

  /**
   * @brief Const accessor for cohesion
   * @return A const reference to arrayView1d<real64 const> containing the
   *         cohesions (at every element).
   */
  real64 const & cohesion() const { return m_cohesion; }

  /**
   * @brief Const accessor for friction angle
   * @return A const reference to arrayView1d<real64 const> containing the
   *         friction coefficient (at every element).
   */
  real64 const & frictionCoefficient() const { return m_frictionCoefficient; }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CoulombFrictionUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates() const;

protected:

  virtual void postInputInitialization() override;

private:

  /// The shear stiffness
  real64 m_shearStiffness;

  /// The cohesion for each upper level dimension (i.e. cell) of *this
  real64 m_cohesion;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  real64 m_frictionCoefficient;

  /// Elastic slip
  array2d< real64 > m_elasticSlip;

/**
 * @struct Set of "char const *" and keys for data specified in this class.
 */
  struct viewKeyStruct : public FrictionBase::viewKeyStruct
  {
    /// string/key for shear stiffness
    static constexpr char const * shearStiffnessString() { return "shearStiffness"; }

    /// string/key for cohesion
    static constexpr char const * cohesionString() { return "cohesion"; }

    /// string/key for friction coefficient
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }

    /// string/key for the elastic slip
    static constexpr char const * elasticSlipString() { return "elasticSlip"; }
  };

};


GEOS_HOST_DEVICE
real64 CoulombFrictionUpdates::computeLimitTangentialTractionNorm( real64 const & normalTraction,
                                                                   real64 & dLimitTangentialTractionNorm_dTraction ) const
{
  dLimitTangentialTractionNorm_dTraction = m_frictionCoefficient;
  return ( m_cohesion - normalTraction * m_frictionCoefficient );
}


GEOS_HOST_DEVICE
inline void CoulombFrictionUpdates::computeShearTraction( localIndex const k,
                                                          arraySlice1d< real64 const > const & oldDispJump,
                                                          arraySlice1d< real64 const > const & dispJump,
                                                          integer const & fractureState,
                                                          arraySlice1d< real64 > const & tractionVector,
                                                          arraySlice2d< real64 > const & dTractionVector_dJump ) const
{
  // Compute the slip
  real64 const slip[2] = { dispJump[1] - oldDispJump[1],
                           dispJump[2] - oldDispJump[2] };


  real64 const tau[2] = { m_shearStiffness * ( slip[0] + m_elasticSlip[k][0] ),
                          m_shearStiffness * ( slip[1] + m_elasticSlip[k][1] ) };

  switch( fractureState )
  {
    case fields::contact::FractureState::Stick:
    {
      // Elastic slip case
      // Tangential components of the traction are equal to tau
      tractionVector[1] = tau[0];
      tractionVector[2] = tau[1];

      dTractionVector_dJump[1][1] = m_shearStiffness;
      dTractionVector_dJump[2][2] = m_shearStiffness;

      // The slip is only elastic: we add the full slip to the elastic one
      LvArray::tensorOps::add< 2 >( m_elasticSlip[k], slip );

      break;
    }
    case fields::contact::FractureState::Slip:
    {
      // Plastic slip case
      real64 dLimitTau_dNormalTraction;
      real64 const limitTau = computeLimitTangentialTractionNorm( tractionVector[0],
                                                                  dLimitTau_dNormalTraction );

      real64 const slipNorm = LvArray::tensorOps::l2Norm< 2 >( slip );

      // Tangential components of the traction computed based on the limitTau
      tractionVector[1] = limitTau * slip[0] / slipNorm;
      tractionVector[2] = limitTau * slip[1] / slipNorm;

      dTractionVector_dJump[1][0] = dTractionVector_dJump[0][0] * dLimitTau_dNormalTraction * slip[0] / slipNorm;
      dTractionVector_dJump[1][1] = limitTau * pow( slip[1], 2 )  / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );
      dTractionVector_dJump[1][2] = limitTau * slip[0] * slip[1] / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );

      dTractionVector_dJump[2][0] = dTractionVector_dJump[0][0] * dLimitTau_dNormalTraction * slip[1] / slipNorm;
      dTractionVector_dJump[2][1] = limitTau * slip[0] * slip[1] / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );
      dTractionVector_dJump[2][2] = limitTau * pow( slip[0], 2 )  / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );

      // Compute elastic component of the slip for this case
      real64 const plasticSlip[2] = { tractionVector[1] / m_shearStiffness,
                                      tractionVector[2] / m_shearStiffness };

      LvArray::tensorOps::copy< 2 >( m_elasticSlip[k], slip );
      LvArray::tensorOps::subtract< 2 >( m_elasticSlip[k], plasticSlip );

      break;
    }
  }
}

GEOS_HOST_DEVICE
inline void CoulombFrictionUpdates::updateFractureState( localIndex const k,
                                                         arraySlice1d< real64 const > const & dispJump,
                                                         arraySlice1d< real64 const > const & tractionVector,
                                                         integer & fractureState ) const
{
  using namespace fields::contact;

  if( dispJump[0] >  -m_displacementJumpThreshold )
  {
    fractureState = FractureState::Open;
    m_elasticSlip[k][0] = 0.0;
    m_elasticSlip[k][1] = 0.0;
  }
  else
  {
    real64 const tau[2] = { tractionVector[1],
                            tractionVector[2] };
    real64 const tauNorm = LvArray::tensorOps::l2Norm< 2 >( tau );

    real64 dLimitTau_dNormalTraction;
    real64 const limitTau = computeLimitTangentialTractionNorm( tractionVector[0],
                                                                dLimitTau_dNormalTraction );

    // Yield function (not necessary but makes it clearer)
    real64 const yield = tauNorm - limitTau;

    fractureState = yield < 0 ? FractureState::Stick : FractureState::Slip;
  }
}

GEOS_HOST_DEVICE
inline void CoulombFrictionUpdates::updateTraction( arraySlice1d< real64 const > const & oldDispJump,
                                                    arraySlice1d< real64 const > const & dispJump,
                                                    arraySlice1d< real64 const > const & penalty,
                                                    arraySlice1d< real64 const > const & traction,
                                                    bool const symmetric,
                                                    bool const fixedLimitTau,
                                                    real64 const normalTractionTolerance,
                                                    real64 const tangentialTractionTolerance,
                                                    real64 ( & dTraction_dDispJump )[3][3],
                                                    real64 ( & tractionNew ) [3],
                                                    integer & fractureState ) const
{

  using namespace fields::contact;

  real64 dLimitTangentialTractionNorm_dTraction = 0.0;
  real64 limitTau = 0.0;

  // Compute the trial traction
  real64 tractionTrial[ 3 ];
  tractionTrial[ 0 ] = traction[0] + penalty[0] * dispJump[0];
  tractionTrial[ 1 ] = traction[1] + penalty[1] * (dispJump[1] - oldDispJump[1]);
  tractionTrial[ 2 ] = traction[2] + penalty[1] * (dispJump[2] - oldDispJump[2]);

  // Compute tangential trial traction norm
  real64 const tau[2] = { tractionTrial[1],
                          tractionTrial[2] };
  real64 const tractionTrialNorm = LvArray::tensorOps::l2Norm< 2 >( tau );

  // If normal tangential trial is positive (opening)
  fractureState = FractureState::Stick;
  if( tractionTrial[ 0 ] > normalTractionTolerance )
  {
    tractionNew[0] = 0.0;
    dTraction_dDispJump[0][0] = 0.0;
    fractureState = FractureState::Open;
  }
  else
  {
    tractionNew[0] = tractionTrial[0];
    dTraction_dDispJump[0][0] = -penalty[ 0 ];
  }

  // Compute limit Tau
  if( fixedLimitTau )
  {
    limitTau = computeLimitTangentialTractionNorm( traction[0],
                                                   dLimitTangentialTractionNorm_dTraction );
  }
  else
  {
    limitTau = computeLimitTangentialTractionNorm( tractionNew[0],
                                                   dLimitTangentialTractionNorm_dTraction );
  }

  if( tractionTrialNorm <= tangentialTractionTolerance )
  {
    // It is needed for the first iteration (both t and u are equal to zero)
    dTraction_dDispJump[1][1] = -penalty[1];
    dTraction_dDispJump[2][2] = -penalty[1];

    tractionNew[1] = tractionTrial[1];
    tractionNew[2] = tractionTrial[2];

    if( fractureState != FractureState::Open )
      fractureState =  FractureState::Stick;
  }
  else if( limitTau <= tangentialTractionTolerance )
  {
    dTraction_dDispJump[1][1] = 0.0;
    dTraction_dDispJump[2][2] = 0.0;

    tractionNew[1] = (fixedLimitTau) ? tractionTrial[1] : 0.0;
    tractionNew[2] = (fixedLimitTau) ? tractionTrial[2] : 0.0;

    if( fractureState != FractureState::Open )
      fractureState =  FractureState::Slip;
  }
  else
  {
    // Compute psi and dpsi
    //real64 const psi = std::tanh( tractionTrialNorm/limitTau );
    //real64 const dpsi = 1.0-std::pow(psi,2);
    real64 const psi = ( tractionTrialNorm > limitTau) ? 1.0 : tractionTrialNorm/limitTau;
    real64 const dpsi = ( tractionTrialNorm > limitTau) ? 0.0 : 1.0;

    if( fractureState != FractureState::Open )
    {
      fractureState = ( tractionTrialNorm > limitTau) ? FractureState::Slip : FractureState::Stick;
    }

    // Two symmetric 2x2 matrices
    real64 dNormTTdgT[ 3 ];
    dNormTTdgT[ 0 ] = tractionTrial[ 1 ] * tractionTrial[ 1 ];
    dNormTTdgT[ 1 ] = tractionTrial[ 2 ] * tractionTrial[ 2 ];
    dNormTTdgT[ 2 ] = tractionTrial[ 1 ] * tractionTrial[ 2 ];

    real64 dTdgT[ 3 ];
    dTdgT[ 0 ] = (tractionTrialNorm * tractionTrialNorm - dNormTTdgT[0]);
    dTdgT[ 1 ] = (tractionTrialNorm * tractionTrialNorm - dNormTTdgT[1]);
    dTdgT[ 2 ] = -dNormTTdgT[2];

    LvArray::tensorOps::scale< 3 >( dNormTTdgT, 1. / std::pow( tractionTrialNorm, 2 ) );
    LvArray::tensorOps::scale< 3 >( dTdgT, 1. / std::pow( tractionTrialNorm, 3 )  );

    // Compute dTdDispJump
    dTraction_dDispJump[1][1] = -penalty[1] * (
      dpsi * dNormTTdgT[0] +
      psi * dTdgT[0] * limitTau );
    dTraction_dDispJump[2][2] = -penalty[1] * (
      dpsi * dNormTTdgT[1] +
      psi * dTdgT[1] * limitTau );
    dTraction_dDispJump[1][2] = -penalty[1] * (
      dpsi * dNormTTdgT[2] +
      psi * dTdgT[2] * limitTau );
    dTraction_dDispJump[2][1] = dTraction_dDispJump[1][2];

    if( !symmetric )
    {
      // Nonsymetric term
      dTraction_dDispJump[1][0] = -dTraction_dDispJump[0][0] * m_frictionCoefficient *
                                  tractionTrial[1] * (psi/tractionTrialNorm - dpsi/limitTau);
      dTraction_dDispJump[2][0] = -dTraction_dDispJump[0][0] * m_frictionCoefficient  *
                                  tractionTrial[2] * (psi/tractionTrialNorm - dpsi/limitTau);
    }

    LvArray::tensorOps::scale< 3 >( tractionTrial, (psi * limitTau)/tractionTrialNorm );
    tractionNew[1] = tractionTrial[1];
    tractionNew[2] = tractionTrial[2];
  }

}

GEOS_HOST_DEVICE
inline void CoulombFrictionUpdates::updateTractionOnly( arraySlice1d< real64 const > const & dispJump,
                                                        arraySlice1d< real64 const > const & deltaDispJump,
                                                        arraySlice1d< real64 const > const & penalty,
                                                        arraySlice1d< real64 const > const & traction,
                                                        arraySlice1d< real64 > const & tractionNew ) const
{

  // TODO: Pass this tol as an argument or define a new class member
  real64 const zero = LvArray::NumericLimits< real64 >::epsilon;

  tractionNew[0] = traction[0] + penalty[0] * dispJump[0];
  tractionNew[1] = traction[1] + penalty[1] * deltaDispJump[1];
  tractionNew[2] = traction[2] + penalty[1] * deltaDispJump[2];

  real64 const tau[2] = { tractionNew[1],
                          tractionNew[2] };
  real64 const currentTau = LvArray::tensorOps::l2Norm< 2 >( tau );

  real64 dLimitTangentialTractionNorm_dTraction = 0.0;
  real64 const limitTau = computeLimitTangentialTractionNorm( tractionNew[0],
                                                              dLimitTangentialTractionNorm_dTraction );

  // Compute psi
  real64 psi;
  if( limitTau < zero )
  {
    psi = 1.0;
  }
  else
  {
    //psi = std::tanh(currentTau / limitTau);
    psi = (currentTau > limitTau ) ? 1.0 : currentTau/limitTau;
  }

  // Compute the new tangential traction
  if( limitTau > zero && currentTau > zero )
  {
    tractionNew[1] *= limitTau * psi / currentTau;
    tractionNew[2] *= limitTau * psi / currentTau;
  }
  else
  {
    tractionNew[1] *= 0.0;
    tractionNew[2] *= 0.0;
  }
}

GEOS_HOST_DEVICE
inline void CoulombFrictionUpdates::constraintCheck( arraySlice1d< real64 const > const & dispJump,
                                                     arraySlice1d< real64 const > const & deltaDispJump,
                                                     arraySlice1d< real64 > const & tractionVector,
                                                     integer const fractureState,
                                                     real64 const normalTractionTolerance,
                                                     real64 const normalDisplacementTolerance,
                                                     real64 const slidingTolerance,
                                                     real64 const slidingCheckTolerance,
                                                     integer & condConv ) const
{

  using namespace fields::contact;

  // Compute the slip
  real64 const deltaDisp[2] = { deltaDispJump[1],
                                deltaDispJump[2] };

  real64 const deltaDispNorm = LvArray::tensorOps::l2Norm< 2 >( deltaDisp );

  // Compute current Tau and limit Tau
  real64 const tau[2] = { tractionVector[1],
                          tractionVector[2] };
  real64 const currentTau = LvArray::tensorOps::l2Norm< 2 >( tau );

  real64 dLimitTangentialTractionNorm_dTraction = 0.0;
  real64 const limitTau = computeLimitTangentialTractionNorm( tractionVector[0],
                                                              dLimitTangentialTractionNorm_dTraction );

  condConv = 0;
  // Case 1: if it is open
  if( tractionVector[0] >= normalTractionTolerance )
  {
    if( fractureState != FractureState::Open )
    {
      condConv = 1;
    }
    tractionVector[0] = 0.0;
    tractionVector[1] = 0.0;
    tractionVector[2] = 0.0;
  }
  else
  {
    // Case 2: compenetration
    if(( LvArray::math::abs( dispJump[0] ) > normalDisplacementTolerance ) &&
       (fractureState != FractureState::Open))
    {
      condConv = 2;
    }
    // Case 3: it is stick and dg is greater than 0
    if( fractureState == FractureState::Stick &&
        deltaDispNorm > slidingTolerance )
    {
      condConv = 3;
    }

    // Case 4: the elastic tangential traction is greater than the limit
    if( currentTau > (LvArray::math::abs( limitTau ) * (1.0 + slidingCheckTolerance)) )
    {
      condConv = 4;
    }
  }
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_COULOMBFRICTION_HPP_ */
