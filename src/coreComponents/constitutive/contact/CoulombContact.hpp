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
 *  @file CoulombContact.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_COULOMBCONTACT_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_COULOMBCONTACT_HPP_

#include "ContactBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class CoulombContactUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class CoulombContactUpdates : public ContactBaseUpdates
{
public:

  CoulombContactUpdates( real64 const & penaltyStiffness,
                         real64 const & shearStiffness,
                         real64 const & displacementJumpThreshold,
                         TableFunction const & apertureTable,
                         real64 const & cohesion,
                         real64 const & frictionCoefficient,
                         arrayView2d< real64 > const & elasticSlip )
    : ContactBaseUpdates( penaltyStiffness, shearStiffness, displacementJumpThreshold, apertureTable ),
    m_cohesion( cohesion ),
    m_frictionCoefficient( frictionCoefficient ),
    m_elasticSlip( elasticSlip )
  {}

  /// Default copy constructor
  CoulombContactUpdates( CoulombContactUpdates const & ) = default;

  /// Default move constructor
  CoulombContactUpdates( CoulombContactUpdates && ) = default;

  /// Deleted default constructor
  CoulombContactUpdates() = delete;

  /// Deleted copy assignment operator
  CoulombContactUpdates & operator=( CoulombContactUpdates const & ) = delete;

  /// Deleted move assignment operator
  CoulombContactUpdates & operator=( CoulombContactUpdates && ) =  delete;

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

private:

  /// The cohesion for each upper level dimension (i.e. cell) of *this
  real64 m_cohesion;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  real64 m_frictionCoefficient;

  arrayView2d< real64 > m_elasticSlip;
};


/**
 * @class CoulombContact
 *
 * Class to provide a CoulombContact friction model.
 */
class CoulombContact : public ContactBase
{
public:

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  CoulombContact( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~CoulombContact() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
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
  using KernelWrapper = CoulombContactUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

protected:

  virtual void postInputInitialization() override;

private:

  /// The cohesion for each upper level dimension (i.e. cell) of *this
  real64 m_cohesion;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  real64 m_frictionCoefficient;

  /// Elastic slip
  array2d< real64 > m_elasticSlip;

/**
 * @struct Set of "char const *" and keys for data specified in this class.
 */
  struct viewKeyStruct : public ContactBase::viewKeyStruct
  {
    /// string/key for cohesion
    static constexpr char const * cohesionString() { return "cohesion"; }

    /// string/key for friction coefficient
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }

    /// string/key for the elastic slip
    static constexpr char const * elasticSlipString() { return "elasticSlip"; }
  };

};


GEOS_HOST_DEVICE
real64 CoulombContactUpdates::computeLimitTangentialTractionNorm( real64 const & normalTraction,
                                                                  real64 & dLimitTangentialTractionNorm_dTraction ) const
{
  dLimitTangentialTractionNorm_dTraction = m_frictionCoefficient;
  return ( m_cohesion - normalTraction * m_frictionCoefficient );
}


GEOS_HOST_DEVICE
inline void CoulombContactUpdates::computeTraction( localIndex const k,
                                                    arraySlice1d< real64 const > const & oldDispJump,
                                                    arraySlice1d< real64 const > const & dispJump,
                                                    integer const & fractureState,
                                                    arraySlice1d< real64 > const & tractionVector,
                                                    arraySlice2d< real64 > const & dTractionVector_dJump ) const
{

  bool const isOpen = fractureState == fields::contact::FractureState::Open;

  // Initialize everyting to 0
  tractionVector[0] = 0.0;
  tractionVector[1] = 0.0;
  tractionVector[2] = 0.0;
  LvArray::forValuesInSlice( dTractionVector_dJump, []( real64 & val ){ val = 0.0; } );
  // If the fracture is open the traction is 0 and so are its derivatives so there is nothing to do
  if( !isOpen )
  {
    // normal component of the traction
    tractionVector[0] = m_penaltyStiffness * dispJump[0];
    // derivative of the normal component w.r.t. to the nomral dispJump
    dTractionVector_dJump[0][0] = m_penaltyStiffness;

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

        dTractionVector_dJump[1][0] = m_penaltyStiffness * dLimitTau_dNormalTraction * slip[0] / slipNorm;
        dTractionVector_dJump[1][1] = limitTau * pow( slip[1], 2 )  / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );
        dTractionVector_dJump[1][2] = limitTau * slip[0] * slip[1] / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );

        dTractionVector_dJump[2][0] = m_penaltyStiffness * dLimitTau_dNormalTraction * slip[1] / slipNorm;
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
}

GEOS_HOST_DEVICE
inline void CoulombContactUpdates::updateFractureState( localIndex const k,
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

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_COULOMBCONTACT_HPP_ */
