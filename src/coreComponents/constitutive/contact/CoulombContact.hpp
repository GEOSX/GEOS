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
  {
  }

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
  void computeTraction( localIndex const k,
                                arraySlice1d< real64 const > const & totalDispJump,
                                integer const & fractureState,
                                arraySlice1d< real64 > const & tractionVector) const;
  
  GEOS_HOST_DEVICE
  inline
  virtual void updateFractureState( localIndex const k,
                                    arraySlice1d< real64 const > const & dispJump,
                                    arraySlice1d< real64 const > const & tractionVector,
                                    integer & fractureState,
                                    real64 const pressure ) const override final;


  GEOS_HOST_DEVICE
  inline 
  void updateFractureState2( localIndex const k,
                             arraySlice1d< real64 const > const & dispJump,
                             arraySlice1d< real64 const > const & oldDispJump,
                             // TO DO tractionVector here is the combined force of contact
                             // fracture pressure, only  contact force is supposed to be here
                             arraySlice1d< real64 const > const & tractionVector,
                             integer & fractureState,
                             integer & oldFractureState,
                             real64 const pressure ) const override final;
  // GEOS_HOST_DEVICE
  // inline
  // void updateFractureStateUsingEffectiveTraction( localIndex const k,
  //                                                 arraySlice1d< real64 const > const & dispJump,
  //                                                 arraySlice1d< real64 const > const & tractionVector,
  //                                                 integer & fractureState ) const;

private:

  /// The cohesion for each upper level dimension (i.e. cell) of *this
  real64 m_cohesion;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  real64 m_frictionCoefficient;

  arrayView2d< real64 > m_elasticSlip;
  //array2d< real64 >  m_plasticSlip, m_slip;
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

// so many tractionVector with different definations, so confused!!!
// tractionVector = penalty*displJump[0]
GEOS_HOST_DEVICE
inline void CoulombContactUpdates::computeTraction( localIndex const k,
                                                    arraySlice1d< real64 const > const & oldDispJump,
                                                    arraySlice1d< real64 const > const & dispJump,
                                                    integer const & fractureState,
                                                    arraySlice1d< real64 > const & tractionVector,
                                                    arraySlice2d< real64 > const & dTractionVector_dJump ) const
{

  std::cout << "before dipljump "  << dispJump[0] << " " << dispJump[1] << " " << dispJump[2] << " " <<
  oldDispJump[0] << " "<< oldDispJump[1] << " " << oldDispJump[2]  << std::endl;
    
  // if (k == 0 || k == 79)
  // {
  //   std::cout << "In CoulombContactUpdates::computeTraction: " << std::endl;
  //   std::cout << "k = " << k << ", dispJump[0] = " << dispJump[0] << ", OldDispJump[0] = " << oldDispJump[0] << std::endl;
  // }

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
    // if (dispJump[0] < 0)
    // {
    //   tractionVector[0] = m_penaltyStiffness * dispJump[0];
    //   dTractionVector_dJump[0][0] = m_penaltyStiffness;
    // }
    // else
    // {
    //   tractionVector[0] = 0;
    //   dTractionVector_dJump[0][0] = 0;
    // }


    // derivative of the normal component w.r.t. to the nomral dispJump
    dTractionVector_dJump[0][0] = m_penaltyStiffness;

    switch( fractureState )
    {
      case fields::contact::FractureState::Stick:
      {
        // Elastic slip case
        // Tangential components of the traction are equal to tau
        // total displacement should be used here
        real64 const tau[2] = { m_shearStiffness * dispJump[1],
                                m_shearStiffness * dispJump[2]};

        tractionVector[1] = tau[0];
        tractionVector[2] = tau[1];

        dTractionVector_dJump[1][1] = m_shearStiffness;
        dTractionVector_dJump[2][2] = m_shearStiffness;

        // The slip is only elastic: we add the full slip to the elastic one
//        LvArray::tensorOps::copy< 2 >( m_elasticSlip[k], slip );


        break;
      }
      case fields::contact::FractureState::Slip:
      {
        // Plastic slip case
        real64 dLimitTau_dNormalTraction;
        real64 const limitTau = computeLimitTangentialTractionNorm( tractionVector[0],
                                                                    dLimitTau_dNormalTraction );

        //directional vector to decompose the tangential traction
        const real64 slip[2] = {dispJump[1], dispJump[2]};
        real64 const slipNorm = LvArray::tensorOps::l2Norm< 2 >( slip );

        // Tangential components of the traction computed based on the limitTau
        tractionVector[1] = limitTau * slip[0] / slipNorm;
        tractionVector[2] = limitTau * slip[1] / slipNorm;
        std::cout << "tractionVector [0][1][2] " << tractionVector[0] << " " << tractionVector[1] << " " << tractionVector[2] << " " <<
                                                    slip[0] << " " << slip[1] << " " << slipNorm << " " << limitTau << std::endl;

        dTractionVector_dJump[1][0] = -m_penaltyStiffness * dLimitTau_dNormalTraction * slip[0] / slipNorm;
        dTractionVector_dJump[1][1] = limitTau * pow( slip[1], 2 )  / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );
        dTractionVector_dJump[1][2] = -limitTau * slip[0] * slip[1] / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );

        dTractionVector_dJump[2][0] = -m_penaltyStiffness * dLimitTau_dNormalTraction * slip[1] / slipNorm;
        dTractionVector_dJump[2][1] = -limitTau * slip[0] * slip[1] / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );
        dTractionVector_dJump[2][2] = limitTau * pow( slip[0], 2 )  / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );

        // Compute elastic component of the slip for this case
//        real64 const plasticSlip[2] = { tractionVector[1] / m_shearStiffness,
//                                        tractionVector[2] / m_shearStiffness };
//
//        LvArray::tensorOps::copy< 2 >( m_elasticSlip[k], slip );
//        LvArray::tensorOps::subtract< 2 >( m_elasticSlip[k], plasticSlip );
        break;
      }
    }
  }
}

GEOS_HOST_DEVICE
inline void CoulombContactUpdates::computeTraction( localIndex const k,
                                                    arraySlice1d< real64 const > const & totalDispJump,
                                                    integer const & fractureState,
                                                    arraySlice1d< real64 > const & tractionVector) const
{

  bool const isOpen = fractureState == fields::contact::FractureState::Open;

  // Initialize everyting to 0
  tractionVector[0] = 0.0;
  tractionVector[1] = 0.0;
  tractionVector[2] = 0.00;
  if( !isOpen )
  {
    tractionVector[0] = m_penaltyStiffness * totalDispJump[0];

    switch( fractureState )
    {
      case fields::contact::FractureState::Stick:
      {
        // Elastic slip case
        // Tangential components of the traction are equal to tau
        // total displacement should be used here
        real64 const tau[2] = { m_shearStiffness * totalDispJump[1],
                                m_shearStiffness * totalDispJump[2]};
        tractionVector[1] = tau[0];
        tractionVector[2] = tau[1];
        break;
      }
      case fields::contact::FractureState::Slip:
      {
        // Plastic slip case
        real64 dLimitTau_dNormalTraction;
        real64 const limitTau = computeLimitTangentialTractionNorm( tractionVector[0],
                                                                    dLimitTau_dNormalTraction );

        //directional vector to decompose the tangential traction
        const real64 slip[2] = {totalDispJump[1], totalDispJump[2]};
        real64 const slipNorm = LvArray::tensorOps::l2Norm< 2 >( slip );

        // Tangential components of the traction computed based on the limitTau
        tractionVector[1] = limitTau * slip[0] / slipNorm;
        tractionVector[2] = limitTau * slip[1] / slipNorm;
        break;
      }
    }
  }
}
// so many tractionVector with different definations!!so confused!!
// tractionVector = penalty*displJump[0] - pressure
GEOS_HOST_DEVICE
inline void CoulombContactUpdates::updateFractureState( localIndex const k,
                                                        arraySlice1d< real64 const > const & dispJump,
                                                        arraySlice1d< real64 const > const & tractionVector,
                                                        integer & fractureState,
                                                        real64 const pressure ) const
{
  using namespace fields::contact;

  if (k == 0)
  {
    std::cout << "In CoulombContactUpdates::updateFractureState, pressure = " << pressure << std::endl;
    // std::cout << "k = " << k << ", dispJump[0] = " << dispJump[0] << std::endl;
  }

  if ( dispJump[0] >  -m_displacementJumpThreshold ) //  original
  {
    fractureState = FractureState::Open;
    m_elasticSlip[k][0] = 0.0;
    m_elasticSlip[k][1] = 0.0;

    std::cout << "k = " << k << ". normal_dispJump = " << dispJump[0] << ", fracture state = Open" << std::endl;
  }
  else
  {
    real64 const tau[2] = { tractionVector[1],
                            tractionVector[2] };
    // real64 const tauNorm = LvArray::tensorOps::l2Norm< 2 >( tau );
    real64 tauNorm = LvArray::tensorOps::l2Norm< 2 >( tau );

    // // add pressure to normal traction
    // tractionVector[0] += pressure;

    real64 dLimitTau_dNormalTraction;
    // real64 const limitTau = computeLimitTangentialTractionNorm( tractionVector[0],
    //                                                             dLimitTau_dNormalTraction );
    real64 const limitTau = computeLimitTangentialTractionNorm( tractionVector[0] + pressure, // pressure added to convert to effective traction, biotCoeff = 1 is assumed in fracture
                                                           dLimitTau_dNormalTraction );

    real64 slidingCheckTolerance = 0.05;
    if( fractureState == FractureState::Stick && tauNorm >= limitTau )
    {
      tauNorm *= (1.0 - slidingCheckTolerance);
    }
    else if( fractureState != FractureState::Stick && tauNorm <= limitTau )
    {
      tauNorm *= (1.0 + slidingCheckTolerance);
    }

    // Yield function (not necessary but makes it clearer)
    real64 const yield = tauNorm - limitTau;

    if (yield < 0)
    {
      // The slip is only elastic: we add the full slip to the elastic one
      std::cout << "k = " << k << ", normal traction = " << tractionVector[0] << ", pressure " << pressure << ", normal_dispJump = " << dispJump[0] << " " << dispJump[1] << " " << dispJump[2] << ", currentTau = " << tauNorm << ", limitTau =" << limitTau <<  ", fracture state = stick" << std::endl;
    }
    else
    {
      std::cout << "k = " << k << ", normal traction = " << tractionVector[0] << ", pressure " << pressure << ", normal_dispJump = " << dispJump[0] <<  " " << dispJump[1] << " " << dispJump[2] << ", currentTau = " << tauNorm << ", limitTau =" << limitTau << ", fracture state = slip" << std::endl;
    }
    fractureState = yield < 0 ? FractureState::Stick : FractureState::Slip;
  }
}

GEOS_HOST_DEVICE
inline void CoulombContactUpdates::updateFractureState2( localIndex const k,
                                                        arraySlice1d< real64 const > const & dispJump,
                                                        arraySlice1d< real64 const > const & oldDispJump,
                                                        // TO DO tractionVector here is the combined force of contact
                                                        // fracture pressure, only  contact force is supposed to be here
                                                        arraySlice1d< real64 const > const & tractionVector,
                                                        integer & fractureState,
                                                        integer & oldFractureState,
                                                        real64 const pressure ) const
{
  using namespace fields::contact;

  array1d< real64 > oldTractionVector(3);
  computeTraction(k, oldDispJump, oldFractureState, oldTractionVector);

  if ( dispJump[0] + oldDispJump[0] >  -m_displacementJumpThreshold ) //  original
  {
    fractureState = FractureState::Open;
    m_elasticSlip[k][0] = 0.0;
    m_elasticSlip[k][1] = 0.0;
  }
  else
  {
    real64 const tau[2] = { tractionVector[1] + oldTractionVector[1],
                            tractionVector[2] + oldTractionVector[2]};

    real64 tauNorm = LvArray::tensorOps::l2Norm< 2 >( tau );

    real64 dLimitTau_dNormalTraction;
    real64 const limitTau = computeLimitTangentialTractionNorm( oldTractionVector[0] + tractionVector[0] + pressure, // pressure added to convert to effective traction, biotCoeff = 1 is assumed in fracture
                                                           dLimitTau_dNormalTraction );

    real64 slidingCheckTolerance = 0.05;
    if( fractureState == FractureState::Stick && tauNorm >= limitTau )
    {
      tauNorm *= (1.0 - slidingCheckTolerance);
    }
    else if( fractureState != FractureState::Stick && tauNorm <= limitTau )
    {
      tauNorm *= (1.0 + slidingCheckTolerance);
    }

    // Yield function (not necessary but makes it clearer)
    real64 const yield = tauNorm - limitTau;

    if (yield < 0)
    {
      // The slip is only elastic: we add the full slip to the elastic one
      std::cout << "k = " << k << ", tractionVec = " << tractionVector[0] << " " << tractionVector[1] << " " << tractionVector[2] << ", pressure " << pressure << ", normal_dispJump = " << dispJump[0] << " " << dispJump[1] << " " << dispJump[2] << ", currentTau = " << tauNorm << ", limitTau =" << limitTau <<  ", fracture state = stick" << std::endl;
      std::cout << oldTractionVector[0] << " " << oldTractionVector[1] << " " << oldTractionVector[2] << std::endl;
    }
    else
    {
      std::cout << "k = " << k << ", tractionVec = " << tractionVector[0] << " " << tractionVector[1] << " " << tractionVector[2] << ", pressure " << pressure << ", normal_dispJump = " << dispJump[0] <<  " " << dispJump[1] << " " << dispJump[2] << ", currentTau = " << tauNorm << ", limitTau =" << limitTau << ", fracture state = slip, old state = " << oldFractureState << std::endl;
      std::cout << oldTractionVector[0] << " " << oldTractionVector[1] << " " << oldTractionVector[2] << std::endl;
    }
    fractureState = yield < 0 ? FractureState::Stick : FractureState::Slip;
  }
}
} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_COULOMBCONTACT_HPP_ */
