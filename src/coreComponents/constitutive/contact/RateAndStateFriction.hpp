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
 *  @file RateAndStateFriction.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_RATEANDSTATEFRICTION_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_RATEANDSTATEFRICTION_HPP_

#include "FrictionBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class RateAndStateFrictionUpdates
 *
 * This class is used for in-kernel contact relation updates
 */


/**
 * @class RateAndStateFriction
 *
 * Class to provide a RateAndStateFriction friction model.
 */
class RateAndStateFriction : public FrictionBase
{
public:

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  RateAndStateFriction( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~RateAndStateFriction() override;

  static string catalogName() { return "RateAndStateFriction"; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override final;

  class KernelWrapper : public FrictionBaseUpdates
  {
public:
    KernelWrapper( real64 const displacementJumpThreshold,
                   arrayView1d< real64 > frictionCoefficient,
                   arrayView1d< real64 const > a,
                   arrayView1d< real64 const > b,
                   arrayView1d< real64 const > Dc,
                   arrayView1d< real64 const > V0,
                   arrayView1d< real64 const > mu0 )
      : FrictionBaseUpdates( displacementJumpThreshold ),
      m_frictionCoefficient( frictionCoefficient ),
      m_a( a ),
      m_b( b ),
      m_Dc( Dc ),
      m_V0( V0 ),
      m_mu0( mu0 )
    {}

    /// Default copy constructor
    KernelWrapper( KernelWrapper const & ) = default;

    /// Default move constructor
    KernelWrapper( KernelWrapper && ) = default;

    /// Deleted default constructor
    KernelWrapper() = delete;

    /// Deleted copy assignment operator
    KernelWrapper & operator=( KernelWrapper const & ) = delete;

    /// Deleted move assignment operator
    KernelWrapper & operator=( KernelWrapper && ) =  delete;

    GEOS_HOST_DEVICE
    real64 getACoefficient( localIndex const k ) const { return m_a[k]; }

    GEOS_HOST_DEVICE
    real64 getBCoefficient( localIndex const k ) const { return m_b[k]; }

    GEOS_HOST_DEVICE
    real64 getDcCoefficient( localIndex const k ) const { return m_Dc[k]; }

    GEOS_HOST_DEVICE
    inline
    virtual void updateFractureState( localIndex const k,
                                      arraySlice1d< real64 const > const & dispJump,
                                      arraySlice1d< real64 const > const & oldDispJump,
                                      arraySlice1d< real64 const > const & tractionVector,
                                      real64 const pressure,
                                      integer & fractureState ) const override final;

    GEOS_HOST_DEVICE
    inline real64 frictionCoefficient( localIndex const k,
                                       real64 const slipRate,
                                       real64 const stateVariable ) const;

    GEOS_HOST_DEVICE
    inline real64 dFrictionCoefficient_dSlipRate( localIndex const k,
                                                  real64 const slipRate,
                                                  real64 const stateVariable ) const;

    GEOS_HOST_DEVICE
    inline real64 dFrictionCoefficient_dStateVariable( localIndex const k,
                                                       real64 const slipRate,
                                                       real64 const stateVariable ) const;

    GEOS_HOST_DEVICE
    inline real64 stateEvolution( localIndex const k,
                                  real64 const slipRate,
                                  real64 const stateVariable ) const;

    GEOS_HOST_DEVICE
    inline real64 dStateEvolution_dStateVariable( localIndex const k,
                                                  real64 const slipRate,
                                                  real64 const stateVariable ) const;

    GEOS_HOST_DEVICE
    inline real64 dStateEvolution_dSlipRate( localIndex const k,
                                             real64 const slipRate,
                                             real64 const stateVariable ) const;
private:
    /// The friction coefficient
    arrayView1d< real64 > m_frictionCoefficient;

    /// Rate and State coefficient a
    arrayView1d< real64 const > m_a;

    /// Rate and State coefficient b
    arrayView1d< real64 const > m_b;

    /// Rate and State characteristic length
    arrayView1d< real64 const > m_Dc;

    /// Rate and State reference velocity
    arrayView1d< real64 const > m_V0;

    /// Rate and State reference friction coefficient
    arrayView1d< real64 const > m_mu0;
  };


  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates() const;

private:

  virtual void postInputInitialization() override;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_frictionCoefficient;

  /// Rate and State coefficient a
  array1d< real64 > m_a;

  /// Rate and State coefficient b
  array1d< real64 > m_b;

  /// Rate and State characteristic length
  array1d< real64 > m_Dc;

  /// Rate and State reference velocity
  array1d< real64 > m_V0;

  /// Rate and State reference friction coefficient
  array1d< real64 > m_mu0;

  ///  Default value of Rate and State coefficient a
  real64 m_defaultA;
  /// Default value of Rate and State coefficient b
  real64 m_defaultB;

  ///  Default value of Rate and State characteristic length
  real64 m_defaultDc;

  ///  Default value of Rate and State reference velocity
  real64 m_defaultV0;

  /// Default value of Rate and State reference friction coefficient
  real64 m_defaultMu0;

/**
 * @struct Set of "char const *" and keys for data specified in this class.
 */
  struct viewKeyStruct : public FrictionBase::viewKeyStruct
  {
    /// string/key for friction coefficient
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }
    /// string/key for Rate and State coefficient a
    static constexpr char const * aCoefficientString() { return "a"; }
    /// string/key for Rate and State coefficient b
    static constexpr char const * bCoefficientString() { return "b"; }
    /// string/key for Rate and State characteristic length
    static constexpr char const * DcCoefficientString() { return "Dc"; }
    /// string/key for reference slip rate
    static constexpr char const * referenceVelocityString() { return "referenceVelocity"; }
    /// string/key for reference friction coefficient
    static constexpr char const * referenceFrictionCoefficientString() { return "referenceFrictionCoefficient"; }
    /// string/key for the default value of Rate and State coefficient a
    static constexpr char const * defaultACoefficientString() { return "defaultA"; }
    /// string/key for the default value of Rate and State coefficient b
    static constexpr char const * defaultBCoefficientString() { return "defaultB"; }
    /// string/key for the default value of Rate and State characteristic length
    static constexpr char const * defaultDcCoefficientString() { return "defaultDc"; }
    /// string/key for the default value ofreference slip rate
    static constexpr char const * defaultReferenceVelocityString() { return "defaultReferenceVelocity"; }
    /// string/key for the default value of reference friction coefficient
    static constexpr char const * defaultReferenceFrictionCoefficientString() { return "defaultReferenceFrictionCoefficient"; }
  };

};

GEOS_HOST_DEVICE
inline void RateAndStateFriction::KernelWrapper::updateFractureState( localIndex const k,
                                                                      arraySlice1d< real64 const > const & dispJump,
                                                                      arraySlice1d< real64 const > const & oldDispJump,
                                                                      arraySlice1d< real64 const > const & tractionVector,
                                                                      real64 const pressure,
                                                                      integer & fractureState ) const
{

  GEOS_UNUSED_VAR( tractionVector, k, oldDispJump, pressure );
  using namespace fields::contact;

  if( dispJump[0] >  -m_displacementJumpThreshold )
  {
    fractureState = FractureState::Open;
  }
  else
  {
    fractureState = FractureState::Slip;
  }
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::frictionCoefficient( localIndex const k,
                                                                        real64 const slipRate,
                                                                        real64 const stateVariable ) const
{

  real64 const arg = ( slipRate / (2 * m_V0[k]) ) * LvArray::math::exp( stateVariable / m_a[k] );
  m_frictionCoefficient[k]  = m_a[k] * LvArray::math::asinh( arg );

  return m_frictionCoefficient[k];
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dFrictionCoefficient_dSlipRate( localIndex const k,
                                                                                   real64 const slipRate,
                                                                                   real64 const stateVariable ) const
{

  real64 const arg = ( slipRate / (2 * m_V0[k]) ) * LvArray::math::exp( stateVariable / m_a[k] );

  return ( m_a[k] * LvArray::math::exp( stateVariable / m_a[k] ) ) / (2 * m_V0[k] * LvArray::math::sqrt( 1 + arg * arg ));
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dFrictionCoefficient_dStateVariable( localIndex const k,
                                                                                        real64 const slipRate,
                                                                                        real64 const stateVariable ) const
{

  real64 const arg = ( slipRate / (2 * m_V0[k]) ) * LvArray::math::exp( stateVariable / m_a[k] );

  return ( slipRate * LvArray::math::exp( stateVariable / m_a[k] ) ) / (2 * m_V0[k] * LvArray::math::sqrt( 1 + arg * arg ));
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::stateEvolution( localIndex const k,
                                                                   real64 const slipRate,
                                                                   real64 const stateVariable ) const
{
  real64 const mu = frictionCoefficient( k, slipRate, stateVariable );

  return -slipRate / m_Dc[k] * (mu - m_mu0[k] + (m_b[k] - m_a[k]) * LvArray::math::log( slipRate / m_V0[k] ));
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dStateEvolution_dStateVariable( localIndex const k,
                                                                                   real64 const slipRate,
                                                                                   real64 const stateVariable ) const
{
  return -slipRate / m_Dc[k] * dFrictionCoefficient_dStateVariable( k, slipRate, stateVariable );
}

GEOS_HOST_DEVICE
inline real64 RateAndStateFriction::KernelWrapper::dStateEvolution_dSlipRate( localIndex const k,
                                                                              real64 const slipRate,
                                                                              real64 const stateVariable ) const
{
  real64 const part1 =  frictionCoefficient( k, slipRate, stateVariable ) - m_mu0[k] + (m_b[k] - m_a[k]) * LvArray::math::log( slipRate / m_V0[k] );

  real64 const part2 = dFrictionCoefficient_dSlipRate( k, slipRate, stateVariable ) * slipRate + (m_b[k] - m_a[k]);

  return -1.0 / m_Dc[k] * ( part1 + part2 );
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_RATEANDSTATEFRICTION_HPP_ */
