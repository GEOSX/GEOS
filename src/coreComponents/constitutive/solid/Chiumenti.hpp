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
 * @file Chiumenti.hpp
 * @brief Linear hardening model for M. Cervera, M Chiumenti (2006)
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_CHIUMENTI_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_CHIUMENTI_HPP_

#include "HyperelasticMMS.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class ChiumentiUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class ChiumentiUpdates : public HyperelasticMMSUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] damage The ArrayView holding the damage for each quardrature point.
   * @param[in] lengthScale The ArrayView holding the length scale for each element.
   * @param[in] strengthScale The ArrayView holding the strength scale for each element/particle.
   * @param[in] criticalLength The critical length.
   * @param[in] failureStrength The failure strength.
   * @param[in] energyReleaseRate The energy release rate.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  ChiumentiUpdates( arrayView2d< real64 > const & damage,
                        arrayView1d< real64 > const & lengthScale,
                        arrayView1d< real64 > const & strengthScale,
                        real64 const & criticalLength,
                        real64 const & failureStrength,
                        real64 const & energyReleaseRate,
                        arrayView1d< real64 const > const & lambda,
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView1d< real64 const > const & thermalExpansionCoefficient,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                        arrayView2d< real64 > const & density,
                        arrayView2d< real64 > const & wavespeed,
                        bool const & disableInelasticity ):
    HyperelasticMMSUpdates( lambda,
                            shearModulus,
                            thermalExpansionCoefficient,
                            newStress,
                            oldStress,
                            density,
                            wavespeed,
                            disableInelasticity ),
    m_damage( damage ),
    m_lengthScale( lengthScale ),
    m_strengthScale( strengthScale ),
    m_criticalLength( criticalLength ),
    m_failureStrength( failureStrength ),
    m_energyReleaseRate( energyReleaseRate )
  {}

  /// Default copy constructor
  ChiumentiUpdates( ChiumentiUpdates const & ) = default;

  /// Default move constructor
  ChiumentiUpdates( ChiumentiUpdates && ) = default;

  /// Deleted default constructor
  ChiumentiUpdates() = delete;

  /// Deleted copy assignment operator
  ChiumentiUpdates & operator=( ChiumentiUpdates const & ) = delete;

  /// Deleted move assignment operator
  ChiumentiUpdates & operator=( ChiumentiUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using HyperelasticMMSUpdates::smallStrainUpdate;
  
  GEOS_HOST_DEVICE
  void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const override final;

  GEOS_HOST_DEVICE
  void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         real64 ( &stiffness )[6][6] ) const override final;

  GEOS_HOST_DEVICE
  void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         DiscretizationOps & stiffness ) const;

  GEOS_HOST_DEVICE
  void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( & beginningRotation )[3][3],
                                             real64 const ( & endRotation )[3][3],
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 ( &stiffness )[6][6] ) const override;

  GEOS_HOST_DEVICE
void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const;

  GEOS_HOST_DEVICE
void getElasticStiffness( localIndex const k,
                                    localIndex const q,
                                    real64 ( &stiffness )[6][6] ) const override;

  GEOS_HOST_DEVICE
void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( &elasticStrain )[6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void viscousStateUpdate( localIndex const k,
                                   localIndex const q,
                                   real64 beta ) const override;

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6],
                            real64 ( & stiffness )[6][6] ) const override final;

private:
  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;

  /// A reference to the ArrayView holding the strength scale.
  arrayView1d< real64 > const m_strengthScale;

  /// The critical length
  real64 const m_criticalLength;

  /// The maximum theoretical strength
  real64 const m_failureStrength;

  // The crack speed
  real64 const m_energyReleaseRate;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( totalStrain );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "smallStrainNoStateUpdate_StressOnly overload not implemented for Chiumenti" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &totalStrain )[6],
                                                 real64 ( &stress )[6],
                                                 real64 ( &stiffness )[6][6] ) const 
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( totalStrain );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "smallStrainNoStateUpdate overload not implemented for Chiumenti" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &totalStrain )[6],
                                                 real64 ( &stress )[6],
                                                 DiscretizationOps & stiffness ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( totalStrain );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "smallStrainNoStateUpdate overload not implemented for Chiumenti" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & timeIncrement,
                                                     real64 const ( &strainIncrement )[6],
                                                     real64 ( &stress )[6] ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_ERROR( "smallStrainUpdate_StressOnly overload not implemented for Chiumenti" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & timeIncrement,
                                                     real64 const ( & beginningRotation )[3][3],
                                                     real64 const ( & endRotation )[3][3],
                                                     real64 const ( &strainIncrement )[6],
                                                     real64 ( &stress )[6] ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( beginningRotation );
    GEOS_UNUSED_VAR( endRotation );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_ERROR( "smallStrainUpdate_StressOnly overload not implemented for Chiumenti" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::smallStrainUpdate( localIndex const k,
                                          localIndex const q,
                                          real64 const & timeIncrement,
                                          real64 const ( &strainIncrement )[6],
                                          real64 ( &stress )[6],
                                          real64 ( &stiffness )[6][6] ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "smallStrainUpdate overload not implemented for Chiumenti" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::smallStrainUpdate( localIndex const k,
                                          localIndex const q,
                                          real64 const & timeIncrement,
                                          real64 const ( &strainIncrement )[6],
                                          real64 ( &stress )[6],
                                          DiscretizationOps & stiffness ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "smallStrainUpdate overload not implemented for Chiumenti" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::getElasticStiffness( localIndex const k,
                                            localIndex const q,
                                            real64 ( &stiffness )[6][6] ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( stiffness );   
    GEOS_ERROR( "getElasticStiffness overload not implemented for Chiumenti" ); 
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::getElasticStrain( localIndex const k,
                                         localIndex const q,
                                         real64 ( &elasticStrain )[6] ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( elasticStrain );   
    GEOS_ERROR( "getElasticStrain overload not implemented for Chiumenti" ); 
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::viscousStateUpdate( localIndex const k,
                                           localIndex const q,
                                           real64 beta ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( beta );   
    GEOS_ERROR( "viscousStateUpdate overload not implemented for Chiumenti" ); 
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::hyperUpdate( localIndex const k,
                                    localIndex const q,
                                    real64 const ( & FminusI )[3][3],
                                    real64 ( & stress )[6] ) const
{
    // Compute hyperelastic trial update of stress
    HyperelasticMMSUpdates::hyperUpdate( k, 
                                         q,
                                         FminusI,
                                         stress );

    if( m_disableInelasticity )
    {
        return;
    }

  real64 failureStrength = m_failureStrength * m_strengthScale[k];
//   failureStrength *= 1.0 - m_damage[k][q];

  real64 youngModulus = conversions::lameConstants::toYoungMod( m_lambda[k], m_shearModulus[k] );
  real64 criticalLength = m_criticalLength * m_lengthScale[k];

  real64 brittlenessFactor = failureStrength * failureStrength / (2.0 * youngModulus * m_energyReleaseRate);
  real64 brittlenessFactorScaled = brittlenessFactor * criticalLength / (1.0 - brittlenessFactor * criticalLength);

  real64 principalStresses[3] = { 0 };
  real64 eigenVectors[3][3] = { { 0 } };
  LvArray::tensorOps::symEigenvectors< 3 >( principalStresses, eigenVectors, stress );

//   std::stringstream ss;
//   ss << "Stress in" << "\n";
//   for(int i = 0; i < 6;  i ++)
//   {
//     ss << stress[i] << ", ";
//   }
//   ss << "\n";

//   ss << "Eigen values: " << principalStresses[0] << ", " << principalStresses[1] << ", " << principalStresses[2] << "\n";
//   ss << "Eigen vectors" << "\n";
//   for(int i = 0; i < 3;  i ++)
//   {
//     for(int j = 0; j < 3; j ++)
//     {
//         ss << eigenVectors[i][j] << ", "; 
//     }
//     ss << "\n";
//   }
//   ss << "\n";

  // Find the largest principalStress
  real64 maximumPrincipalStress = 0.0;
  for( localIndex i = 0; i < 3; ++i )
  {
    maximumPrincipalStress = LvArray::math::max( principalStresses[i], maximumPrincipalStress );
  }

  if( maximumPrincipalStress > failureStrength )
  {
    real64 newDamage = 1.0;
    if( maximumPrincipalStress <= failureStrength*(1.0+1.0/brittlenessFactorScaled) )
    {
        newDamage = (1.0+brittlenessFactorScaled) * (1.0-failureStrength / maximumPrincipalStress);
    }   
    m_damage[k][q] = LvArray::math::max(0.0, LvArray::math::min(1.0, LvArray::math::max(m_damage[k][q], newDamage)));
  }
//   ss << "Updated damage" << m_damage[k][q] << "\n";

  // If damage is > 0.0, scale principal stresses
  if( m_damage[k][q] > 0.0 )
  {
    for( int i = 0; i < 3; i++ )
    {
        principalStresses[i] *= ( principalStresses[i] > 0.0 ) ? (1.0-m_damage[k][q]) : 1.0;
    }
  }

//   ss << "Dense stress" << "\n";
  real64 denseStress[3][3] = { { 0.0 } };
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
    { 
        denseStress[i][j] = principalStresses[0] * eigenVectors[0][i] * eigenVectors[0][j] + principalStresses[1] * eigenVectors[1][i] * eigenVectors[1][j] + principalStresses[2] * eigenVectors[2][i] * eigenVectors[2][j]; // Change magnitude and orientation of eigenvectors
        // ss << denseStress[i][j] << ", ";
    }
    // ss << "\n";
  }
//   ss << "\n";

  LvArray::tensorOps::denseToSymmetric< 3 >( stress, denseStress );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ChiumentiUpdates::hyperUpdate( localIndex const k,
                                          localIndex const q,
                                          real64 const ( & FminusI )[3][3],
                                          real64 ( & stress )[6],
                                          real64 ( & stiffness )[6][6] ) const
{
  hyperUpdate(k, q, FminusI, stress);
  getElasticStiffness( k, q, stiffness );
}

/**
 * @class Chiumenti
 *
 * Ceramic damage material model.
 */
class Chiumenti : public HyperelasticMMS
{
public:

  /// @typedef Alias for ChiumentiUpdates
  using KernelWrapper = ChiumentiUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Chiumenti( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Chiumenti() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Chiumenti";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }
    
    /// string/key for strength scale value
    static constexpr char const * strengthScaleString() { return "strengthScale"; }

    /// string/key for failure strength
    static constexpr char const * criticalLengthString() { return "criticalLength"; }

    /// string/key for failure strength
    static constexpr char const * failureStrengthString() { return "failureStrength"; }

    /// string/key for energy release rate
    static constexpr char const * energyReleaseRateString() { return "energyReleaseRate"; }
  };

  /**
   * @brief Create a instantiation of the ChiumentiUpdate class that refers to the data in this.
   * @return An instantiation of ChiumentiUpdate.
   */
  ChiumentiUpdates createKernelUpdates() const
  {
    return ChiumentiUpdates( m_damage,
                             m_lengthScale,
                             m_strengthScale,
                             m_criticalLength,
                             m_failureStrength,
                             m_energyReleaseRate,
                             m_lambda,
                             m_shearModulus,
                             m_thermalExpansionCoefficient,
                             m_newStress,
                             m_oldStress,
                             m_density,
                             m_wavespeed,
                             m_disableInelasticity );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_damage,
                          m_lengthScale,
                          m_strengthScale,
                          m_criticalLength,
                          m_failureStrength,
                          m_energyReleaseRate,
                          m_lambda,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Material parameter: The strength scale values
  array1d< real64 > m_strengthScale;

  /// Material parameter: The value of critical length
  real64 m_criticalLength;

  /// Material parameter: The value of failure strength
  real64 m_failureStrength;

  /// Material parameter: The value of energyReleaseRate
  real64 m_energyReleaseRate;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_CHIUMENTI_HPP_ */
