/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Graphite.hpp
 * @brief Simple damage model for modeling material failure in brittle materials.
 *
 * This damage model is intended for use with damage-field partitioning (DFG) within the
 * MPM solver, but can also be used without DFG by any solver. It is only appropriate for
 * schemes implementing explicit time integration. The model is really a hybrid plasticity/
 * damage model in the sense that we assume damaged material behaves like granular material
 * and hence follows a modified Mohr-Coulomb law. The modifications are that at low pressures,
 * the shape of the yield surface is modified to resemble a maximum principal stress criterion,
 * and at high pressures the shape converges on the von Mises yield surface. The damage
 * parameter results in softening of the deviatoric response i.e. causes the yield surface to
 * contract. Furthermore, damage is used to scale back tensile pressure: p = (1 - d) * pTrial.
 * pTrial is calculatd as pTrial = -k * log(J), where the Jacobian J of the material motion is
 * integrated and tracked by this model.
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_GRAPHITE_HPP_
#define GEOS_CONSTITUTIVE_SOLID_GRAPHITE_HPP_

#include "SolidBase.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class GraphiteUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class GraphiteUpdates : public SolidBaseUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] velocityGradient The ArrayView holding the velocity gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point.
   * @param[in] relaxation The ArrayView holding the relaxation for each quadrature point.
   * @param[in] damage The ArrayView holding the damage for each quardrature point.
   * @param[in] jacobian The ArrayView holding the jacobian for each quardrature point.
   * @param[in] materialDirection The ArrayView holding the material direction for each element/particle.
   * @param[in] lengthScale The ArrayView holding the length scale for each element.
   * @param[in] strengthScale The ArrayView holding the strength scale for each element.
   * @param[in] failureStrength The failure strength.
   * @param[in] crackSpeed The crack speed velocity.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  GraphiteUpdates( real64 const & defaultYoungModulusTransverse,
                   real64 const & defaultYoungModulusAxial,
                   real64 const & defaultPoissonRatioTransverse,
                   real64 const & defaultPoissonRatioAxialTransverse,
                   real64 const & defaultShearModulusAxialTransverse,
                   arrayView1d< real64 > const & effectiveBulkModulus,
                   arrayView1d< real64 > const & effectiveShearModulus,
                   arrayView2d< real64 > const & materialDirection,
                   real64 const defaultYoungModulusTransversePressureDerivative,
                   real64 const defaultYoungModulusAxialPressureDerivative,
                   real64 const defaultShearModulusAxialTransversePressureDerivative,
                   arrayView3d< real64 > const & velocityGradient,
                   arrayView3d< real64 > const & plasticStrain,
                   arrayView2d< real64 > const & relaxation,
                   arrayView2d< real64 > const & damage,
                   arrayView2d< real64 > const & jacobian,
                   arrayView1d< real64 > const & lengthScale,
                   arrayView1d< real64 > const & strengthScale,
                   real64 const & failureStrength,
                   real64 const & crackSpeed,
                   real64 const & damagedMaterialFrictionalSlope,
                   real64 const & distortionShearResponseX2,
                   real64 const & distortionShearResponseY1,
                   real64 const & distortionShearResponseY2,
                   real64 const & distortionShearResponseM1,
                   real64 const & inPlaneShearResponseX2,
                   real64 const & inPlaneShearResponseY1,
                   real64 const & inPlaneShearResponseY2,
                   real64 const & inPlaneShearResponseM1,
                   real64 const & coupledShearResponseX2,
                   real64 const & coupledShearResponseY1,
                   real64 const & coupledShearResponseY2,
                   real64 const & coupledShearResponseM1,
                   real64 const & maximumPlasticStrain,
                   arrayView1d< real64 > const & thermalExpansionCoefficient,
                   arrayView3d< real64, solid::STRESS_USD > const & newStress,
                   arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                   arrayView2d< real64 > const & density,
                   arrayView2d< real64 > const & wavespeed,
                   bool const & disableInelasticity ):
     SolidBaseUpdates( newStress, 
                       oldStress,
                       density,
                       wavespeed,
                       thermalExpansionCoefficient, 
                       disableInelasticity ),
    m_defaultYoungModulusTransverse( defaultYoungModulusTransverse ),
    m_defaultYoungModulusAxial( defaultYoungModulusAxial ),
    m_defaultPoissonRatioTransverse( defaultPoissonRatioTransverse ),
    m_defaultPoissonRatioAxialTransverse( defaultPoissonRatioAxialTransverse ),
    m_defaultShearModulusAxialTransverse( defaultShearModulusAxialTransverse ),
    m_effectiveBulkModulus( effectiveBulkModulus ),
    m_effectiveShearModulus( effectiveShearModulus ),
    m_materialDirection( materialDirection ),
    m_defaultYoungModulusTransversePressureDerivative( defaultYoungModulusTransversePressureDerivative ),
    m_defaultYoungModulusAxialPressureDerivative( defaultYoungModulusAxialPressureDerivative ),
    m_defaultShearModulusAxialTransversePressureDerivative( defaultShearModulusAxialTransversePressureDerivative ),                   
    m_velocityGradient( velocityGradient ),
    m_plasticStrain( plasticStrain ),
    m_relaxation( relaxation ),
    m_damage( damage ),
    m_jacobian( jacobian ),
    m_lengthScale( lengthScale ),
    m_strengthScale( strengthScale ),
    m_failureStrength( failureStrength ),
    m_crackSpeed( crackSpeed ),
    m_damagedMaterialFrictionalSlope( damagedMaterialFrictionalSlope ),
    m_distortionShearResponseX2( distortionShearResponseX2 ),
    m_distortionShearResponseY1( distortionShearResponseY1 ),
    m_distortionShearResponseY2( distortionShearResponseY2 ),
    m_distortionShearResponseM1( distortionShearResponseM1 ),
    m_inPlaneShearResponseX2( inPlaneShearResponseX2 ),
    m_inPlaneShearResponseY1( inPlaneShearResponseY1 ),
    m_inPlaneShearResponseY2( inPlaneShearResponseY2 ),
    m_inPlaneShearResponseM1( inPlaneShearResponseM1 ),
    m_coupledShearResponseX2( coupledShearResponseX2 ),
    m_coupledShearResponseY1( coupledShearResponseY1 ),
    m_coupledShearResponseY2( coupledShearResponseY2 ),
    m_coupledShearResponseM1( coupledShearResponseM1 ),
    m_maximumPlasticStrain( maximumPlasticStrain )
  {}

  /// Default copy constructor
  GraphiteUpdates( GraphiteUpdates const & ) = default;

  /// Default move constructor
  GraphiteUpdates( GraphiteUpdates && ) = default;

  /// Deleted default constructor
  GraphiteUpdates() = delete;

  /// Deleted copy assignment operator
  GraphiteUpdates & operator=( GraphiteUpdates const & ) = delete;

  /// Deleted move assignment operator
  GraphiteUpdates & operator=( GraphiteUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using SolidBaseUpdates::smallStrainUpdate;

  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 ( &stiffness )[6][6] ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
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
  void smallStrainUpdateHelper( localIndex const k,
                                localIndex const q,
                                real64 const timeIncrement,
                                real64 const ( & beginningRotation )[3][3],
                                real64 const ( & endRotation )[3][3],
                                real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  void computeTransverselyIsotropicTrialStress( const real64 timeIncrement,      
                                                const real64 Ez,                 
                                                const real64 Ep,                 
                                                const real64 nuzp,               
                                                const real64 nup,                
                                                const real64 Gzp,                
                                                real64 const (& materialDirection)[3], 
                                                real64 const (& oldStress)[6],
                                                real64 const (& D)[6],           
                                                real64 (& newStress) [6] ) const;

  GEOS_HOST_DEVICE
  void computeTransverselyIsotropicPlasticStrainIncrement( real64 const ( & velocityGradient )[3][3],
                                                           real64 const ( & oldStress )[6],
                                                           real64 const ( & newStress )[6],
                                                           const real64 Ez,             
                                                           const real64 Ep,             
                                                           const real64 nuzp,           
                                                           const real64 nup,            
                                                           const real64 Gzp,            
                                                           real64 const ( & materialDirection )[3],	
                                                           const real64 timeIncrement,
                                                           real64 ( & plasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB1( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB2( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB3( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB4( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB5( real64 const (& materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;

  GEOS_HOST_DEVICE
  real64 delta( const int i,
                const int j)  const;

  GEOS_HOST_DEVICE
  real64 slopePoint0( const real64 x,
                      const real64 x1,
                      const real64 x2,
                      const real64 y1,
                      const real64 y2,
                      const real64 m1 ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    SolidBaseUpdates::saveConvergedState( k, q );
  }

private:
  real64 const & m_defaultYoungModulusTransverse;
  
  real64 const & m_defaultYoungModulusAxial;
  
  real64 const & m_defaultPoissonRatioTransverse;

  real64 const & m_defaultPoissonRatioAxialTransverse;
  
  real64 const & m_defaultShearModulusAxialTransverse;

  arrayView1d< real64 > const m_effectiveBulkModulus;
  
  arrayView1d< real64 > const m_effectiveShearModulus;
  
  arrayView2d< real64 > const m_materialDirection;
  
  real64 const m_defaultYoungModulusTransversePressureDerivative;
  
  real64 const m_defaultYoungModulusAxialPressureDerivative; 
  
  real64 const m_defaultShearModulusAxialTransversePressureDerivative;

  /// A reference to the ArrayView holding the velocity gradient for each element/particle.
  arrayView3d< real64 > const m_velocityGradient;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_relaxation;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;

  /// Discretization-sized variable: The strength scale for each element/particle
  arrayView1d< real64 > const m_strengthScale;

  /// The maximum theoretical strength
  real64 const m_failureStrength;

  // The time to failure
  real64 const m_crackSpeed;

  // The damaged material frictional slope
  real64 const m_damagedMaterialFrictionalSlope;

  /// Material parameters: The values controlling the failure envelope
  real64 const m_distortionShearResponseX2;
  real64 const m_distortionShearResponseY1;
  real64 const m_distortionShearResponseY2;
  real64 const m_distortionShearResponseM1;

  real64 const m_inPlaneShearResponseX2;
  real64 const m_inPlaneShearResponseY1;
  real64 const m_inPlaneShearResponseY2;
  real64 const m_inPlaneShearResponseM1;

  real64 const m_coupledShearResponseX2;
  real64 const m_coupledShearResponseY1;
  real64 const m_coupledShearResponseY2;
  real64 const m_coupledShearResponseM1;
   
  real64 const m_maximumPlasticStrain;

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const & timeIncrement,
                                         real64 const ( &strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for Graphite model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const & timeIncrement,
                                         real64 const ( &strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, 
                     q, 
                     timeIncrement, 
                     strainIncrement, 
                     stress, 
                     stiffness.m_c );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                         localIndex const q,
                                                         real64 const & timeIncrement,
                                                         real64 const ( & strainIncrement )[6],
                                                         real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "smallStrainUpdateStressOnly overload not implemented for CeramicDamage model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const & timeIncrement,
                                                    real64 const ( & beginningRotation )[3][3],
                                                    real64 const ( & endRotation )[3][3],
                                                    real64 const ( & strainIncrement )[6],
                                                    real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  // ElasticTransverseIsotropicPressureDependentUpdates::smallStrainUpdate_StressOnly( k,
  //                                                                                   q,
  //                                                                                   timeIncrement,
  //                                                                                   beginningRotation,
  //                                                                                   endRotation,
  //                                                                                   strainIncrement,
  //                                                                                   stress );

  m_jacobian[k][q] *= LvArray::math::exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  // if( m_disableInelasticity )
  // {
  //   return;
  // }

  // call the constitutive model
  GraphiteUpdates::smallStrainUpdateHelper( k, 
                                            q, 
                                            timeIncrement,                
                                            beginningRotation,
                                            endRotation,
                                            stress );

  // save new stress and return
  saveStress( k, q, stress );

  // CC: debug
  // GEOS_LOG_RANK( "Particle " << k << ", Saved stress: {" << stress[0] << ", " << stress[1] << ", " << stress[2] << ", " << stress[3] << ", " << stress[4] << ", " << stress[5] << "}" );
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdateHelper( localIndex const k,
                                               localIndex const q,
                                               real64 const timeIncrement,  
                                               real64 const ( & beginningRotation )[3][3],
                                               real64 const ( & endRotation )[3][3],
                                               real64 ( & stress )[6] ) const
{   
    // CC: debug
    // GEOS_LOG_RANK( "Particle " << k << ", Stress in: {" << stress[0] << ", " << stress[1] << ", " << stress[2] << ", " << stress[3] << ", " << stress[4] << ", " << stress[5] << "}" );
    real64 oldStress[6] = { 0 };
    LvArray::tensorOps::copy< 6 >( oldStress, m_oldStress[k][q]); //stress );
    // CC: debug
    // GEOS_LOG_RANK( "Particle " << k << ", Old stress copy: {" << oldStress[0] << ", " << oldStress[1] << ", " << oldStress[2] << ", " << oldStress[3] << ", " << oldStress[4] << ", " << oldStress[5] << "}" );
    
    real64 rotationTranspose[3][3] = { { 0 } };
    LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation ); 

    real64 tempMat[ 3 ][ 3 ]= { { 0 } };
    
    real64 unrotatedVelocityGradient[3][3] = { { 0 } };
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( tempMat, beginningRotation, m_velocityGradient[k] );
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedVelocityGradient, tempMat, beginningRotation );

    real64 unrotatedVelocityGradientTranspose[3][3] = { { 0 } };
    LvArray::tensorOps::transpose< 3, 3 >( unrotatedVelocityGradientTranspose, unrotatedVelocityGradient );

    // CC: Is there an LvArray operation to get the symmetric part of a matrix?
    real64 denseD[3][3] = { { 0 } };
    LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
    LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
    LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );

    real64 D[6] = { 0 };
    LvArray::tensorOps::denseToSymmetric<3>( D, denseD );

    // make sure material direction is normalized.
    real64 materialDirection[3] = { 0 };
    LvArray::tensorOps::copy< 3 >( materialDirection, m_materialDirection[k] );
    LvArray::tensorOps::normalize< 3 >( materialDirection );

    // Unrotate material direction
    real64 unrotatedMaterialDirection[3] = { 0 };
    LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( unrotatedMaterialDirection, rotationTranspose, materialDirection);

    // Use beginning of step normal stress to compute stress dependence of Ez
    real64 temp[3] = { 0 };
    int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
    LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, oldStress, unrotatedMaterialDirection );
    real64 oldPlaneNormalStress = LvArray::tensorOps::AiBi< 3 >( unrotatedMaterialDirection, temp );  // CC: Unused?

    // Beginning of step pressure to compute pressure-dependence of elastic moduli
    real64 oldPressure = (-1.0/3.0)*( oldStress[0] + oldStress[1] + oldStress[2] ); // CC: Unused?

    // This is a transversely isotropic material for graphite-like crystals having
    // some weak plane with plane-normal-stress- and pressure-dependent elastic properties:

    // CC: in old geos the elastic on two different lines of the same code, Mike used planeNormalStress in one but not the other
    // need to ask him about that
    real64 Ez = m_defaultYoungModulusAxial + m_defaultYoungModulusAxialPressureDerivative * LvArray::math::max( 0.0, -0.5*oldPlaneNormalStress + 0.5*oldPressure);
    real64 Ep = m_defaultYoungModulusTransverse + m_defaultYoungModulusTransversePressureDerivative * LvArray::math::max( 0.0, oldPressure);
    real64 Gzp  = m_defaultShearModulusAxialTransverse + m_defaultShearModulusAxialTransversePressureDerivative * LvArray::math::max( 0.0, oldPressure );
    real64 nuzp = m_defaultPoissonRatioAxialTransverse;
    real64 nup = m_defaultPoissonRatioTransverse;

    // Update effective elastic properties
    m_effectiveBulkModulus[k] = -Ep*Ez/(2*Ez*(nup+nuzp-1) + Ep*(2*nuzp-1));
    m_effectiveShearModulus[k] = 0.6*m_effectiveBulkModulus[k];
    m_wavespeed[k][0] = LvArray::math::sqrt( ( m_effectiveBulkModulus[k] + (4.0/3.0) * m_effectiveShearModulus[k] ) / m_density[k][0] );

    // CC: debug
    // GEOS_LOG_RANK( "Particle " << k << ":\n"
    //                 "Old pressure: " << oldPressure << "\n"
    //                 "\tEz=" << Ez << "\n" << 
    //                 "\tEp=" << Ep << "\n" << 
    //                 "\tGzp=" << Gzp << "\n" << 
    //                 "\tnuzp=" << nuzp << "\n" << 
    //                 "\tnup=" << nup << "\n" << 
    //                 "\teffK=" << m_effectiveBulkModulus[k] << "\n" << 
    //                 "\teffG=" << m_effectiveShearModulus[k] << "\n");

    // real64 Ez = c33 - c13 * c13 / ( c11 - c66 );
    // real64 Ep = 4 * c66 * ( c11 * c33 - c66 * c33 - c13 * c13 ) / ( c11 * c33 - c13 * c13 );
    // real64 Gzp = c44 / 2.0;
    // real64 nuzp = c13 / ( 2 * ( c11 - c66 ) );
    // real64 nup = Ep / c66 - 1; //4 * ( c11 * c33 - c66 * c33 - c13 * c13 ) / ( c11 * c33 - c13 * c13 ) - 1;

    // Hypoelastic trial stress update.
    computeTransverselyIsotropicTrialStress( timeIncrement,      // time step
                                             Ez,                 // Elastic modulus preferred direction
                                             Ep,                 // Elastic modulus transverse plane
                                             nuzp,               // Poisson ratio coupled
                                             nup,                // Poisson ratio transverse plane
                                             Gzp,                // Shear modulus coupled plane
                                             unrotatedMaterialDirection,	 // preferred direction
                                             oldStress,          // stress at start of step
                                             D,                  // D=sym(L)
                                             stress );           // stress at end of step

    // CC: debug
    // GEOS_LOG_RANK( "Particle " << k << ", Trial stress: {" << stress[0] << ", " << stress[1] << ", " << stress[2] << ", " << stress[3] << ", " << stress[4] << ", " << stress[5] << "}" );

    // m_jacobian[k][q] *= LvArray::math::exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

    if( m_disableInelasticity )
    {
      return;
    }

    // Decompose stress tensor into pieces.
    real64 sigma1Dense [3][3] = { { 0 } };
    real64 sigma2Dense [3][3] = { { 0 } };
    real64 sigma4Dense [3][3] = { { 0 } };
    real64 sigma5Dense [3][3] = { { 0 } };
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            for(int p=0; p<3; p++)
            {
                for(int w=0; w<3; w++)
                {
                    sigma1Dense[i][j] += transverselyIsotropicB1(unrotatedMaterialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                    sigma2Dense[i][j] += transverselyIsotropicB2(unrotatedMaterialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                    sigma4Dense[i][j] += transverselyIsotropicB4(unrotatedMaterialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                    sigma5Dense[i][j] += transverselyIsotropicB5(unrotatedMaterialDirection,i,j,p,w)*stress[voigtMap[p][w]];
                }
            }
        }
    }

    real64 sigma1[6] = {0}; // axial
    real64 sigma2[6] = {0}; // in-plane normal
    real64 sigma4[6] = {0}; // in-plane total stress
    real64 sigma5[6] = {0}; // weak plane - shear
    LvArray::tensorOps::denseToSymmetric< 3 >( sigma1, sigma1Dense );
    LvArray::tensorOps::denseToSymmetric< 3 >( sigma2, sigma2Dense );
    LvArray::tensorOps::denseToSymmetric< 3 >( sigma4, sigma4Dense );
    LvArray::tensorOps::denseToSymmetric< 3 >( sigma5, sigma5Dense );

    // Trial pressure to compute pressure-dependence of strength
    real64 pressure = (-1.0/3.0)*( stress[0] + stress[1] + stress[2] );

    // CC: debug
    // GEOS_LOG_RANK( "Particle " << k << ":\n" <<
    //                 "\toldStress: {" << oldStress[0] << ", " << oldStress[1] << ", " << oldStress[2] << ", " << oldStress[3] << ", " << oldStress[4] << ", " << oldStress[5] << "}" << 
    //                 "\tnewStress: {" << stress[0] << ", " << stress[1] << ", " << stress[2] << ", " << stress[3] << ", " << stress[4] << ", " << stress[5] << "}" << 
    //                 "\tpressure: " << pressure);


    // Check for tensile failure in preferred direction
    // real64 temp[3] = { 0 };
    LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, stress, unrotatedMaterialDirection );
    real64 planeNormalStress = LvArray::tensorOps::AiBi< 3 >( unrotatedMaterialDirection, temp );

    real64 failureStrength = m_failureStrength * m_strengthScale[k];

    // increment damage, but enforce 0<=d<=1
    if ( planeNormalStress > failureStrength )
    {
        real64 timeToFailure = m_lengthScale[k] / m_crackSpeed;
        m_damage[k][q] = LvArray::math::min( m_damage[k][q] + timeIncrement / timeToFailure, 1.0 );
    }

    // strength scale factor, combining plastic softening and damage
    real64 fac = (1.0 - m_damage[k][q])*(1.0 - m_relaxation[k][q]);

    // CC: debug
    // GEOS_LOG_RANK( "Particle " << k << ": dmg: " << m_damage[k][q] << ", Fac: " << fac );

    // Enforce damage, no tensile stress on plane, and frictional response to shear.
    if(m_damage[k][q] >= 1.0) //Some compilers complain about == comparison with floating point numbers (use tolerance check?)
    {
        if ( planeNormalStress > 0 )
        {
            LvArray::tensorOps::scale< 6 >(sigma1, 0.0);
        }
    }

    // find in-plane isotropic and deviatoric stress
    real64 inPlaneIso[6];
    LvArray::tensorOps::copy< 6 >( inPlaneIso, sigma2 );
    LvArray::tensorOps::scale< 6 >( inPlaneIso, 0.5 );

    real64 inPlaneDev[6];
    LvArray::tensorOps::copy< 6 >( inPlaneDev, sigma4 );
    LvArray::tensorOps::subtract< 6 >( inPlaneDev, inPlaneIso );

    // find distortion stress from in-plane isotropic and plane-normal stress
    real64 distortion[6] = { 0 };
    LvArray::tensorOps::copy< 6 >(distortion, inPlaneIso);
    LvArray::tensorOps::add< 6 >(distortion, sigma1);

    real64 distortionMeanStress = ( distortion[0] + distortion[1] + distortion[2] ) / 3.0;

    real64 distortion_iso[6] = {0};
    distortion_iso[0] = distortionMeanStress;
    distortion_iso[1] = distortionMeanStress;
    distortion_iso[2] = distortionMeanStress;

    real64 distortion_dev[6] = { 0 };
    LvArray::tensorOps::copy< 6 >( distortion_dev, distortion );
    LvArray::tensorOps::subtract< 6 >( distortion_dev, distortion_iso );

    // Define 3 pressure-dependent shear strengths for different shear modes.
    real64 inPlaneShearStrength, coupledYieldStrength, totalShearStrength;
    real64 x1, x2, y1, y2, m1;

    // -------------------------------
    // "distortion" shear relating sigma normal and in-plane iso"
    x1 = 0;
    x2 = m_distortionShearResponseX2;
    y1 = m_distortionShearResponseY1 * m_strengthScale[k];
    y2 = m_distortionShearResponseY2 * m_strengthScale[k];
    m1 = m_distortionShearResponseM1;

    // damage or softening reduces cohesion and reduces slope to failed value.
    y1 *= fac;
    m1 = fac*m1 + (1.0 - fac)*LvArray::math::max( m_damagedMaterialFrictionalSlope, (y2-y1)/(x2-x1) );
    if(pressure < x1)
    {
        totalShearStrength=LvArray::math::max(0.0,y1-(x1-pressure)*m1);
        // totalShearStrength=LvArray::math::max(0.0,y1-(x2-pressure)*m1);
    }
    else if(pressure < x2)
    {
        totalShearStrength=slopePoint0(pressure,x1,x2,y1,y2,m1);
    }
    else
    {
        totalShearStrength=y2;
    }

    // -------------------------------------------
    // Coupled shear response (slip on weak plane)
    x1 = 0;
    x2 = m_coupledShearResponseX2;
    y1 = m_coupledShearResponseY1 * m_strengthScale[k];
    y2 = m_coupledShearResponseY2 * m_strengthScale[k];
    m1 = m_coupledShearResponseM1;

    // damage or softening reduces cohesion and reduces slope to failed value.
    y1 *= fac;
    m1 = fac*m1 + (1. - fac)*LvArray::math::max( m_damagedMaterialFrictionalSlope, (y2-y1)/(x2-x1) );
    if(pressure<x1)
    {
      coupledYieldStrength=LvArray::math::max(0.0,y1-(x1-pressure)*m1);
      // coupledYieldStrength=LvArray::math::max(0.0,y1-(x2-pressure)*m1);
    }
    else if(pressure<x2)
    {
        coupledYieldStrength=slopePoint0(pressure,x1,x2,y1,y2,m1);
    }
    else
    {
        coupledYieldStrength=y2;
    }

    // -------------------------------
    // In-plane shear response
    x1 = 0;
    x2 = m_inPlaneShearResponseX2; 
    y1 = m_inPlaneShearResponseY1 * m_strengthScale[k];
    y2 = m_inPlaneShearResponseY2 * m_strengthScale[k];
    m1 = m_inPlaneShearResponseM1;

    // damage or softening reduces cohesion and reduces slope to failed value.
    y1 *= fac;
    m1 = fac*m1 + (1. - fac)*LvArray::math::max( m_damagedMaterialFrictionalSlope, (y2-y1)/(x2-x1) );
    if( pressure < x1 )
    {
      inPlaneShearStrength=LvArray::math::max(0.0, y1  -( x1 - pressure ) * m1 );
      // inPlaneShearStrength=LvArray::math::max(0.0, y1  -( x2 - pressure ) * m1 );
    }
    else if( pressure < x2 )
    {
        inPlaneShearStrength=slopePoint0( pressure, x1, x2, y1, y2, m1);
    }
    else
    {
        inPlaneShearStrength = y2;
    }

    // Enforce strength ---------------------------------------------

    // CC: add check to make sure strength is some small value to avoid floating point error for dmg = 1.0 initially in strength computations

    // flags indicating whether plastic strain needs to be updated.
    bool plastic = false;

    // Enforce distortion strain yield
    //For symmetric matrix need to double off diagonal elements for l2Norm?
    distortion_dev[3] *= 1.41421356237;
    distortion_dev[4] *= 1.41421356237;
    distortion_dev[5] *= 1.41421356237;
    real64 totalShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( distortion_dev );
    distortion_dev[3] /= 1.41421356237;
    distortion_dev[4] /= 1.41421356237;
    distortion_dev[5] /= 1.41421356237;
    if ( totalShearStress > totalShearStrength )
    {
      LvArray::tensorOps::scale< 6 >( distortion_dev, totalShearStrength / totalShearStress );
      plastic = true;
    }

    // enforce in-plane yield
    //For symmetric matrix need to double off diagonal elements for l2Norm?
    inPlaneDev[3] *= 1.41421356237;
    inPlaneDev[4] *= 1.41421356237;
    inPlaneDev[5] *= 1.41421356237;
    real64 inPlaneShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( inPlaneDev ) ;
    inPlaneDev[3] /= 1.41421356237;
    inPlaneDev[4] /= 1.41421356237;
    inPlaneDev[5] /= 1.41421356237;
    if ( inPlaneShearStress > inPlaneShearStrength )
    {
      LvArray::tensorOps::scale< 6 >( inPlaneDev, inPlaneShearStrength / inPlaneShearStress );
      plastic = true;
    }

    // enforce coupled yield
    //For symmetric matrix need to double off diagonal elements for l2Norm?
    sigma5[3] *= 1.41421356237;
    sigma5[4] *= 1.41421356237;
    sigma5[5] *= 1.41421356237;
    real64 coupledShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( sigma5 );
    sigma5[3] /= 1.41421356237;
    sigma5[4] /= 1.41421356237;
    sigma5[5] /= 1.41421356237;
    if ( coupledShearStress > coupledYieldStrength )
    {
      LvArray::tensorOps::scale< 6 >( sigma5, coupledYieldStrength / coupledShearStress );
      plastic = true;
    }

    // reassemble end-of-step stress  
    real64 newStress[6] = {0};
    LvArray::tensorOps::add< 6 >(newStress, distortion_iso);
    LvArray::tensorOps::add< 6 >(newStress, distortion_dev);
    LvArray::tensorOps::add< 6 >(newStress, inPlaneDev);
    LvArray::tensorOps::add< 6 >(newStress, sigma5);

    // Copy the new stress 
    LvArray::tensorOps::copy< 6 >( stress, newStress );

    //////////////////////////////////////////////////////
    // Evolve state variables.
    //
    // For inelastic response evolve damage if pressure is below the
    // brittle ductile transition pressure and evolve plastic strain
    // if pressure is above brittle ductile transition pressure.
    // plastic strain doesn't currently affect material response, but
    // it may be useful for plotting regions of high-pressure yield.
    if( plastic )
    {
      // CC: debug
        // GEOS_LOG_RANK( "Particle " << k << ": PlaneNormalStress: " << planeNormalStress << ", " << 
        //                                      "TotalShearStress: " << totalShearStress << " (" << totalShearStrength << ")" << ", " << 
        //                                      "inPlaneShearStress: " << inPlaneShearStress << " (" << inPlaneShearStrength << ")" << ", " << 
        //                                      "coupledShearStress: " << coupledShearStress << " (" << coupledYieldStrength << ")");

        // increment plastic strain
        real64 plasticStrainIncrement[6] = {0};
        computeTransverselyIsotropicPlasticStrainIncrement( unrotatedVelocityGradient, // Velocity gradient tensor
                                                            oldStress,                 // Stress at start of step
                                                            stress,                    // Stress at end of step
                                                            Ez,                        // Elastic modulus preferred direction
                                                            Ep,                        // Elastic modulus transverse plane
                                                            nuzp,                      // Poisson ratio coupled
                                                            nup,                       // Poisson ratio transverse plane
                                                            Gzp,                       // Shear modulus coupled plane
                                                            unrotatedMaterialDirection,		     // material direction (unit vector, plane normal)
                                                            timeIncrement,             // timeStep
                                                            plasticStrainIncrement );

        // Unrotate old plastic strain
        real64 oldPlasticStrain[6] = { 0 };
        LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] );
        oldPlasticStrain[3] *= 0.5;
        oldPlasticStrain[4] *= 0.5;
        oldPlasticStrain[5] *= 0.5;

        real64 unrotatedOldPlasticStrain[6] = { 0 };
        LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( unrotatedOldPlasticStrain, rotationTranspose, oldPlasticStrain );
        unrotatedOldPlasticStrain[3] *= 2.0;
        unrotatedOldPlasticStrain[4] *= 2.0;
        unrotatedOldPlasticStrain[5] *= 2.0; 

        // GEOS_LOG_RANK( "Particle: " << k << "\n " <<
        //               "\tOld Stress: {" << oldStress[0] << ", " << oldStress[1] << ", "<< oldStress[2] << ", "<< oldStress[3] << ", "<< oldStress[4] << ", "<< oldStress[5] << "}\n " << 
        //               "\tNew Stress: {" << stress[0] << ", " << stress[1] << ", "<< stress[2] << ", "<< stress[3] << ", "<< stress[4] << ", "<< stress[5] << "}\n " <<
        //               "\tStress Incr: {" << stress[0] - oldStress[0] << ", " << stress[1] - oldStress[1] << ", "<< stress[2] - oldStress[2] << ", "<< stress[3] - oldStress[3] << ", "<< stress[4] - oldStress[4] << ", "<< stress[5] - oldStress[5] << "}\n " <<
        //               "\tOld plastic strain: {" << unrotatedOldPlasticStrain[0] << ", " << unrotatedOldPlasticStrain[1] << ", "<< unrotatedOldPlasticStrain[2] << ", "<< unrotatedOldPlasticStrain[3] << ", "<< unrotatedOldPlasticStrain[4]<< ", "<< unrotatedOldPlasticStrain[5] << "}"
        //               "\tPlastic strain increment: {" << plasticStrainIncrement[0] << ", " << plasticStrainIncrement[1] << ", " << plasticStrainIncrement[2] << ", " << plasticStrainIncrement[3] << ", "<< plasticStrainIncrement[4]<< ", " << plasticStrainIncrement[5] << "}" );

        // Add plastic strain increment
        real64 unrotatedNewPlasticStrain[6] = { 0 };
        LvArray::tensorOps::copy< 6 >(unrotatedNewPlasticStrain, unrotatedOldPlasticStrain);
        LvArray::tensorOps::add< 6 >(unrotatedNewPlasticStrain, plasticStrainIncrement);

        // Rotate new plastic strain
        unrotatedNewPlasticStrain[3] *= 0.5;
        unrotatedNewPlasticStrain[4] *= 0.5;
        unrotatedNewPlasticStrain[5] *= 0.5;
        
        real64 newPlasticStrain[ 6 ] = { 0 }; 
        LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( newPlasticStrain, endRotation, unrotatedNewPlasticStrain );
        newPlasticStrain[3] *= 2.0;
        newPlasticStrain[4] *= 2.0;
        newPlasticStrain[5] *= 2.0;

        // Assign new plastic strain to state variable
        LvArray::tensorOps::copy< 6 >(m_plasticStrain[k][q], newPlasticStrain);

        // increment relaxation
        //For symmetric matrix need to double off diagonal elements for l2Norm?
        plasticStrainIncrement[3] *= 1.41421356237;
        plasticStrainIncrement[4] *= 1.41421356237;
        plasticStrainIncrement[5] *= 1.41421356237;
        m_relaxation[k][q] += LvArray::tensorOps::l2Norm< 6 >( plasticStrainIncrement ) / m_maximumPlasticStrain;
        m_relaxation[k][q] = LvArray::math::min(1.0, m_relaxation[k][q]);
    }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::computeTransverselyIsotropicTrialStress(const real64 timeIncrement,      // time step
                                                              const real64 Ez,                 // Elastic modulus preferred direction
                                                              const real64 Ep,                 // Elastic modulus transverse plane
                                                              const real64 nuzp,               // Poisson ratio coupled
                                                              const real64 nup,                // Poisson ratio transverse plane
                                                              const real64 Gzp,                // Shear modulus coupled plane
                                                              real64 const (& materialDirection)[3], // preferred direction
                                                              real64 const (& oldStress)[6],   // stress at start of step.
                                                              real64 const (& D)[6],           // D=sym(L)
                                                              real64 (& newStress) [6]         // stress at end of step
) const
{
  // These are the TI elastic stiffness coefficients using Brannon's TI basis tensors:
	real64 h1 = ( Ez*Ez*(-1 + nup) ) / ( Ez*(-1 + nup) + 2*Ep*nuzp*nuzp );
	real64 h2 = -( ( Ep*(Ez*nup + Ep*nuzp*nuzp) ) / ( (1 + nup)*(Ez*(-1 + nup) + 2*Ep*nuzp*nuzp ) ) );
	real64 h3 = ( Ep*Ez*nuzp ) / ( Ez - Ez*nup - 2*Ep*nuzp*nuzp );
	real64 h4 = Ep/( 1 + nup );
	real64 h5 = 2*Gzp;

  real64 stressIncrementDense[3][3] = { { 0 } };
  int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int p=0; p<3; p++)
			{
				for(int w=0; w<3; w++)
				{
					// newStress[voigtMap[i][j]] 
          stressIncrementDense[i][j] += ( h1*transverselyIsotropicB1(materialDirection,i,j,p,w) +
                                          h2*transverselyIsotropicB2(materialDirection,i,j,p,w) +
                                          h3*transverselyIsotropicB3(materialDirection,i,j,p,w) +
                                          h4*transverselyIsotropicB4(materialDirection,i,j,p,w) +
                                          h5*transverselyIsotropicB5(materialDirection,i,j,p,w))*D[voigtMap[p][w]]*timeIncrement;
				}
			}
		}
	}

  real64 stressIncrement[6] = { 0 };
  LvArray::tensorOps::denseToSymmetric< 3 >( stressIncrement, stressIncrementDense );

  LvArray::tensorOps::copy< 6 >(newStress, oldStress);
  LvArray::tensorOps::add< 6 >( newStress, stressIncrement );

// CC: debug
  // GEOS_LOG_RANK( "h constants: {" << h1 << ", " << h2 << ", " << h3 << ", " << h4 << ", " << h5 << "}\n" <<
  //                "Mat dir: {" << materialDirection[0] << ", " << materialDirection[1] << ", " << materialDirection[2] << "}\n" << 
  //                "D: {" << D[0] << ", " << D[1] << ", " << D[2] << ", " << D[3] << ", " << D[4] << ", " << D[5] << "}\n" <<
  //                "Ddt: {" << D[0]*timeIncrement << ", " << D[1]*timeIncrement << ", " << D[2]*timeIncrement << ", " << D[3]*timeIncrement << ", " << D[4]*timeIncrement << ", " << D[5]*timeIncrement << "}\n" <<
  //                "Old stress: {" << oldStress[0] << ", " << oldStress[1] << ", " << oldStress[2] << ", " << oldStress[3] << ", " << oldStress[4] << ", " << oldStress[5] << "}\n" <<
  //                "Stress Incr: {" << stressIncrement[0] << ", " << stressIncrement[1] << ", " << stressIncrement[2] << ", " << stressIncrement[3] << ", " << stressIncrement[4] << ", " << stressIncrement[5] << "}\n" <<  
  //                "New stress: {" << newStress[0] << ", " << newStress[1] << ", " << newStress[2] << ", " << newStress[3] << ", " << newStress[4] << ", " << newStress[5] << "}\n");
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::computeTransverselyIsotropicPlasticStrainIncrement( real64 const ( & velocityGradient )[3][3],          // Velocity gradient tensor
                                                                          real64 const ( & oldStress )[6], // Stress at start of step
                                                                          real64 const ( & newStress )[6], // Stress at end of step
                                                                          const real64 Ez,             // Elastic modulus preferred direction
                                                                          const real64 Ep,             // Elastic modulus transverse plane
                                                                          const real64 nuzp,           // Poisson ratio coupled
                                                                          const real64 nup,            // Poisson ratio transverse plane
                                                                          const real64 Gzp,            // Shear modulus coupled plane
                                                                          real64 const ( & materialDirection )[3],			// material direction (unit vector)
                                                                          const real64 timeIncrement,
                                                                          real64 ( & plasticStrainIncrement )[6] ) const
{  // For hypo-elastic transversely isotropic models we compute the increment in plastic strain assuming
	// for some increment in total strain and stress and elastic properties.

	// New stress minus old stress
  real64 stressIncrement[6] = { 0 };
  LvArray::tensorOps::copy< 6 >( stressIncrement, newStress );
  LvArray::tensorOps::subtract< 6 >( stressIncrement, oldStress );

	// These are the TI elastic compliance coefficients using Brannon's TI basis tensors:
	real64 s1 = 1.0/Ez;
	real64 s2 = -nup/Ep;
	real64 s3 = -nuzp/Ez;
	real64 s4 = (1+nup)/Ep;
	real64 s5 = 1.0/(2.0*Gzp);

//Could be issue with using symmetric voigt vector for strain here, try using dense tensor
	// real64 elasticStrainIncrement[6] = { 0 };
  real64 elasticStrainIncrementDense[3][3] = { { 0 } };
  int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int p=0; p<3; p++)
			{
				for(int w=0; w<3; w++)
				{
					// elasticStrainIncrement[voigtMap[i][j]] 
          elasticStrainIncrementDense[i][j] += ( s1 * transverselyIsotropicB1(materialDirection,i,j,p,w) +
                                                 s2 * transverselyIsotropicB2(materialDirection,i,j,p,w) +
                                                 s3 * transverselyIsotropicB3(materialDirection,i,j,p,w) +
                                                 s4 * transverselyIsotropicB4(materialDirection,i,j,p,w) +
                                                 s5 * transverselyIsotropicB5(materialDirection,i,j,p,w) ) * stressIncrement[voigtMap[p][w]];
				}
			}
		}
	}

  // GEOS_LOG_RANK( "Elastic Strain Incr: {" << elasticStrainIncrement[0] << ", " << elasticStrainIncrement[1] << ", " << elasticStrainIncrement[2] << ", " << elasticStrainIncrement[3] << ", " << elasticStrainIncrement[4] << ", " << elasticStrainIncrement[5] << "}" );

	// plastic strain increment = Total strain increment - elastic strain increment
  real64 plasticStrainIncrementDense [3][3] = { { 0 } };
	for( int i = 0; i < 3; ++i )
	{
		for( int j = 0; j < 3; ++j )
		{
      //Old GEOS code did not index plasticStrainIncrement, how did this not throw an error
      //does the voigt notation include the factore of 2 for off axis entries of the tensor?
			// plasticStrainIncrement[voigtMap[i][j]] += 0.5 * ( velocityGradient[i][j] + velocityGradient[j][i] ) * timeIncrement - elasticStrainIncrement[voigtMap[i][j]];
      plasticStrainIncrementDense[i][j] += 0.5 * ( velocityGradient[i][j] + velocityGradient[j][i] ) * timeIncrement - elasticStrainIncrementDense[i][j];
		}
	}

   LvArray::tensorOps::denseToSymmetric<3>( plasticStrainIncrement, plasticStrainIncrementDense );
   // Strain off diagonal elements stored in voigt notation (2x symmetric elements)
   plasticStrainIncrement[3] *= 2.0;
   plasticStrainIncrement[4] *= 2.0;
   plasticStrainIncrement[5] *= 2.0;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB1( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B1_ijkl component, of the transversely isotropic basis tensor B1
  // described in Brannon's rotation/tensor book.

	real64 B1 = materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B1;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB2( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B2_ijkl component, of the transversely isotropic basis tensor B2
  // described in Brannon's rotation/tensor book.

	real64 B2 = delta(i,j)*delta(p,w) - 
               ( materialDirection[i]*materialDirection[j]*delta(p,w) + 
                 delta(i,j)*materialDirection[p]*materialDirection[w] ) + 
                 materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B2;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB3( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B3_ijkl component, of the transversely isotropic basis tensor B3
  // described in Brannon's rotation/tensor book.

	real64 B3 = materialDirection[i]*materialDirection[j]*delta(p,w) + 
                delta(i,j)*materialDirection[p]*materialDirection[w] - 
                2.0*materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B3;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB4( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B4_ijkl component, of the transversely isotropic basis tensor B4
  // described in Brannon's rotation/tensor book.

	real64 B4 = 0.5*( delta(i,p)*delta(j,w) + delta(i,w)*delta(j,p) ) -
              0.5*( delta(i,w)*materialDirection[j]*materialDirection[p] + 
                    materialDirection[i]*delta(j,p)*materialDirection[w] + 
                    materialDirection[i]*delta(j,w)*materialDirection[p] + 
                    delta(i,p)*materialDirection[j]*materialDirection[w] ) + 
              materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
	return B4;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB5( real64 const (& materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B5_ijkl component, of the transversely isotropic basis tensor B5
  // described in Brannon's rotation/tensor book.

	real64 B5 = 0.5*( delta(i,p)*materialDirection[j]*materialDirection[w] + 
                    delta(i,w)*materialDirection[j]*materialDirection[p] + 
                    materialDirection[i]*delta(j,w)*materialDirection[p] + 
                    materialDirection[i]*delta(j,p)*materialDirection[w] ) - 
              2.0*materialDirection[i]*materialDirection[p]*materialDirection[j]*materialDirection[w];
	return B5;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::delta( const int i,
		                       const int j ) const
{ // Dirac delta function that returns a real
	return i == j ? 1.0 : 0.0;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::slopePoint0( const real64 x,
                                     const real64 x1,
                                     const real64 x2,
                                     const real64 y1,
                                     const real64 y2,
                                     const real64 m1 ) const
{
  // returns a function h(x), for x1 <= x <= x2, with h(x1)=y1, h(x2)=y2, h'(x1)=m1, x'(x2)=0
  //realT beta = (m1*x1 - m1*x2 - y1 + y2)/(y1 - y2);
  //return -(((m1*(x1 - x2) + m1*powf((x - x2)/(x1 - x2),beta)*(-x + x2) + (m1*(-x1 + x2)*y1)/(y1 - y2))*(y1 - y2))/(m1*(x1 - x2)));
  // This should be the same, but simplified:
  return y2 + (y1-y2)*powf( (x-x2)/(x1-x2) , m1*(x1-x2)/(y1-y2) );
}


/**
 * @class Graphite
 *
 * Graphite material model.
 */
class Graphite : public SolidBase
{
public:

  /// @typedef Alias for GraphiteUpdates
  using KernelWrapper = GraphiteUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Graphite( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Graphite() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Graphite";

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
    /// string/key for transverse Young's modulus
    static constexpr char const * defaultYoungModulusTransverseString() { return "defaultYoungModulusTransverse"; }

    /// string/key for axial Young's modulus
    static constexpr char const * defaultYoungModulusAxialString() { return "defaultYoungModulusAxial"; }

    /// string/key for transverse Poisson's Ratio
    static constexpr char const * defaultPoissonRatioTransverseString() { return "defaultPoissonRatioTransverse"; }

    /// string/key for axial Poisson's Ratio
    static constexpr char const * defaultPoissonRatioAxialTransverseString() { return "defaultPoissonRatioAxialTransverse"; }

    /// string/key for transverse shear modulus
    static constexpr char const * defaultShearModulusAxialTransverseString() { return "defaultShearModulusAxialTransverse"; }

        /// string/key for default transverse Young's modulus presssure derivative
    static constexpr char const * defaultYoungModulusTransversePressureDerivativeString() { return "defaultYoungModulusTransversePressureDerivative"; }

    /// string/key for default axial Young's modulus presssure derivative
    static constexpr char const * defaultYoungModulusAxialPressureDerivativeString() { return "defaultYoungModulusAxialPressureDerivative"; }

    /// string/key for default axial Young's modulus presssure derivative
    static constexpr char const * defaultShearModulusAxialTransversePressureDerivativeString() { return "defaultShearModulusAxialTransversePressureDerivative"; }

    //string/key for element/particle velocityGradient value
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for quadrature point plasticStrain value 
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point relaxation value
    static constexpr char const * relaxationString() { return "relaxation"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }

    /// string/key for element/particle length scale
    static constexpr char const * strengthScaleString() { return "strengthScale"; }

    /// string/key for maximum strength
    static constexpr char const * failureStrengthString() { return "failureStrength"; }

    /// string/key for crack speed
    static constexpr char const * crackSpeedString() { return "crackSpeed"; }

    /// string/key for damaged material frictional slope
    static constexpr char const * damagedMaterialFrictionalSlopeString() { return "damagedMaterialFrictionalSlope"; }
    
    // string/key for distortion shear response parameter x2
    static constexpr char const * distortionShearResponseX2String() { return "distortionShearResponseX2"; }
    
    // string/key for distortion shear response parameter y1
    static constexpr char const * distortionShearResponseY1String() { return "distortionShearResponseY1"; }
    
    // string/key for distortion shear response parameter y2
    static constexpr char const * distortionShearResponseY2String() { return "distortionShearResponseY2"; }
    
    // string/key for distortion shear response parameter m1
    static constexpr char const * distortionShearResponseM1String() { return "distortionShearResponseM1"; }
    
    // string/key for in plane shear response x2
    static constexpr char const * inPlaneShearResponseX2String() { return "inPlaneShearResponseX2"; }
    
    // string/key for in plane shear response y1
    static constexpr char const * inPlaneShearResponseY1String() { return "inPlaneShearResponseY1"; }
    
    // string/key for in plane shear response y2
    static constexpr char const * inPlaneShearResponseY2String() { return "inPlaneShearResponseY2"; }
    
    // string/key for in plane shear response m1
    static constexpr char const * inPlaneShearResponseM1String() { return "inPlaneShearResponseM1"; }

    // string/key for coupled shear response x2
    static constexpr char const * coupledShearResponseX2String() { return "coupledShearResponseX2"; }
    
    // string/key for coupled shear response y1 
    static constexpr char const * coupledShearResponseY1String() { return "coupledShearResponseY1"; }
    
    // string/key for coupled shear response y2 
    static constexpr char const * coupledShearResponseY2String() { return "coupledShearResponseY2"; }
    
    // string/key for coupled shear response m1 
    static constexpr char const * coupledShearResponseM1String() { return "coupledShearResponseM1"; }
    
    // string/key for maximum plastic strain
    static constexpr char const * maximumPlasticStrainString() { return "maximumPlasticStrain"; }
    
    /// string/key for effective bulk modulus 
    static constexpr char const * effectiveBulkModulusString() { return "effectiveBulkModulus"; }

    /// string/key for effective shear modulus 
    static constexpr char const * effectiveShearModulusString() { return "effectiveShearModulus"; }

    /// string/key for material direction value
    static constexpr char const * materialDirectionString() { return "materialDirection"; }
  };

  /**
   * @brief Create a instantiation of the GraphiteUpdate class that refers to the data in this.
   * @return An instantiation of GraphiteUpdate.
   */
  GraphiteUpdates createKernelUpdates() const
  {
    return GraphiteUpdates( m_defaultYoungModulusTransverse,
                            m_defaultYoungModulusAxial,
                            m_defaultPoissonRatioTransverse,
                            m_defaultPoissonRatioAxialTransverse,
                            m_defaultShearModulusAxialTransverse,
                            m_effectiveBulkModulus,
                            m_effectiveShearModulus,
                            m_materialDirection,
                            m_defaultYoungModulusTransversePressureDerivative,
                            m_defaultYoungModulusAxialPressureDerivative,
                            m_defaultShearModulusAxialTransversePressureDerivative,
                            m_velocityGradient,
                            m_plasticStrain,
                            m_relaxation,
                            m_damage,
                            m_jacobian,
                            m_lengthScale,
                            m_strengthScale,
                            m_failureStrength,
                            m_crackSpeed,
                            m_damagedMaterialFrictionalSlope,
                            m_distortionShearResponseX2,
                            m_distortionShearResponseY1,
                            m_distortionShearResponseY2,
                            m_distortionShearResponseM1,
                            m_inPlaneShearResponseX2,
                            m_inPlaneShearResponseY1,
                            m_inPlaneShearResponseY2,
                            m_inPlaneShearResponseM1,
                            m_coupledShearResponseX2,
                            m_coupledShearResponseY1,
                            m_coupledShearResponseY2,
                            m_coupledShearResponseM1,
                            m_maximumPlasticStrain,
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
                          m_defaultYoungModulusTransverse,
                          m_defaultYoungModulusAxial,
                          m_defaultPoissonRatioTransverse,
                          m_defaultPoissonRatioAxialTransverse,
                          m_defaultShearModulusAxialTransverse,
                          m_effectiveBulkModulus,
                          m_effectiveShearModulus,
                          m_materialDirection,
                          m_defaultYoungModulusTransversePressureDerivative,
                          m_defaultYoungModulusAxialPressureDerivative,
                          m_defaultShearModulusAxialTransversePressureDerivative,
                          m_velocityGradient,
                          m_plasticStrain,
                          m_relaxation,
                          m_damage,
                          m_jacobian,
                          m_lengthScale,
                          m_strengthScale,
                          m_failureStrength,
                          m_crackSpeed,
                          m_damagedMaterialFrictionalSlope,
                          m_distortionShearResponseX2,
                          m_distortionShearResponseY1,
                          m_distortionShearResponseY2,
                          m_distortionShearResponseM1,
                          m_inPlaneShearResponseX2,
                          m_inPlaneShearResponseY1,
                          m_inPlaneShearResponseY2,
                          m_inPlaneShearResponseM1,
                          m_coupledShearResponseX2,
                          m_coupledShearResponseY1,
                          m_coupledShearResponseY2,
                          m_coupledShearResponseM1,
                          m_maximumPlasticStrain,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }

  /**
   * @brief Getter for default transverse Young's modulus
   * @return The value of the default transverse Young's modulus.
   */
  real64 getDefaultYoungModulusTransverse() const
  {
    return m_defaultYoungModulusTransverse;
  }

  /**
   * @brief Setter for the default transverse Young's modulus.
   * @param[in] input New value for the default transverse Young's modulus
   */
  void setDefaultYoungModulusTransverse( real64 const input )
  {
    m_defaultYoungModulusTransverse = input;
  }

  /**
   * @brief Getter for default axial Young's modulus
   * @return The value of the default axial Young's modulus.
   */
  real64 getDefaultYoungModulusAxial() const
  {
    return m_defaultYoungModulusAxial;
  }

  /**
   * @brief Setter for the default axial Young's modulus.
   * @param[in] input New value for the default axial Young's modulus
   */
  void setDefaultYoungModulusAxial( real64 const input )
  {
    m_defaultYoungModulusAxial = input;
  }

  /**
   * @brief Getter for default transverse Poisson's ratio
   * @return The value of the default transverse Poisson's ratio.
   */
  real64 getDefaultPoissonRatioTransverse() const
  {
    return m_defaultPoissonRatioTransverse;
  }

  /**
   * @brief Setter for the default transverse Poisson's ratio.
   * @param[in] input New value for the default transverse Poisson's ratio
   */
  void setDefaultPoissonRatioTransverse( real64 const input )
  {
    m_defaultPoissonRatioTransverse = input;
  }

  /**
   * @brief Getter for default axial Poisson's ratio
   * @return The value of the default axial/transverse Poisson's modulus.
   */
  real64 getDefaultPoissonRatioAxialTransverse() const
  {
    return m_defaultPoissonRatioAxialTransverse;
  }

  /**
   * @brief Setter for the default axial Poisson's modulus.
   * @param[in] input New value for the default axial/transverse Poisson's
   *             modulus
   */
  void setDefaultPoissonRatioAxialTransverse( real64 const input )
  {
    m_defaultPoissonRatioAxialTransverse = input;
  }

  /**
   * @brief Getter for default axial/transverse Shear modulus
   * @return The value of the default axial/transverse Shear modulus.
   */
  real64 getDefaultShearModulusAxialTransverse() const
  {
    return m_defaultShearModulusAxialTransverse;
  }

  /**
   * @brief Setter for the default axial/transverse Shear modulus.
   * @param[in] input New value for the default axial/transverse Shear modulus
   */
  void setDefaultShearModulusAxialTransverse( real64 const input )
  {
    m_defaultShearModulusAxialTransverse = input;
  }

  /**
   * @brief Accessor for effective bulk modulus
   * @return A const reference to arrayView1d<real64> containing the effective bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const effectiveBulkModulus() { return m_effectiveBulkModulus; }

  /**
   * @brief Const accessor for effective bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         effective bulk modulus (at every element).
   */
  arrayView1d< real64 const > const effectiveBulkModulus() const { return m_effectiveBulkModulus; }

 /**
   * @brief Accessor for effective bulk modulus
   * @return A const reference to arrayView1d<real64> containing the effective bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const effectiveShearModulus() { return m_effectiveShearModulus; }

  /**
   * @brief Const accessor for effective shear modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         effective shear modulus (at every element).
   */
  arrayView1d< real64 const > const effectiveShearModulus() const { return m_effectiveShearModulus; }

  /**
   * @brief Getter for effective bulk modulus.
   * @return reference to mutable effective bulk modulus.
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getEffectiveBulkModulus() const { return m_effectiveBulkModulus; }

  /**
  * @brief Getter for effective shear modulus.
  * @return reference to mutable effective shear modulus.
  */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getEffectiveShearModulus() const { return m_effectiveShearModulus; }


protected:
  virtual void postInputInitialization() override;

  /// The default value of the transverse Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungModulusTransverse;

  /// The default value of the axial Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungModulusAxial;

  /// The default value of the transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonRatioTransverse;

  /// The default value of the axial/transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonRatioAxialTransverse;

  /// The default value of the axial/transverse Shear modulus for any new
  /// allocations.
  real64 m_defaultShearModulusAxialTransverse;

  /// The effective bulk modulus.
  array1d< real64 > m_effectiveBulkModulus;

  /// The effective shear modulus
  array1d< real64 > m_effectiveShearModulus;
    
  /// State variable: The material direction for each element/particle
  array2d< real64 > m_materialDirection;

  /// The default value of the transverse Young's modulus pressure derivative for new allocations.
  real64 m_defaultYoungModulusTransversePressureDerivative;

  /// The default value of the axial Young's modulus pressure derivative for new allocations.
  real64 m_defaultYoungModulusAxialPressureDerivative;

  /// The default value of the axial transverse Shear modulus pressure derivative for new allocations.
  real64 m_defaultShearModulusAxialTransversePressureDerivative;

  ///State variable: The velocity gradient for each element/particle
  array3d< real64 > m_velocityGradient;

  ///State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  /// State variable: The relaxation values for each quadrature point
  array2d< real64 > m_relaxation;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Discretization-sized variable: The strength scale for each element/particle
  array1d< real64 > m_strengthScale;

  /// Material parameter: The value of the failure strength
  real64 m_failureStrength;

  /// Material parameter: The value of crack speed
  real64 m_crackSpeed;

  /// Material parameter: The value of the damaged material frictional slope
  real64 m_damagedMaterialFrictionalSlope;

  /// Material parameters: The values controlling the failure envelope
  real64 m_distortionShearResponseX2;
  real64 m_distortionShearResponseY1;
  real64 m_distortionShearResponseY2;
  real64 m_distortionShearResponseM1;

  real64 m_inPlaneShearResponseX2;
  real64 m_inPlaneShearResponseY1;
  real64 m_inPlaneShearResponseY2;
  real64 m_inPlaneShearResponseM1;

  real64 m_coupledShearResponseX2;
  real64 m_coupledShearResponseY1;
  real64 m_coupledShearResponseY2;
  real64 m_coupledShearResponseM1;
  
  real64 m_maximumPlasticStrain;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_GRAPHITE_HPP_ */
