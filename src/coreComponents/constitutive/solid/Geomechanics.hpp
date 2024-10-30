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
 * @file Geomechanics.hpp
 * @brief To Do
 *
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_GEOMECHANICS_HPP_
#define GEOS_CONSTITUTIVE_SOLID_GEOMECHANICS_HPP_

#include "SolidBase.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class GeomechanicsUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class GeomechanicsUpdates : public SolidBaseUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] b0 The value for the tangent elastic bulk modulus paramter 0
   * @param[in] b1 The value for the tangent elastic bulk modulus paramter 1
   * @param[in] b2 The value for the tangent elastic bulk modulus paramter 2
   * @param[in] b3 The value for the tangent elastic bulk modulus paramter 3
   * @param[in] b4 The value for the tangent elastic bulk modulus paramter 4
   * @param[in] g0 The value for the tangent elastic shear modulus paramter 0
   * @param[in] g1 The value for the tangent elastic shear modulus paramter 1
   * @param[in] g2 The value for the tangent elastic shear modulus paramter 2
   * @param[in] g3 The value for the tangent elastic shear modulus paramter 3
   * @param[in] g4 The value for the tangent elastic shear modulus paramter 4
   * @param[in] p0 The value for the crush curve paramter 0
   * @param[in] p1 The value for the crush curve paramter 1
   * @param[in] p2 The value for the crush curve paramter 2
   * @param[in] p3 The value for the crush curve paramter 3
   * @param[in] p4 The value for the crush curve paramter 4
   * @param[in] peakT1 The value for the peak T1 shear limit parameter
   * @param[in] fSlope The value for the F slope shear limit parameter
   * @param[in] stren The value for the stren shear limit paramter
   * @param[in] ySlope The value for the Y Slope shear limit parameter
   * @param[in] beta The value for the nonassociativity parameter
   * @param[in] t1RateDependence The value for the rate dependence parameter 1
   * @param[in] t2RateDependence The value for the rate dependence parameter 2
   * @param[in] fractureEnergyReleaseRate The value for the fracture energy release rate
   * @param[in] cr The value for the cap shape paramter
   * @param[in] fluidBulkModulus The value for the fluid bulk modulus
   * @param[in] fluidInitialPressure The value for the fluid initial pressure
   * @param[in] bulkModulus The ArrayView holding the bulk modulus for each element/particle.
   * @param[in] shearModulus The ArrayView holding the shear modulus for each element/particle.
   * @param[in] velocityGradient The ArrayView holding the velocity gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point.
   * @param[in] porosity The ArrayView holding the porosity for each quadrature point.
   * @param[in] damage The ArrayView holding the damage for each quardrature point.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   * @param[in] density The ArrayView holding the density data for each quadrature point.
   * @param[in] wavespeed The ArrayView holding the wavespeed data for each quadrature point.
   * @param[in] disableInelasticity Flag to disable plastic response/
   */
  GeomechanicsUpdates( real64 const & b0,
                       real64 const & b1,
                       real64 const & b2,
                       real64 const & b3,
                       real64 const & b4,
                       real64 const & g0,
                       real64 const & g1,
                       real64 const & g2,
                       real64 const & g3,
                       real64 const & g4,  
                       real64 const & p0,
                       real64 const & p1,
                       real64 const & p2,
                       real64 const & p3,
                       real64 const & p4,
                       real64 const & peakT1,
                       real64 const & fSlope,
                       real64 const & stren,
                       real64 const & ySlope,
                       real64 const & beta,
                       real64 const & t1RateDependence,
                       real64 const & t2RateDependence,
                       real64 const & fractureEnergyReleaseRate,
                       real64 const & cr,
                       real64 const & fluidBulkModulus,
                       real64 const & fluidInitialPressure,
                       int const & creep,
                       real64 const & creepC0,
                       real64 const & creepC1,
                       real64 const & creepA,
                       real64 const & creepB,
                       real64 const & creepC,
                       real64 const & strainHardeningN,
                       real64 const & strainHardeningK,
                       arrayView1d< real64 > const bulkModulus,
                       arrayView1d< real64 > const shearModulus,
                       arrayView3d< real64 const > const velocityGradient,
                       arrayView3d< real64 > const plasticStrain,
                       arrayView2d< real64 > const porosity,
                       arrayView2d< real64 > const damage,
                       arrayView1d< real64 > const lengthScale,
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
    m_b0( b0 ),
    m_b1( b1 ),
    m_b2( b2 ),
    m_b3( b3 ),
    m_b4( b4 ),
    m_g0( g0 ),
    m_g1( g1 ),
    m_g2( g2 ),
    m_g3( g3 ),
    m_g4( g4 ),
    m_p0( p0 ),
    m_p1( p1 ),
    m_p2( p2 ),
    m_p3( p3 ),
    m_p4( p4 ),
    m_peakT1( peakT1 ),
    m_fSlope( fSlope ),
    m_stren( stren ),
    m_ySlope( ySlope ),
    m_beta( beta ),
    m_t1RateDependence( t1RateDependence ),
    m_t2RateDependence( t2RateDependence ),
    m_fractureEnergyReleaseRate( fractureEnergyReleaseRate ),
    m_cr( cr ),
    m_fluidBulkModulus( fluidBulkModulus ),
    m_fluidInitialPressure( fluidInitialPressure ),
    m_creep( creep ),
    m_creepC0( creepC0 ),
    m_creepC1( creepC1 ),
    m_creepA( creepA ),
    m_creepB( creepB ),
    m_creepC( creepC ),
    m_strainHardeningN( strainHardeningN ),
    m_strainHardeningK( strainHardeningK ),
    m_bulkModulus( bulkModulus ),
    m_shearModulus( shearModulus ),
    m_velocityGradient( velocityGradient ),
    m_plasticStrain( plasticStrain ),
    m_porosity( porosity ),
    m_damage( damage ),
    m_lengthScale( lengthScale )
  {
    GEOS_UNUSED_VAR( m_g3 );
    GEOS_UNUSED_VAR( m_g4 );  
    GEOS_UNUSED_VAR( m_p2 );
    GEOS_UNUSED_VAR( m_p4 );
    GEOS_UNUSED_VAR( m_t1RateDependence );
    GEOS_UNUSED_VAR( m_t2RateDependence );
  }

  /// Default copy constructor
  GeomechanicsUpdates( GeomechanicsUpdates const & ) = default;

  /// Default move constructor
  GeomechanicsUpdates( GeomechanicsUpdates && ) = default;

  /// Deleted default constructor
  GeomechanicsUpdates() = delete;

  /// Deleted copy assignment operator
  GeomechanicsUpdates & operator=( GeomechanicsUpdates const & ) = delete;

  /// Deleted move assignment operator
  GeomechanicsUpdates & operator=( GeomechanicsUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic; // TODO: typo in anistropic (fix in DiscOps PR)

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
  int computeStep( real64 const ( & D )[6],              
                   const real64 & Dt,
                   const real64 & lch,                    
                   const real64 & Zeta_n,                
                   const real64 & coher_n,               
                   const real64 & porosity_n,            
                   real64 const ( & sigma_n )[6],       
                   real64 const ( & ep_n )[6],          
                   real64 & Zeta_p,                      
                   real64 & coher_p,                     
                   real64 & porosity_p,                  
                   real64 ( & sigma_p )[6],   
                   real64 ( & ep_p )[6] ) const;  

  GEOS_HOST_DEVICE
  void computeElasticProperties( real64 & bulk,
                                 real64 & shear ) const;

  GEOS_HOST_DEVICE
  void computeElasticProperties( real64 const ( &stress )[6],
                                 real64 const ( &ep )[6],
                                 const real64 & p3_crush_curve,
                                 const real64 & fluid_pressure_initial,
                                 const real64 & Km,     
                                 const real64 & Kf,     
                                 const real64 & C1,     
                                 const real64 & ev0,    
                                 const real64 & phi_i,
                                 real64 & bulk,
                                 real64 & shear ) const;
  GEOS_HOST_DEVICE
  void computeHypoelasticTrialStress( const real64 dt,                 
                                      const real64 bulk,               
                                      const real64 shear,              
                                      real64 const ( & stress_old )[6],
                                      real64 const ( & D )[6],         
                                      real64 ( & stress_new )[6] ) const;    

  GEOS_HOST_DEVICE
  void computeInvariants( real64 const ( & stress )[6],
                          real64 ( & S )[6],
                          real64 & I1,
                          real64 & J2,
                          real64 & rJ2 ) const;

  GEOS_HOST_DEVICE
  void computeCoher( const real64 & lch,
                     const real64 & I1_trial, 
			               const real64 & I1_0,     
			               const real64 & d_evp,    
			               const real64 & Gf,       
			               const real64 & coher_old,
			               real64 & coher_new ) const;    

  GEOS_HOST_DEVICE
  real64 computeX( const real64 & evp,
		               const real64 & phi_i,
		               const real64 & Km,
		               const real64 & Kf,
		               const real64 & C1,
		               const real64 & ev0,
		               const real64 & fluid_pressure_initial ) const;

  GEOS_HOST_DEVICE
  int computeStepDivisions( const real64 & X,
		                        const real64 & Zeta,
		                        real64 const ( & ep )[6],
		                        real64 const ( & sigma_n )[6],
		                        real64 const ( & sigma_trial )[6],
		                        const real64 & fluid_B0,
		                        const real64 & fluid_pressure_initial,
		                        const real64 & Km,       
		                        const real64 & Kf,       
		                        const real64 & C1,       
		                        const real64 & ev0,      
		                        const real64 & phi_i ) const;   

  GEOS_HOST_DEVICE
  int computeSubstep( real64 const ( & D )[6],        
                      const real64 & dt,
                      const real64 & lch,              
                      real64 const ( & sigma_old )[6],
                      real64 const ( & ep_old )[6],   
                      const real64 & X_old,
                      const real64 & Zeta_old,
                      const real64 & coher_old,       
                      const real64 & fluid_pressure_initial,
                      const real64 & Km,
                      const real64 & Kf,
                      const real64 & C1,
                      const real64 & ev0,
                      const real64 & phi_i,                   
                      real64 ( & sigma_new )[6],           
                      real64 ( & ep_new )[6],              
                      real64 & X_new,                      
                      real64 & Zeta_new,                   
                      real64 & coher_new ) const;	                

  GEOS_HOST_DEVICE
  int nonHardeningReturn( const real64 & I1_trial,            
                          const real64 & rJ2_trial,          
                          real64 const ( & S_trial )[6],       
                          const real64 & I1_old,               
                          const real64 & rJ2_old,          
                          real64 const ( & S_old )[6],        
                          real64 const ( & d_e )[6],           
                          const real64 & X,                    
                          const real64 & Zeta,                 
                          const real64 & coher,     
                          const real64 & hardening,          
                          const real64 & bulk,                 
                          const real64 & shear,
                          real64 & I1_new,                     
                          real64 & rJ2_new,
                          real64 ( & S_new )[6],
                          real64 ( & d_ep_new )[6] ) const;          

  GEOS_HOST_DEVICE
  void transformedBisection( real64 & z_0,
                             real64 & r_0,
                             const real64 & z_trial,
                             const real64 & r_trial,
                             const real64 & X,
                             const real64 & Zeta,
                             const real64 & coher,
                             const real64 & a1,
                             const real64 & a2,
                             const real64 & a3,
                             const real64 & a4,
                             const real64 & r_to_rJ2 ) const;

  GEOS_HOST_DEVICE
  int transformedYieldFunction( const real64 & z,
                                const real64 & r,
                                const real64 & X,
                                const real64 & Zeta,
										            const real64 & coher,
										            const real64 & a1,
										            const real64 & a2,
										            const real64 & a3,
										            const real64 & a4,
                                const real64 & r_to_rJ2 ) const;

  GEOS_HOST_DEVICE
  int computeYieldFunction( const real64 & I1,
                            const real64 & rJ2,
                            const real64 & X,
                            const real64 & Zeta,
									          const real64 & coher,
									          const real64 & a1,
									          const real64 & a2,
									          const real64 & a3,
									          const real64 & a4 ) const;

  GEOS_HOST_DEVICE
  real64 computedZetadevp( real64 const & fluid_pressure_initial,
                           real64 const & Km,
                           real64 const & Kf,
                           real64 const & ev0,
                           real64 const & phi_i,
                           real64 Zeta,
                           real64 evp ) const;

  GEOS_HOST_DEVICE
  void computeLimitParameters( real64 & a1,
		                           real64 & a2,
		                           real64 & a3,
		                           real64 & a4,
		                           const real64 & coher,
                               const real64 & hardening ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    SolidBaseUpdates::saveConvergedState( k, q );
  }

private:
  real64 const & m_b0;
  real64 const & m_b1;
  real64 const & m_b2;
  real64 const & m_b3;
  real64 const & m_b4;
  real64 const & m_g0;
  real64 const & m_g1;
  real64 const & m_g2;
  real64 const & m_g3;
  real64 const & m_g4;  
  real64 const & m_p0;
  real64 const & m_p1;
  real64 const & m_p2;
  real64 const & m_p3;
  real64 const & m_p4;
  real64 const & m_peakT1;
  real64 const & m_fSlope;
  real64 const & m_stren;
  real64 const & m_ySlope;
  real64 const & m_beta;
  real64 const & m_t1RateDependence;
  real64 const & m_t2RateDependence;
  real64 const & m_fractureEnergyReleaseRate;
  real64 const & m_cr;
  real64 const & m_fluidBulkModulus;
  real64 const & m_fluidInitialPressure;
  int const & m_creep;
  real64 const & m_creepC0;
  real64 const & m_creepC1;
  real64 const & m_creepA;
  real64 const & m_creepB;
  real64 const & m_creepC;
  real64 const & m_strainHardeningN;
  real64 const & m_strainHardeningK;

  /// A reference to the ArrayView holding the bulk modulus for each element/particle.
  arrayView1d< real64 > const m_bulkModulus;
  
  /// A reference to the ArrayView holding the shear modulus for each element/particle.+
  arrayView1d< real64 > const m_shearModulus;

  /// A reference to the ArrayView holding the velocity gradient for each element/particle.
  arrayView3d< real64 const > const m_velocityGradient;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the porosity for each quadrature point.
  arrayView2d< real64 > const m_porosity;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;
};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdate not implemented for Geomechanics model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( & stress )[6],
                                             DiscretizationOps & stiffness ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for Geomechanics model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdateStressOnly overload not implemented for Geomechanics model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const & timeIncrement,
                                                    real64 const ( & beginningRotation )[3][3],
                                                    real64 const ( & endRotation )[3][3],
                                                    real64 const ( & strainIncrement )[6],
                                                    real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( strainIncrement );

  // call the constitutive model
  GeomechanicsUpdates::smallStrainUpdateHelper( k, 
                                                q, 
                                                timeIncrement,                
                                                beginningRotation,
                                                endRotation,
                                                stress );

  // save new stress and return
  saveStress( k, q, stress );
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::smallStrainUpdateHelper( localIndex const k,
                                                   localIndex const q,
                                                   real64 const timeIncrement,  
                                                   real64 const ( & beginningRotation )[3][3],
                                                   real64 const ( & endRotation )[3][3],
                                                   real64 ( & stress )[6] ) const
{   
  // GEOS_UNUSED_VAR( timeIncrement );
  // GEOS_UNUSED_VAR( endRotation );
  
  //  temp matrix for computing rotations.
  real64 tempMat[ 3 ][ 3 ]= { { 0 } };
  
  // Unrotate velocity gradient:
  real64 unrotatedVelocityGradient[3][3] = { { 0 } };
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( tempMat, beginningRotation, m_velocityGradient[k] );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedVelocityGradient, tempMat, beginningRotation );
  real64 unrotatedVelocityGradientTranspose[3][3] = { { 0 } };
  LvArray::tensorOps::transpose< 3, 3 >( unrotatedVelocityGradientTranspose, unrotatedVelocityGradient );

  // CC: Is there an LvArray operation to get the symmetric part of a matrix?
  // Symmetric part of unrotated velocity gradient, the strain "rate" for the step
  real64 denseD[3][3] = { { 0 } };
  LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
  LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
  LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );
  real64 D[6] = { 0 };
  LvArray::tensorOps::denseToSymmetric<3>( D, denseD );

  // unrotated beginning-of-step stress:  sigma_old_unrotated = R_old^T*sigma_old*R_old
  real64 oldStress[6] = { 0 };
  real64 denseStress[3][3] = { { 0 } };
  LvArray::tensorOps::copy< 6 >( oldStress, m_oldStress[k][q] ); 
  LvArray::tensorOps::symmetricToDense<3>( denseStress, oldStress );
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( tempMat, beginningRotation, denseStress );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( denseStress, tempMat, beginningRotation );
  LvArray::tensorOps::denseToSymmetric<3>( oldStress, denseStress );
  
  // unrotated beginning-of-step plastic strain
  real64 oldPlasticStrain[6] = {0.};
  real64 densePlasticStrain[3][3] = { { 0 } };
  LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] ); 
  LvArray::tensorOps::symmetricToDense<3>( densePlasticStrain, oldPlasticStrain );
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( tempMat, beginningRotation, densePlasticStrain );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( densePlasticStrain, tempMat, beginningRotation );
  LvArray::tensorOps::denseToSymmetric<3>( oldPlasticStrain, densePlasticStrain );
   
  // Update effective elastic properties
  m_bulkModulus[k] = m_b0 + m_b1;
  m_shearModulus[k] = m_g0;
  
  // MH:TODO set wavespeed based on actual pressure-dependent bulk modulus, not the high-pressure limit.
  // m_wavespeed[k][0] = LvArray::math::sqrt( ( m_bulkModulus[k] + (4.0/3.0) * m_shearModulus[k] ) / m_density[k][0] );
  m_wavespeed[k][0] = LvArray::math::sqrt( ( m_b0 + m_b1 + (4.0/3.0) * m_g0 ) / m_density[k][0] );  //MH: should this be [0] or [q]?

  real64 characteristicLength = m_lengthScale[k];
  real64 oldZeta = 0.0;
  real64 oldCoher = 1.0 - m_damage[k][q];
  real64 oldPorosity = m_porosity[k][q];
  real64 newZeta;
  real64 newCoher;
  real64 newPorosity;
  real64 newStress[6] = {0.};
  real64 newPlasticStrain[6] = {0.};

  int errorFlag;
  errorFlag = computeStep( D,               // strain "rate"
               timeIncrement,                     // time step (s)
               characteristicLength,                    // length scale
               oldZeta,                 // trace of isotropic backstress at start of step(t_n)
               oldCoher,                // scalar-valued coherence at start of step(t_n)
						   oldPorosity,             // scalar-valued coherence at start of step(t_n)
						   oldStress,         // unrotated stress at start of step(t_n)
               oldPlasticStrain,            // plastic strain at start of step(t_n)
               newZeta,                        // trace of isotropic backstress at end of step(t_n+1)
               newCoher,                       // scalar-valued coherence at end of step(t_n+1)
						   newPorosity,                    // scalar-valued porosity at end of step(t_n+1)
						   newStress,         // unrotated stress at end of step(t_n+1)
               newPlasticStrain );           // plastic strain at end of step(t_n+1)

  if (errorFlag == 1)
  {
    GEOS_ERROR( "geomechanics model failed to converge." );
  }

  // re-rotate end-of-step stress: sigma_new = R_new*sigma_new*R_new^T
  LvArray::tensorOps::symmetricToDense<3>( denseStress, newStress );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( tempMat, endRotation, denseStress );
  LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >( denseStress, tempMat, endRotation );
  LvArray::tensorOps::denseToSymmetric<3>( newStress, denseStress );
  LvArray::tensorOps::copy< 6 >( stress, newStress ); // Copy to material data:

  // re-rotate end-of-step plastic strain: ep_new = R_new*ep_new*R_new^T
  LvArray::tensorOps::symmetricToDense<3>( densePlasticStrain, newPlasticStrain );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( tempMat, endRotation, densePlasticStrain );
  LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >( densePlasticStrain, tempMat, endRotation );
  LvArray::tensorOps::denseToSymmetric<3>( newPlasticStrain, densePlasticStrain );
  LvArray::tensorOps::copy< 6 >( m_plasticStrain[k][q], newPlasticStrain ); // Copy to material data:

  m_damage[k][q] = 1. - newCoher; // Copy to material data:
  m_porosity[k][q] = newPorosity; // Copy to material data:

  //return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeStep( real64 const ( & D )[6],               // strain "rate"
                                      const real64 & Dt,                     // time step (s)
                                      const real64 & lch,                    // length scale
                                      const real64 & Zeta_n,                 // trace of isotropic backstress at start of step(t_n)
                                      const real64 & coher_n,                // scalar-valued coherence at start of step(t_n)
						                          const real64 & porosity_n,             // scalar-valued coherence at start of step(t_n)
						                          real64 const ( & sigma_n )[6],         // unrotated stress at start of step(t_n)
                                      real64 const ( & ep_n )[6],            // plastic strain at start of step(t_n)
                                      real64 & Zeta_p,                        // trace of isotropic backstress at end of step(t_n+1)
                                      real64 & coher_p,                       // scalar-valued coherence at end of step(t_n+1)
						                          real64 & porosity_p,                    // scalar-valued porosity at end of step(t_n+1)
						                          real64 ( & sigma_p )[6],         // unrotated stress at end of step(t_n+1)
                                      real64 ( & ep_p )[6]            // plastic strain at end of step(t_n+1)
) const
{
	
  // Figure out where to initialize and set these and use pressure dependence to get value for each
  // particle correctly.

  real64 fluid_B0 = m_fluidBulkModulus; // Fluid bulk modulus (K_f)
	real64 fluid_pressure_initial = m_fluidInitialPressure; // Zero strain Fluid Pressure (Pf0)

	// These  variables are computed from input parameters and are used throughout the code
	// The are evaluates here to avoid repeated computation, or to simplify expressions.

	// This phi_i value is not modified by the disaggregation strain, because
	// it is the same for all particles.  Thus disaggregation strain is
	// not supported when there is a pore fluid.
	real64 phi_i = 1.0 - LvArray::math::exp( -m_p3 );                                 // Initial porosity (inferred from crush curve, used for fluid model/
	real64 Km = m_b0 + m_b1;                                           // Matrix bulk modulus
	real64 Kf = fluid_B0;                                          // Fluid bulk modulus
	real64 C1 = Kf * ( 1.0 - phi_i ) + Km * ( phi_i );                     // Term to simplify the fluid model expressions
	real64 ev0 = Kf > 0 ? C1 * fluid_pressure_initial / ( Kf * Km ) : 0.0; // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)

  // All stress values within computeStep are quasistatic.
  int n,
      chimax = 1,                                   // max allowed subcycle multiplier
      chi = 1,                                      // subcycle multiplier
      stepFlag,                                     // 0/1 good/bad step
      substepFlag;                                  // 0/1 good/bad substep

  // int m_mat_geo_subcycling_characteristic_number = 256; // allowable subcycles

  real64 evp_n = LvArray::tensorOps::symTrace< 3 >( ep_n ); // ep_n.Trace(); // Volumetric plastic strain at start of step.
  
  // hydrostatic compressive strength at start of step(t_n)
  real64 X_n = computeX( evp_n,
                          phi_i,
                          Km,
                          Kf,
                          C1,
                          ev0,
                          fluid_pressure_initial );
  
  real64 dt,                                        // substep time increment
          X_old,                                     // X at start of substep
          X_new,                                     // X at end of substep
          Zeta_old = Zeta_n,                         // Zeta at start of substep
          Zeta_new,                                  // Zeta at end of substep
          coher_old, 							                  // coher at start of substep
          coher_new;  								                // coher at end of substep


  real64 sigma_old[6];                                // sigma at start of substep
  real64 sigma_new[6];                                // sigma at end of substep
  real64 ep_old[6];                                   // plastic strain at start of substep
  real64 ep_new[6];                                   // plastic strain at end of substep

  real64 bulk, shear;

  // (1) Define conservative elastic properties based on the high-pressure
  // limit of the bulk modulus function. The bulk modulus function can be
  // without stress and plastic strain arguments to return the
  // high pressure limit.  These values are used only for computing number
  // of substeps and stable time steps.
  computeElasticProperties( bulk, 
                            shear ); 
  //Compute the trial stress: [sigma_trial] = computeTrialStress(sigma_old,d_e,K,G)
  real64 sigma_trial[6] = { 0 };
  computeHypoelasticTrialStress( Dt,              // full time step
                                 bulk,            // bulk modulus
                                 shear,           // shear modulus
                                 sigma_n,         // stress at start of step.
                                 D,               // D=sym(L)
                                 sigma_trial );   // stress at end of step

  real64 I1_trial,
        J2_trial,
        rJ2_trial;

  real64 d_e[6],
         S_trial[6];

  computeInvariants( sigma_trial,
                     S_trial,
                     I1_trial,
                     J2_trial,
                     rJ2_trial);

  // (2) Determine the number of substeps (nsub) based on the magnitude of
  // the trial stress increment relative to the characteristic dimensions
  // of the yield surface.  Also compare the value of the pressure dependent
  // elastic properties as sigma_old and sigma_trial and adjust nsub if
  // there is a large change to ensure an accurate solution for nonlinear
  // elasticity even with fully elastic loading.
  int nsub = computeStepDivisions( X_n,
                                  Zeta_n,
                                  ep_n,
                                  sigma_n,
                                  sigma_trial,
                                  fluid_B0,
                                  fluid_pressure_initial,
                                  Km,
                                  Kf,
                                  C1,
                                  ev0,
                                  phi_i);

  if ( nsub < 0 ) { // nsub > m_mat_geo_subcycling_characteristic_number[p_mat]. Delete particle
    goto failedStep;
  }

  // (3) Compute a subdivided time step:
  stepDivide:
    dt = Dt / ( chi * nsub );
    // d_e = (D*dt);
    // d_e = D;
    // d_e *= dt;
    LvArray::tensorOps::copy< 6 >( d_e, D );
    LvArray::tensorOps::scale< 6 >( d_e, dt );

  // (4) Set {sigma_old,X_old,Zeta_old,ep_old}={sigma_n,X_n,Zeta_n,ep_n} to their
  // values at the start of the step, and set n = 1:
  X_old = X_n;
  Zeta_old = Zeta_n;
  coher_old = coher_n;
  // sigma_old = sigma_n;
  // ep_old    = ep_n;
  LvArray::tensorOps::copy< 6 >( sigma_old, sigma_n );
  LvArray::tensorOps::copy< 6 >( ep_old, ep_n );
  n = 1;

	// -------------------------------------------------------------------------------
	// Apply creep to relax deviatoric stress for the whole step, this will be
  // the starting point for the plasticity solution.
	if ( m_creep == 1 )
	{
		real64 c0 = m_creepC0;
		real64 c1 = m_creepC1;
    real64 A = m_creepA;  // volumetric creep rate parameter
    real64 B = m_creepB;  // volumetric creep rate parameter
    real64 C = m_creepC;  // volumetric creep rate parameter

    real64 stress_iso[6] = { 0 };
    real64 stress_dev[6] = { 0 };
    real64 iso_old;
    real64 J2_old;
    twoInvariant::stressDecomposition( sigma_old,
                                       iso_old,
                                       J2_old,
                                       stress_dev); //This gives unit vector in direction of dev stress
    // rescale to have actual deviator
    LvArray::tensorOps::scale< 6 >( stress_dev, LvArray::math::sqrt(2/3) * J2_old);
                                    
    stress_iso[0] = iso_old;
    stress_iso[1] = iso_old;
    stress_iso[2] = iso_old;

    real64 ep_iso[6] = { 0 };
    real64 ep_dev[6] = { 0 };
    real64 ep_iso_old;
    real64 ep_J2_old;
    twoInvariant::stressDecomposition( ep_old,
                                       ep_iso_old,
                                       ep_J2_old,
                                       ep_dev); //This gives unit vector in direction of dev p. strain 
    LvArray::tensorOps::scale< 6 >( ep_dev, LvArray::math::sqrt(2/3) * ep_J2_old);                                       
    ep_iso[0] = ep_iso_old;
    ep_iso[1] = ep_iso_old;
    ep_iso[2] = ep_iso_old;

   	real64 sigma_vm_old = LvArray::math::sqrt( 3.0 * J2_old );
		computeElasticProperties( bulk,
                              shear );
		real64 elasticVMShearStrain = sigma_vm_old / ( 3 * shear ); // This is equivalent to sqrt(2/3) * J2 invariant of sigma_dev/(2*shear)

		if ( elasticVMShearStrain > 1.e-12 )
		{  // only apply creep if there is elastic strain

      real64 equilibriumVMShearStrain = 0.0; // this could be some other function of pressure, granular density, etc.
      real64 creepVMStrainIncrement = c0 * powf( elasticVMShearStrain - equilibriumVMShearStrain, c1 ) * Dt;
			real64 devStressCreepScale = LvArray::math::max( elasticVMShearStrain - creepVMStrainIncrement , 0.0 ) / elasticVMShearStrain;

			// increment plastic strain such that total strain is constant:  (0.5/g0)*stress_dev*( 1.0 - devStressCreepScale );
			real64 creepStrainIncrement[6];  
			// creepStrainIncrement = stress_dev; // this is temporarily the dev stress
			// creepStrainIncrement *= 0.5/g0*( 1.0 - devStressCreepScale ); // this is now the change in elastic deviatoric strain
			// ep_old += creepStrainIncrement; // increment plastic strain such that total strain is constant
      LvArray::tensorOps::copy< 6 >( creepStrainIncrement, stress_dev );
      LvArray::tensorOps::scale< 6 >( creepStrainIncrement,  0.5/shear*( 1.0 - devStressCreepScale ) ); // CC: TODO check if new to double off diagonal components
      LvArray::tensorOps::add< 6 >( ep_dev, creepStrainIncrement );

			// relax elastic deviatoric stress due to creep
			// stress_dev *= devStressCreepScale;
      LvArray::tensorOps::scale< 6 >( stress_dev, devStressCreepScale );
		}

	  // volumetric creep ----------------


    real64 p = -iso_old;  // hydrostatic pressure, positive in compression

    // volumetric plastic strain, negative in compression
 	  real64 evp = ep_old[0] + ep_old[1] + ep_old[2]; 
    // unloaded porosity at the start of the step.
    real64 phi_p = LvArray::math::max( 1.e-10 , 1.0 + LvArray::math::exp(-evp)*( phi_i - 1 ) ); 
    // equilibrium porosity at the start of the step.
		real64 phi_e = LvArray::math::max(1.e-10 , A*LvArray::math::exp(-p/B));    

    // uncomment for debugging:
    //		  std::cout<<"evp = "<<evp<<", phi_p = "<<phi_p<<", phi_e = "<<phi_e<<std::endl;

    if ( (phi_p > phi_e) && (p > 1.e-12) && (C > 0) )
 		  {  // creep compaction
 			  real64 dphidt = -1.0*p*C*( phi_p - phi_e );  // creep compaction rate:
 			  real64 phi_c = LvArray::math::max( phi_e, phi_p + dphidt*dt ); // unloaded porosity after creep, don't let it go below equilibrium level
 			  real64 evp_c = LvArray::math::log( (phi_i - 1. ) / ( phi_c - 1. ) ); // vol. strain after creep.
 			  real64 devp = evp_c - evp;  // creep vol. strain increment
 			  real64 p_c = LvArray::math::max( 0., p + bulk*devp );  // relaxed pressure after creep.

        // uncomment for debugging:
        // real64 evp_0 = LvArray::math::log( (phi_i - 1. ) / ( phi_p - 1. ) ); // vol. strain before creep, computed from porosity (uncomment for debugging)
        // std::cout<<"creep compaction:phi_p = "<<phi_p<<", phi_e = "<<phi_e<<", phi_c = "<<phi_c<<", dphiDt = "<<dphidt
        // <<", evp_0 = "<<evp_0<<", evp = "<<evp<<", devp = "<<devp<<", bulk = "<<bulk
        // <<", p_old = "<<p<<", p_c = "<<p_c<<std::endl;

        // update stress and plastic strain values after creep
        stress_iso[0] = -1.0*p_c;
			  stress_iso[1] = -1.0*p_c;
        stress_iso[2] = -1.0*p_c;
 			  ep_iso[0] = evp_c/3.;
 			  ep_iso[1] = evp_c/3.;
 			  ep_iso[2] = evp_c/3.;

        // update cap function to be consistent with new vol. plastic
        X_old = computeX( evp_c,
        phi_i, // Initial porosity (inferred from crush curve, used for fluid model/
        bulk, // matrix bulk modulus
        0.0, // Fluid bulk modulus
        0.0, // Term to simplify the fluid model expressions
        0.0, // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)
        0.0 );
 		  }

	    // Reassemble the stress with a scaled deviatoric stress tensor
      LvArray::tensorOps::copy< 6 >( sigma_old, stress_dev );
      LvArray::tensorOps::add< 6 >( sigma_old, stress_iso );

      LvArray::tensorOps::copy< 6 >( ep_old, ep_dev );
      LvArray::tensorOps::add< 6 >( ep_old, ep_iso );
	}
	// -------------------------------------------------------------------------------

  // (5) Call substep function {sigma_new,ep_new,X_new,Zeta_new}
  //                               = computeSubstep(D,dt,sigma_old,ep_old,X_old,Zeta_old)
  computeSubstep:
    substepFlag = computeSubstep( D,
                                  dt,
                                  lch,
                                  sigma_old,
                                  ep_old,
                                  X_old,
                                  Zeta_old,
                                  coher_old,
                                  fluid_pressure_initial,
                                  Km,
                                  Kf,
                                  C1,
                                  ev0,
                                  phi_i,
                                  sigma_new,
                                  ep_new,
                                  X_new,
                                  Zeta_new,
                                  coher_new );


  // (6) Check error flag from substep calculation:
  if (substepFlag == 0) { // no errors in substep calculation
    if (n < (chi*nsub)) { // update and keep substepping
      X_old     = X_new;
      Zeta_old  = Zeta_new;
      coher_old = coher_new;
      // sigma_old = sigma_new;
      // ep_old    = ep_new;
      LvArray::tensorOps::copy< 6 >( sigma_old, sigma_new );
      LvArray::tensorOps::copy< 6 >( ep_old, ep_new );
      n++;
      goto computeSubstep;
    }
    else goto successfulStep; // n = chi*nsub, Step is done
  }
  else
  {  // failed substep
    if (chi < chimax)   {       // errors in substep calculation
      chi = 2*chi;
      goto stepDivide;
    }
    else goto failedStep;     // bad substep and chi>=chimax, Step failed even with substepping
  }

  // (7) Successful step, set value at end of step to value at end of last substep.
  successfulStep:
    Zeta_p    = Zeta_new;
    coher_p = coher_new;
    // sigma_p   = sigma_new;
    // ep_p      = ep_new;
    LvArray::tensorOps::copy< 6 >( sigma_p, sigma_new );
    LvArray::tensorOps::copy< 6 >( ep_p, ep_new );
    // porosity_p = 1 + LvArray::math::exp( ep_new.Trace() )*( phi_i - 1 );
    porosity_p = 1 + LvArray::math::exp( LvArray::tensorOps::symTrace< 3 >( ep_new ) )*( phi_i - 1 );

    stepFlag  = 0;
    return stepFlag;

  // (8) Failed step, Send ParticleDelete Flag to Host Code, Store Inputs to particle data:
  failedStep:
    // input values for sigma_new,X_new,Zeta_new,ep_new, along with error flag
    Zeta_p    = Zeta_n;
    coher_p = coher_n;
    // sigma_p   = sigma_n;
    // ep_p      = ep_n;
    LvArray::tensorOps::copy< 6 >( sigma_p, sigma_n );
    LvArray::tensorOps::copy< 6 >( ep_p, ep_n );
    porosity_p = porosity_n;
    stepFlag  = 1;

    return stepFlag;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeHypoelasticTrialStress( const real64 dt,                    // time step
                                                         const real64 bulk,                  // bulk modulus
                                                         const real64 GEOS_UNUSED_PARAM(shear),                 // shear modulus
                                                         real64 const ( & stress_old )[6],   // stress at start of step.
                                                         real64 const ( & D )[6],            // D=sym(L)
                                                         real64 ( & stress_new )[6]         // stress at end of step
) const
{
  // Hypoelastic, isotropic, linear elasticity.
  real64 Diso[6] = { 0 }, 
         Ddev[6] = { 0 };

  // isotropicDeviatoricDecomposition( D, Diso, Ddev );
  LvArray::tensorOps::symAddIdentity< 3 >( Diso, LvArray::tensorOps::symTrace< 3 >( D )/3.0 );
  LvArray::tensorOps::copy< 6 >( Ddev, D );
  LvArray::tensorOps::subtract< 6 >( Ddev, Diso );

  // Cauchy trial stress:  sigma_n+1^trial = sigma_n + [ 3*K*iso(D) + 2*G*dev(D) ]*dt
  // Diso *= 3.*bulk*dt; // this is now the increment in isotropic stress
  // Ddev *= 2.*shear*dt; // this is now the increment in deviatoric stress
  // CC: TODO check if need to double off diagonal components for symmetric calculations
  LvArray::tensorOps::scale< 6 >( Diso, 3.*bulk*dt );
  LvArray::tensorOps::scale< 6 >( Ddev, 3.*bulk*dt );

  // stress_new = stress_old;
  // stress_new += Diso;
  // stress_new += Ddev;
  LvArray::tensorOps::copy< 6 >( stress_new, stress_old );
  LvArray::tensorOps::add< 6 >( stress_new, Diso );
  LvArray::tensorOps::add< 6 >( stress_new, Ddev );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeElasticProperties( real64 & bulk,
                                                    real64 & shear 
) const
{
    // When computeElasticProperties() is called with two real64s as arguments, it
    // computes the high pressure limit tangent elastic shear and bulk modulus
    // This is used to estimate wave speeds and make conservative estimates of substepping.
    shear = m_g0;       // Linear elastic shear Modulus
    bulk  = m_b0 + m_b1;  // Bulk Modulus
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeElasticProperties( real64 const ( &stress )[6],
                                                    real64 const ( &ep )[6],
		                                                const real64 & GEOS_UNUSED_PARAM( p3_crush_curve ),
		                                                const real64 & GEOS_UNUSED_PARAM( fluid_pressure_initial ),
		                                                const real64 & Km,             // Matrix bulk modulus
		                                                const real64 & Kf,             // Fluid bulk modulus
		                                                const real64 & C1,             // Term to simplify the fluid model expressions
		                                                const real64 & ev0,            // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)
		                                                const real64 & phi_i,			    // Initial porosity (inferred from crush curve, used for fluid model/
		                                                real64 & bulk,
		                                                real64 & shear
) const
{
  // Compute the nonlinear elastic tangent stiffness as a function of the pressure
  // plastic strain, and fluid parameters.
	real64 I1 = LvArray::tensorOps::symTrace< 3 >( stress ), //stress.Trace(),
		     evp = LvArray::tensorOps::symTrace< 3 >( ep ); //ep.Trace();

  // ..........................................................Undrained
  // The low pressure bulk and shear moduli are also used for the tensile response.
	bulk = m_b0;
	shear = m_g0;

	// To be thermodynamically consistent, the shear modulus in an isotropic model
	// must be constant, but the bulk modulus can depend on pressure.  However, this
	// leads to a Poisson's ratio that approaches 0.5 at high pressures, which is
	// inconsistent with experimental data for the Poisson's ratio, inferred from the
	// Young's modulus.  Induced anisotropy is likely the cause of the discrepency,
	// but it may be better to allow the shear modulus to vary so the Poisson's ratio
	// remains reasonable.
	//
	// If the user has specified a nonzero value of G1 and G2, the shear modulus will
	// vary with pressure so the drained Poisson's ratio transitions from G1 to G1+G2 as
	// the bulk modulus varies from B0 to B0+B1.  The fluid model further affects the
	// bulk modulus, but does not alter the shear modulus, so the pore fluid does
	// increase the Poisson's ratio.
	if( evp <= 0.0 )
    {
		if ( I1 < 0.0 )
        {
			real64 expb2byI1 = LvArray::math::exp( m_b2 / I1 );
			bulk = bulk + m_b1 * expb2byI1;
			if(m_g1 != 0.0 && m_g2 != 0.0){
				real64 nu = m_g1 + m_g2 * expb2byI1;
				shear = 1.5 * bulk * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu );
			}
		}

		// Elastic-plastic coupling
		if ( evp < 0.0 )
        {
            bulk = bulk - m_b3 * LvArray::math::exp( m_b4 / evp );
        }
	}


    // In  compression, or with fluid effects if the strain is more compressive
    // than the zero fluid pressure volumetric strain:
	if ( evp <= ev0 && Kf != 0.0 )
    {   // ..........................................................Undrained
		// Compute the porosity from the strain using Homel's simplified model, and
		// then use this in the Biot-Gassmann formula to compute the bulk modulus.

		// The dry bulk modulus, taken as the low pressure limit of the nonlinear
		// formulation:
		real64 Kd = m_b0;
		if ( evp < 0.0 )
        { 
            Kd = m_b0 - m_b3 * LvArray::math::exp( m_b4 / evp );
        }

		// Current unloaded porosity (phi):
		real64 C2 = LvArray::math::exp( evp * Km / C1 ) * phi_i;
		real64 phi = C2 / ( -LvArray::math::exp( evp * Kf / C1 ) * ( phi_i - 1.0) + C2 );

		// Biot-Gassmann formula for the saturated bulk modulus, evaluated at the
		// current porosity.  This introduces some error since the Kd term is a
		// function of the initial porosity, but since large strains would also
		// modify the bulk modulus through damage
		real64 oneminusKdbyKm = 1.0 - Kd / Km;
		bulk = Kd + oneminusKdbyKm * oneminusKdbyKm / ( ( oneminusKdbyKm - phi ) / Km + ( 1.0 / Kf - 1.0 / Km ) * phi );
	}
}

// [nsub] = computeStepDivisions(X,Zeta,ep,sigma_n,sigma_trial)
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeStepDivisions( const real64 & X,
		                                           const real64 & GEOS_UNUSED_PARAM( Zeta ),
		                                           real64 const ( & ep )[6],
		                                           real64 const ( & sigma_n )[6],
		                                           real64 const ( & sigma_trial )[6],
		                                           const real64 & GEOS_UNUSED_PARAM( fluid_B0 ),
		                                           const real64 & fluid_pressure_initial,
		                                           const real64 & Km,          // Matrix bulk modulus
		                                           const real64 & Kf,          // Fluid bulk modulus
		                                           const real64 & C1,          // Term to simplify the fluid model expressions
		                                           const real64 & ev0,         // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)
		                                           const real64 & phi_i			// Initial porosity (inferred from crush curve, used for fluid model/
) const
{
  // compute the number of step divisions (substeps) based on a comparison
  // of the trial stress relative to the size of the yield surface, as well
  // as change in elastic properties between sigma_n and sigma_trial.
  int subcycling_characteristic_number = 256;
  int nmax = LvArray::math::ceil(subcycling_characteristic_number);

  // Compute change in bulk modulus:
  real64 bulk_n,
         shear_n,
         bulk_trial,
         shear_trial;

  computeElasticProperties( sigma_n,
                            ep,
                            m_p3,
                            fluid_pressure_initial,
                            Km,
                            Kf,
                            C1,
                            ev0,
                            phi_i,
                            bulk_n,
                            shear_n );

  computeElasticProperties( sigma_trial,
                            ep,
                            m_p3,
                            fluid_pressure_initial,
                            Km,
                            Kf,
                            C1,
                            ev0,
                            phi_i,
                            bulk_trial,
                            shear_trial );
  int n_bulk = LvArray::math::ceil( fabs( bulk_n - bulk_trial ) / bulk_n );

  // Compute trial stress increment relative to yield surface size:
  real64 d_sigma[6] = { 0 };
  // d_sigma = sigma_trial; //sigma_trial - sigma_n;
  // d_sigma -= sigma_n;
  LvArray::tensorOps::copy< 6 >( d_sigma, sigma_trial );
  LvArray::tensorOps::subtract< 6 >( d_sigma, sigma_n );

  real64 size = 0.5*(m_peakT1 - X);

  if( m_stren > 0.0 )
  {
	  size = LvArray::math::min( size, m_stren );
  }

  d_sigma[3] *= 2.0;
  d_sigma[4] *= 2.0;
  d_sigma[5] *= 2.0;
  int n_yield = int( LvArray::math::ceil( 1.0e-4 * LvArray::tensorOps::l2Norm< 6 >( d_sigma ) / size ));
  d_sigma[3] *= 0.5;
  d_sigma[4] *= 0.5;
  d_sigma[5] *= 0.5;

  // nsub is the maximum of the two values.above.  If this exceeds allowable,
  // throw warning and delete particle.
  int nsub = LvArray::math::max( n_bulk, n_yield);

  if( nsub > subcycling_characteristic_number )
  {
    nsub = -1;
  }
  else
  {
    nsub = LvArray::math::min( LvArray::math::max( nsub, 1 ), nmax  );
  }
  return nsub;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeInvariants( real64 const ( & stress )[6],
                                             real64 ( & S )[6],
                                             real64 & I1,
                                             real64 & J2,
                                             real64 & rJ2 
) const
{
  // Compute the first invariants
//   I1 = LvArray::tensorOps::symTrace< 3 >( stress ); //stress.Trace();  //Pa

  // Compute the deviatoric part of the tensor
//   S = dev(stress); // stress - one_third*Identity*I1

//   real64 deviator[6] = { 0 };
  twoInvariant::stressDecomposition( stress,
                                     I1,
                                     J2, //this is rootJ2
                                     S); //this gives a unit verctor
    I1 *= 3.0;
    LvArray::tensorOps::scale< 6 >( S, LvArray::math::sqrt(2.0 / 3.0) * J2 );
    J2 = 3 * J2 * J2;  //MH: why are we squaring and then taking a square root?


//   // Compute the second invariant
//   J2 = computeJ2fromDevStress(S);  // 0.5*S.Contract(S);

  if( J2 < 1e-16*( I1 * I1 + J2 ) )
  {
    J2 = 0.0;
  }
  rJ2 = LvArray::math::sqrt(J2);
}

// update the coherence (1-damage) variable based on dilational plastic work
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeCoher( const real64 & lch,       // length scale
                                        const real64 & I1_trial,  // trial value of I1
			                                  const real64 & I1_0,      // I1 value on yield surface
			                                  const real64 & d_evp,     // increment in vol plastic strain
			                                  const real64 & Gf,        // fracture energy per unit area
			                                  const real64 & coher_old, // old coherence = 1-damage
			                                  real64 & coher_new      // OUPUT: new value of coher
) const
{
	coher_new = coher_old;
	if( Gf > 0 )
	{
		// real64 d_I1 = I1_trial - I1_0; // Seemed unused
		if ( d_evp > 0 && ( I1_trial - I1_0 ) > 0 )
		{
			// increment in work per unit volume.
			real64 d_dilationalPlasticWork = d_evp*0.5*(I1_trial - I1_0);
			// particle length scale
			// real64 lch = m_lengthScale[k]; //m_plane_strain ? LvArray::math::sqrt( m_dx * m_dx + m_dy * m_dy ) : LvArray::math::sqrt( m_dx * m_dx + m_dy * m_dy + m_dz * m_dz );
			// fracture energy
			real64 d_damage = ( d_dilationalPlasticWork*lch ) / Gf;
			coher_new = LvArray::math::max( 0.0, coher_old - d_damage );
		}
	}
}

// Compute state variable X, the Hydrostatic Compressive strength (cap position)
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GeomechanicsUpdates::computeX( const real64 & evp,
		                                  const real64 & phi_i, // Initial porosity (inferred from crush curve, used for fluid model/
		                                  const real64 & Km, // matrix bulk modulus
		                                  const real64 & Kf, // Fluid bulk modulus
		                                  const real64 & C1, // Term to simplify the fluid model expressions
		                                  const real64 & ev0, // Zero fluid pressure vol. strain.  (will equal zero if pfi=0)
		                                  const real64 & fluid_pressure_initial
) const
{
  // X is the value of (I1 - Zeta) at which the cap function crosses
  // the hydrostat. For the drained material in compression. X(evp)
  // is derived from the emprical Kayenta crush curve, but with p2 = 0.
  // In tension, M. Homel's piecewise formulation is used.

  // define and initialize some variables
  real64 X;
  real64 identity[6] = {1.0,1.0,1.0,0.0,0.0,0.0};
  // R2Tensor Identity = R2Tensor( 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

  if( evp <= -m_p3 ) { // ------------Plastic strain exceeds allowable limit--------------------------
    // The plastic strain for this iteration has exceed the allowable
    // value.  X is not defined in this region, so we set it to a large
    // negative number.
    //
    // The code should never have evp<-p3, but will have evp=-p3 if the
    // porosity approaches zero (within the specified tolerance).  By setting
    // X=1e12*p0, the material will respond as though there is no porosity.
    X = 1.0e12 * m_p0;
  }
  else { // ------------------Plastic strain is within allowable domain------------------------
    // We first compute the drained response.  If there are fluid effects, this value will
    // be used in determining the elastic volumetric strain to yield.
    if( evp <= 0.0 ){
      X = ( m_p0 * m_p1 + LvArray::math::log( ( evp + m_p3 ) / m_p3 ) ) / m_p1;
    }
    else{
      X = m_p0 * LvArray::math::pow( 1.0 + evp, 1.0 / ( m_p0 * m_p1 * m_p3 ) );
    }

    if( Kf !=0.0 && evp <= ev0 ) { // ------------------------------------------- Fluid Effects
      // First we evaluate the elastic volumetric strain to yield from the
      // empirical crush curve (Xfit) and bulk modulus (Kfit) formula for
      // the drained material.  Xfit was computed as X above.

      // Kfit is the drained bulk modulus evaluated at evp, and for I1 = Xdry/2.
      real64 Kdry = m_b0 + m_b1 * LvArray::math::exp( 2.0 * m_b2 / X );
      if ( evp < 0.0 )
      {
        Kdry = Kdry - m_b3 * LvArray::math::exp( m_b4 / evp );
      }

      // Now we use our engineering model for the bulk modulus of the
      // saturated material (Keng) to compute the stress at our elastic strain to yield.
      // Since the stress and plastic strain tensors are not available in this scope, we call the
      // computeElasticProperties function with isotropic matrices that will have the
      // correct values of evp. (The saturated bulk modulus doesn't depend on I1).
      real64 Ksat, 
             Gsat;       // Not used, but needed to call computeElasticProperties()
      // This needs to be evaluated at the current value of pressure.

      real64 hydroStress[6] = { 0 }, 
             hydroStrain[6] = { 0 };

      // hydroStress = Identity; //(1./6.)*X*Identity;
      // hydroStress *= ( 1.0 / 6.0 ) * X;
      LvArray::tensorOps::copy< 6 >( hydroStress, identity );
      LvArray::tensorOps::scale< 6 >( hydroStress, ( 1.0 / 6.0) * X );

      // hydroStrain = Identity; //(1./3.)*evp*Identity;
      // hydroStrain *= (1./3.) * evp;
      LvArray::tensorOps::copy< 6 >( hydroStrain, identity );
      LvArray::tensorOps::scale< 6 >( hydroStrain, ( 1.0 / 3.0) * evp );

      computeElasticProperties( hydroStress,
			                          hydroStrain,
			                          m_p3,
			                          fluid_pressure_initial,
                                Km,
                                Kf,
                                C1,
                                ev0,
                                phi_i,
			                          Ksat, 
                                Gsat ); // OUTPUT: saturated bulk and shear modulus
			
      // Compute the stress to hydrostatic yield.
      // We are only in this loop if(evp <= ev0)
      X = X * Ksat / Kdry;   // This is X_sat = K_sat*eve = K_sat*(X_dry/K_dry)
    } // End fluid effects
  } // End plastic strain in allowable domain
  return X;
}


// Computes the updated stress state for a substep
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeSubstep( real64 const ( & D )[6],         // strain "rate" for substep D=sym(L)
		                                     const real64 & dt,               // substep time interval
                                         const real64 & lch,              // length scale
		                                     real64 const ( & sigma_old )[6], // stress at start of substep
		                                     real64 const ( & ep_old )[6],    // plastic strain at start of substep
		                                     const real64 & X_old,            // hydrostatic compressive strength at start of substep
		                                     const real64 & Zeta_old,         // trace of isotropic backstress at start of substep
		                                     const real64 & coher_old,        // scalar valued coherence = 1-Damage
		                                     const real64 & fluid_pressure_initial,
		                                     const real64 & Km,
		                                     const real64 & Kf,
		                                     const real64 & C1,
		                                     const real64 & ev0,
		                                     const real64 & phi_i,
		                                     real64 ( & sigma_new )[6],            // stress at end of substep
		                                     real64 ( & ep_new )[6],               // plastic strain at end of substep
		                                     real64 & X_new,                       // hydrostatic compressive strength at end of substep
		                                     real64 & Zeta_new,                    // trace of isotropic backstress at end of substep
		                                     real64 & coher_new	                 // coherence at the end of step.
) const
{
  // Computes the updated stress state for a substep that may be either elastic, plastic, or
  // partially elastic.   Returns an integer flag 0/1 for a good/bad update.
  int substepFlag,
		  returnFlag;

  real64 one_third = 1.0 / 3.0;
  real64 identity[6] = {1.0,1.0,1.0,0.0,0.0,0.0};

  // (1)  Compute the elastic properties based on the stress and plastic strain at
  // the start of the substep.  These will be constant over the step unless elastic-plastic
  // is used to modify the tangent stiffness in the consistency bisection iteration.
  real64 bulk,
         shear;
  computeElasticProperties( sigma_old,
                            ep_old,
                            m_p3,
                            fluid_pressure_initial,
                            Km,
                            Kf,
                            C1,
                            ev0,
                            phi_i,
                            bulk,
                            shear);

  // (3) Compute the trial stress: [sigma_trail] = computeTrialStress(sigma_old,d_e,K,G)
  real64 sigma_trial[6],
         S_trial[6];

  computeHypoelasticTrialStress( dt,              // time step
  		                           bulk,            // bulk modulus
  		                           shear,           // shear modulus
  		                           sigma_old,       // stress at start of step.
  		                           D,               // D=sym(L)
  		                           sigma_trial );   // stress at end of step

  real64 I1_trial,
         J2_trial,
         rJ2_trial;
  computeInvariants( sigma_trial, 
                     S_trial,
                     I1_trial,
                     J2_trial,
                     rJ2_trial );

  // (4) Evaluate the yield function at the trial stress:
  // Compute the limit parameters based on the value of coher.  These are then passed down
  // to the computeYieldFunction, to avoid the expense of repeatedly computing a3

  // Compute strain hardening applied to the STREN parameter
  real64 hardening = 0.0;
   if (m_strainHardeningK > 0)
   {

	   // This experssion is incomplete until the (2./3.)*sqrt(equilivantPlasticStrain) is applied below.
	  //  real64 equilivantPlasticStrain = ep_old(0,0)*ep_old(0,0) + ep_old(1,1)*ep_old(1,1) + ep_old(2,2)*ep_old(2,2)
		// 	   + 3*ep_old(0,1)*ep_old(0,1) + 3*ep_old(0,2)*ep_old(0,2) + 3*ep_old(1,2)*ep_old(1,2)
		// 	   - ep_old(1,1)*ep_old(2,2) - ep_old(0,0)*ep_old(1,1) - ep_old(0,0)*ep_old(2,2);

   	real64 ep_J2 = ep_old[0]*ep_old[0] + ep_old[1]*ep_old[1] + ep_old[2]*ep_old[2]
		 	   + 3*ep_old[3]*ep_old[3] + 3*ep_old[4]*ep_old[4] + 3*ep_old[5]*ep_old[5]
		 	   - ep_old[1]*ep_old[2] - ep_old[0]*ep_old[1] - ep_old[0]*ep_old[2]; 

		if (ep_J2 > 0.)
		{
			real64 equilivantPlasticStrain = (2./3.)*LvArray::math::sqrt(ep_J2);
			hardening = m_strainHardeningK*( 1.0 - LvArray::math::exp(-1.0*m_strainHardeningN*equilivantPlasticStrain) );
		}
   }

  real64 a1,
         a2,
         a3,
         a4;
  computeLimitParameters( a1,
                          a2,
                          a3,
                          a4,
                          coher_old,
                          hardening );

  int YIELD = computeYieldFunction( I1_trial,
                                    rJ2_trial,
                                    X_old,
                                    Zeta_old,
                                    coher_old,
                                    a1,
                                    a2,
                                    a3,
                                    a4 );

  if (YIELD == -1)
  { // elastic substep
    X_new = X_old;
    Zeta_new = Zeta_old;
    coher_new = coher_old;
    // sigma_new = sigma_trial;
    // ep_new = ep_old;
    LvArray::tensorOps::copy< 6 >( sigma_new, sigma_trial );
    LvArray::tensorOps::copy< 6 >( ep_new, ep_old );
    substepFlag = 0;
    goto successfulSubstep;
  }

  if (YIELD == 1) {  // elastic-plastic or fully-plastic substep
  // (5) Compute non-hardening return to initial yield surface:
  //     [sigma_0,d_e_p,0] = (nonhardeningReturn(sigma_trial,sigma_old,X_old,Zeta_old,K,G)
    real64 I1_0,       // I1 at stress update for non-hardening return
           rJ2_0,      // rJ2 at stress update for non-hardening return
           TOL = 1e-4; // bisection convergence tolerance on eta (if changed, change imax)
    real64 S_0[6],        // S (deviator) at stress update for non-hardening return
           d_ep_0[6],     // increment in plastic strain for non-hardening return
           S_old[6];
    real64 I1_old,
           J2_old,
           rJ2_old,
           evp_old = LvArray::tensorOps::symTrace< 3 >( ep_old ); //ep_old.Trace();
    computeInvariants( sigma_old,
                       S_old,
                       I1_old,
                       J2_old,
                       rJ2_old );

    // returnFlag would be != 0 if there was an error in the nonHardeningReturn call, but
    // there are currently no tests in that function that could detect such an error.
    real64 d_e[6];
    // real64 D_copy[6]; // = D;
 
    // d_e = D; // d_e=D*dt
    // d_e *= dt;
    LvArray::tensorOps::copy< 6 >( d_e, D );
    LvArray::tensorOps::scale< 6 >( d_e, dt);

    returnFlag = nonHardeningReturn( I1_trial,
                                     rJ2_trial,
                                     S_trial,
									                   I1_old,
                                     rJ2_old,
                                     S_old,
									                   d_e,
                                     X_old,
                                     Zeta_old,
                                     coher_old,
                                     hardening,
                                     bulk,
                                     shear,
									                   I1_0,
                                     rJ2_0,
                                     S_0,
                                     d_ep_0 );

	if (returnFlag!=0){
		goto failedSubstep;
	}

	// If there is no porosity (p3=0) and no fluid effects (Kf=0) then the nonhardening
	// return will be the solution
	if ( (m_p3 == 0.0)&&(Kf==0.0) ){
        Zeta_new = Zeta_old,
		X_new = X_old;

		// ep_new = ep_old; // ep_old + d_ep_0;
		// ep_new += d_ep_0;
    LvArray::tensorOps::copy< 6 >( ep_new, ep_old );
    LvArray::tensorOps::add< 6 >( ep_new, d_ep_0 );

		// sigma_new = Identity;  // one_third*I1_0*Identity + S_0;
		// sigma_new *= one_third*I1_0;
		// sigma_new += S_0;
    LvArray::tensorOps::copy< 6 >( sigma_new, identity );
    LvArray::tensorOps::scale< 6 >( sigma_new, one_third*I1_0 );
    LvArray::tensorOps::add< 6 >( sigma_new, S_0 );

		goto successfulSubstep;
	}

    real64 d_evp_0 = LvArray::tensorOps::symTrace<3>( d_ep_0 ); //d_ep_0.Trace();

  // (6) Iterate to solve for plastic volumetric strain consistent with the updated
  //     values for the cap (X) and isotropic backstress (Zeta).  Use a bisection method
  //     based on the multiplier eta,  where  0<eta<1

    real64 eta_out = 1.0,
           eta_in = 0.0,
           eta_mid,
           d_evp;
    int i = 0,
        imax = 93;  // imax = LvArray::math::ceil(-10.0*LvArray::math::log(TOL)); // Update this if TOL changes

    real64 dZetadevp = computedZetadevp(fluid_pressure_initial,Km,Kf,ev0,phi_i,Zeta_old,evp_old);

  // (7) Update Internal State Variables based on Last Non-Hardening Return:
  //
  updateISV:
    i++;
    eta_mid   = 0.5 * ( eta_out + eta_in );
    d_evp     = eta_mid * d_evp_0;

    // Update X exactly
    real64 evp_new = evp_old + d_evp;

    X_new = computeX( evp_new,
                      phi_i,
                      Km,
                      Kf,
                      C1,
                      ev0,
                      fluid_pressure_initial );

    // update the damage variable.  this is untested, and only active if Gf > 0
    computeCoher( lch,
                  I1_trial,
                  I1_0,
                  d_evp,
                  m_fractureEnergyReleaseRate,
                  coher_old,
                  coher_new );

    // Update zeta. min() eliminates tensile fluid pressure from explicit integration error
    Zeta_new = LvArray::math::min( Zeta_old + dZetadevp * d_evp, 0.0 );

    // (8) Check if the updated yield surface encloses trial stres.  If it does, there is too much
    //     plastic strain for this iteration, so we adjust the bisection parameters and recompute
    //     the state variable update.

    if( computeYieldFunction( I1_trial,
                              rJ2_trial,
                              X_new,
                              Zeta_new,
                              coher_new,
                              a1,
                              a2,
                              a3,
                              a4 ) !=1 )
    {
      eta_out = eta_mid;
      if( i >= imax ){
        // solution failed to converge within the allowable iterations, which means
        // the solution requires a plastic strain that is less than TOL*d_evp_0
        // In this case we are near the zero porosity limit, so the response should
        // be that of no porosity. By setting eta_out=eta_in, the next step will
        // converge with the cap position of the previous iteration.  In this case,
        // we set evp=-p3 (which corresponds to X=1e12*p0) so subsequent compressive
        // loading will respond as though there is no porosity.  If there is dilatation
        // in subsequent loading, the porosity will be recovered.
        eta_out = eta_in;
      }
      goto updateISV;
    }

    // (9) Recompute the elastic properties based on the midpoint of the updated step:
    //     [K,G] = computeElasticProperties( (sigma_old+sigma_new)/2,ep_old+d_ep/2 )
    //     and compute return to updated surface.
    //    MH! change this when elastic-plastic coupling is used.
    //    R2Tensor sigma_new = ...,
    //            ep_new    = ...,
    //    computeElasticProperties((sigma_old+sigma_new)/2,(ep_old+ep_new)/2,bulk,shear);
    real64 I1_new,
           rJ2_new,
           d_evp_new;
    real64 S_new[6],
           d_ep_new[6];
    returnFlag = nonHardeningReturn( I1_trial,
                                     rJ2_trial,
                                     S_trial,
                                     I1_old,
                                     rJ2_old,
                                     S_old,
                                     d_e,
                                     X_new,
                                     Zeta_new,
                                     coher_new,
                                     hardening,
                                     bulk,
                                     shear,
                                     I1_new,
                                     rJ2_new,
                                     S_new,
                                     d_ep_new );

    if( returnFlag !=0 )
    {
      goto failedSubstep;
    }

    // (10) Check whether the isotropic component of the return has changed sign, as this
    //      would indicate that the cap apex has moved past the trial stress, indicating
    //      too much plastic strain in the return.

    //if(fabs(I1_trial - I1_new)>(m_mat_geo_B0[p_mat]*TOL) && Sign(I1_trial - I1_new)!=Sign(I1_trial - I1_0)){
    real64 sgnI1tmn = ( (I1_trial - I1_new) < 0.0 ) ? ( -1.0 ) : ( 1.0 );
    real64 sgnI1tm0 = ( (I1_trial - I1_0) < 0.0 ) ? ( -1.0 ) : ( 1.0 );
    if( abs( sgnI1tmn - sgnI1tm0 ) > 1e-12 ){
      eta_out = eta_mid;
      if( i >= imax ){
        // solution failed to converge within the allowable iterations, which means
        // the solution requires a plastic strain that is less than TOL*d_evp_0
        // In this case we are near the zero porosity limit, so the response should
        // be that of no porosity. By setting eta_out=eta_in, the next step will
        // converge with the cap position of the previous iteration.  In this case,
        // we set evp=-p3 (which corresponds to X=1e12*p0) so subsequent compressive
        // loading will respond as though there is no porosity.  If there is dilatation
        // in subsequent loading, the porosity will be recovered.
        eta_out = eta_in;
      }
      goto updateISV;
    }

    // Compare magnitude of plastic strain with prior update
    d_evp_new = LvArray::tensorOps::symTrace< 3 >( d_ep_new ); //d_ep_new.Trace();   // Increment in vol. plastic strain for return to new surface
    // ep_new = ep_old; //ep_old + d_ep_new;
    // ep_new += d_ep_new;
    LvArray::tensorOps::copy< 6 >( ep_new, ep_old );
    LvArray::tensorOps::add< 6 >( ep_new, d_ep_new );

    // Check for convergence
    if( fabs(eta_out-eta_in) < TOL ){ // Solution is converged
      // sigma_new = Identity; // one_third*I1_new*Identity + S_new;
      // sigma_new *= one_third*I1_new;
      // sigma_new += S_new;
      LvArray::tensorOps::copy< 6 >( sigma_new, identity );
      LvArray::tensorOps::scale< 6 >( sigma_new, one_third*I1_new );
      LvArray::tensorOps::add< 6 >( sigma_new, S_new );

      // If out of range, scale back isotropic plastic strain.
      if( LvArray::tensorOps::symTrace< 3 >( ep_new ) < -m_p3 ) //ep_new.Trace()<-p3){
      {
        d_evp_new = -m_p3 - LvArray::tensorOps::symTrace< 3 >( ep_old ); //ep_old.Trace();

        real64 d_ep_new_iso[6],
		           d_ep_new_dev[6];

        // d_ep_new_iso = Identity; // one_third*d_ep_new.Trace()*Identity;
        // d_ep_new_iso *= one_third*d_ep_new.Trace();
        LvArray::tensorOps::copy< 6 >( d_ep_new_iso, identity );
        LvArray::tensorOps::scale< 6 >( d_ep_new_iso, one_third*LvArray::tensorOps::symTrace< 3 >( d_ep_new ) );

        // d_ep_new_dev = d_ep_new; //d_ep_new - d_ep_new_iso;
        // d_ep_new_dev -= d_ep_new_iso;
        LvArray::tensorOps::copy< 6 >( d_ep_new_dev, d_ep_new );
        LvArray::tensorOps::subtract< 6 >( d_ep_new_dev, d_ep_new_iso );

        // ep_new = Identity; // ep_old + d_ep_new_dev + one_third*d_evp_new*Identity;
        // ep_new *= one_third*d_evp_new;
        // ep_new += ep_old;
        // ep_new += d_ep_new_dev;
        LvArray::tensorOps::copy< 6 >( ep_new, identity );
        LvArray::tensorOps::scale< 6 >( ep_new, one_third*d_evp_new );
        LvArray::tensorOps::add< 6 >( ep_new, ep_old );
        LvArray::tensorOps::add< 6 >( ep_new, d_ep_new_dev );
      }

      // Update X exactly
      real64 evp;
      evp = LvArray::tensorOps::symTrace< 3 >( ep_new ); //ep_new.Trace();

      X_new = computeX(	evp,
        					      phi_i, 
                        Km, 
                        Kf, 
                        C1, 
                        ev0,
        					      fluid_pressure_initial );

      // Update zeta. min() eliminates tensile fluid pressure from explicit integration error
      Zeta_new = LvArray::math::min(Zeta_old + dZetadevp*d_evp_new,0.0);

      goto successfulSubstep;
    }
    if( i >= imax ){
      // Solution failed to converge but not because of too much plastic strain
      // (which would have been caught by the checks above).  In this case we
      // go to the failed substep return, which will trigger subcycling and
      // particle deletion (if subcycling doesn't work).
      //
      // This code was never reached in testing, but is here to catch
      // unforseen errors.
      goto failedSubstep;
    }

// (11) Compare magnitude of the volumetric plastic strain and bisect on eta
//
    if( fabs(d_evp_new) > eta_mid*fabs(d_evp_0) ){
      eta_in = eta_mid;
    }
    else {
      eta_out = eta_mid;
    }
    goto updateISV;
  }

// (12) Compute softening, collapse yield surface, and compute nonhardening return
// to collapse yield surface.
//
// This is not rigorous, since the treatment of softening is not included as a
// hardenining mechanism in the computation of the increment in plastic strain.
//
// Figure out where to put this:
// coher = max(coher - d_ep_new.Norm()/m_mat_geo_failure_strain[p_mat]);
//  //
//  returnFlag = nonHardeningReturn(I1_trial,rJ2_trial,S_trial,
//								  I1_old,rJ2_old,S_old,
//								  d_e,X_new,Zeta_new,coher,bulk,shear,
//								  I1_new,rJ2_new,S_new,d_ep_new);
//  if (returnFlag!=0){
//	  goto failedSubstep;



// (13) Return updated values for successful/unsuccessful steps
//
successfulSubstep:
  substepFlag = 0;
  return substepFlag;

failedSubstep:
  // sigma_new = sigma_old;
  // ep_new = ep_old;
  LvArray::tensorOps::copy< 6 >( sigma_new, sigma_old );
  LvArray::tensorOps::copy< 6 >( ep_new, ep_old );
  X_new = X_old;
  Zeta_new = Zeta_old;
  substepFlag = 1;
  return substepFlag;
}

// Compute nonhardening return from trial stress to some yield surface
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::nonHardeningReturn( const real64 & I1_trial,              // Trial Stress
                                             const real64 & rJ2_trial,          
                                             real64 const ( & S_trial )[6],          
                                             const real64 & I1_old,                // Stress at start of subtep
                                             const real64 & GEOS_UNUSED_PARAM( rJ2_old ),          
                                             real64 const ( & S_old )[6],          
                                             real64 const ( & d_e )[6],            // increment in total strain
                                             const real64 & X,                     // cap position
                                             const real64 & Zeta,                  // isotropic bacstress
                                             const real64 & coher,    
                                             const real64 & hardening,          
                                             const real64 & bulk,                  // elastic bulk modulus
                                             const real64 & shear,                 // elastic shear modulus
                                             real64 & I1_new,                      // New stress state on yield surface
                                             real64 & rJ2_new,
                                             real64 ( & S_new )[6],
                                             real64 ( & d_ep_new )[6]            // increment in plastic strain for return
) const
{
  // Computes a non-hardening return to the yield surface in the meridional profile
  // (constant Lode angle) based on the current values of the internal state variables
  // and elastic properties.  Returns the updated stress and  the increment in plastic
  // strain corresponding to this return.
  //
  // NOTE: all values of r and z in this function are transformed!

	real64 one_sqrt_three = 0.5773502691896258;
	real64 sqrt_two = 1.414213562373095;
	real64 one_third = 0.3333333333333333;
	real64 sqrt_three = 1.732050807568877;
	real64 one_sixth = 0.1666666666666667;
	real64 one_ninth = 0.1111111111111111;

	real64 identity[6] = {1.0,1.0,1.0,0.0,0.0,0.0};

  const int nmax = 19;  // If this is changed, more entries may need to be added to sinV cosV.
  int n = 0,
	    returnFlag = 0,   // error flag = 0 for successful return.
      interior;

  // (1) Define an interior point, (nominally I1_0 = Zeta, also, J2_0 = 0 but no need to  create this variable.)
  real64 I1_0,
		     I1trialMinusZeta = I1_trial-Zeta;

  // It may be better to use an interior point at the center of the yield surface, rather than at zeta, in particular
  // when PEAKI1=0.  Picking the midpoint between PEAKI1 and X would be problematic when the user has specified
  // some no porosity condition (e.g. p0=-1e99)
  if( I1trialMinusZeta>= coher * m_peakT1 ) // Trial is past vertex
  { 
	  real64 lTrial = LvArray::math::sqrt(I1trialMinusZeta * I1trialMinusZeta + rJ2_trial * rJ2_trial),
			     lYield = 0.5 * (coher * m_peakT1 - X);
	  I1_0 = Zeta + coher * m_peakT1 - LvArray::math::min(lTrial, lYield);
  }
  else if( (I1trialMinusZeta < coher * m_peakT1) && (I1trialMinusZeta > X) ){ // Trial is above yield surface
	  I1_0 = I1_trial;
  }
  else if( I1trialMinusZeta <= X ) // Trial is past X, use yield midpoint as interior point
  {
	  I1_0 = Zeta + 0.5 * (coher * m_peakT1 + X);
  }
  else
  { // Shouldn't get here
	  I1_0 = Zeta;
  }

  // (2) Transform the trial and interior points as follows where beta defines the degree
  //  of non-associativity.
  // multiplier to compute Lode R to sqrt(J2)
  real64 rJ2_to_r = sqrt_two * m_beta * LvArray::math::sqrt(1.5 * bulk / shear);
  // multiplier to compute sqrt(J2) to Lode R
  real64 r_to_rJ2 = 1.0 / rJ2_to_r;
  real64 r_trial = rJ2_to_r * rJ2_trial,
         z_trial = I1_trial * one_sqrt_three,
         z_test,
         r_test,
         r_0 = 0.0,
         z_0 = I1_0 * one_sqrt_three;

  // Lookup tables for computing the sin() and cos() of th rotation angle.
  real64 sinV[]={0.7071067811865475,-0.5,0.3420201433256687,-0.2306158707424402,0.1545187928078405,
                  -0.1032426220806015,0.06889665647555759,-0.04595133277786571,0.03064021661344469,
                  -0.02042858745187096,0.01361958465478159,-0.009079879062402308,0.006053298918749807,
                  -0.004035546304539714,0.002690368259933135,-0.001793580042002626,0.001195720384163988,
                  -0.0007971470283055577,0.0005314313834717263,-0.00035428759824575,0.0002361917349088998};
  real64 cosV[]={0.7071067811865475,0.8660254037844386,0.9396926207859084,0.9730448705798238,
                  0.987989849476809,0.9946562024066014,0.9976238022052647,0.9989436796015769,
                  0.9995304783376449,0.9997913146325693,0.999907249155556,0.9999587770484402,
                  0.9999816786182636,0.999991857149859,0.9999963809527642,0.9999983915340229,
                  0.9999992851261259,0.9999996822782572,0.9999998587903324,0.9999999372401469,
                  0.9999999721067318};
  real64 sinTheta = sinV[0],
         cosTheta = cosV[0];

  // Compute the a1,a2,a3,a4 parameters from FSLOPE,YSLOPE,STREN and PEAKI1,
  // which are perturbed by variability according to coher.  These are then
  // passed down to the computeYieldFunction, to avoid the expense of computing a3
  real64 a1,
         a2,
         a3,
         a4;
  computeLimitParameters( a1,
                          a2,
                          a3,
                          a4,
                          coher,
                          hardening );

  // (3) Perform Bisection between in transformed space, to find the new point on the
  //  yield surface: [znew,rnew] = transformedBisection(z0,r0,z_trial,r_trial,X,Zeta,K,G)
  int k = 0;
  while ( (n < nmax)&&(k < 10*nmax) )
  {
    // transformed bisection to find a new interior point, just inside the boundary of the
    // yield surface.  This function overwrites the inputs for z_0 and r_0
    //  [z_0,r_0] = transformedBisection(z_0,r_0,z_trial,r_trial,X_Zeta,bulk,shear)
    transformedBisection( z_0,
                          r_0,
                          z_trial,
                          r_trial,
                          X,
                          Zeta,
                          coher,
                          a1,
                          a2,
                          a3,
                          a4,
                          r_to_rJ2 );

    // (4) Perform a rotation of {z_new,r_new} about {z_trial,r_trial} until a new interior point
    // is found, set this as {z0,r0}
    interior = 0;
    n = LvArray::math::max(n-4,0);  //
    // (5) Test for convergence:
    while ( (interior==0)&&(n < nmax) )
    {
		  k++;
      
      // To avoid the cost of computing pow() to get theta, and then sin(), cos(),
      // we use a lookup table defined above by sinV and cosV.
      //
      // theta = pi_fourth*Pow(-two_third,n);
      // z_test = z_trial + cos(theta)*(z_0-z_trial) - sin(theta)*(r_0-r_trial);
      // r_test = r_trial + sin(theta)*(z_0-z_trial) + cos(theta)*(r_0-r_trial);
      sinTheta = sinV[n];
      cosTheta = cosV[n];
      z_test = z_trial + cosTheta*(z_0-z_trial) - sinTheta*(r_0-r_trial);
      r_test = r_trial + sinTheta*(z_0-z_trial) + cosTheta*(r_0-r_trial);

      if ( transformedYieldFunction( z_test,
                                     r_test,
                                     X,
                                     Zeta,
                                     coher,
                                     a1,
                                     a2,
                                     a3,
                                     a4,
                                     r_to_rJ2) == -1 ) 
      { // new interior point
        interior = 1;
        z_0 = z_test;
        r_0 = r_test;
      }
      else { n++; }
    }
  }

  if (k>=10*nmax)
  {
	  returnFlag = 1;
  }

// (6) Solution Converged, Compute Untransformed Updated Stress:
  I1_new = sqrt_three*z_0;
  rJ2_new = r_to_rJ2*r_0;

  LvArray::tensorOps::copy< 6 >( S_new, S_trial );
  if ( rJ2_trial != 0.0 )
  {
	  // S_new = S_trial; //S_trial*rJ2_new/rJ2_trial;
	  // S_new *= rJ2_new/rJ2_trial;
    LvArray::tensorOps::scale< 6 >( S_new, rJ2_new / rJ2_trial );
  }
  // else
  // {
	//   S_new = S_trial;
  // }

  real64 sigma_new[6],
         sigma_old[6],
         d_sigma[6];
  // sigma_new = Identity; // I1_new*one_third*Identity + S_new;
  // sigma_new *= I1_new*one_third;
  // sigma_new += S_new;
  LvArray::tensorOps::copy< 6 >( sigma_new, identity );
  LvArray::tensorOps::scale< 6 >( sigma_new, I1_new*one_third );
  LvArray::tensorOps::add< 6 >( sigma_new, S_new );

  // sigma_old = Identity; // I1_old*one_third*Identity + S_old;
  // sigma_old *= I1_old*one_third;
  // sigma_old += S_old;
  LvArray::tensorOps::copy< 6 >( sigma_old, identity );
  LvArray::tensorOps::scale< 6 >( sigma_old, I1_old*one_third );
  LvArray::tensorOps::add< 6 >( sigma_old, S_old );

  // d_sigma = sigma_new; // sigma_new - sigma_old;
  // d_sigma -= sigma_old;
  LvArray::tensorOps::copy< 6 >( d_sigma, sigma_new );
  LvArray::tensorOps::subtract< 6 >( d_sigma, sigma_old );

  // (7) Compute increment in plastic strain for return:
  //  d_ep0 = d_e - [C]^-1:(sigma_new-sigma_old)
  real64 d_ee[6], 
         d_ee_2[6];

  // d_ee = Identity; // 0.5*d_sigma/shear + (one_ninth/bulk - one_sixth/shear)*d_sigma.Trace()*Identity;
  // d_ee *= (one_ninth/bulk - one_sixth/shear)*d_sigma.Trace();
  LvArray::tensorOps::copy< 6 >( d_ee, identity );
  LvArray::tensorOps::scale< 6 >( d_ee, (one_ninth / bulk - one_sixth / shear) * LvArray::tensorOps::symTrace< 3 >( d_sigma ) );

  // d_ee_2 = d_sigma;
  // d_ee_2 *= (0.5/shear);
  // d_ee += d_ee_2;
  LvArray::tensorOps::copy< 6 >( d_ee_2, d_sigma );
  LvArray::tensorOps::scale< 6 >( d_ee_2, 0.5 / shear );
  LvArray::tensorOps::add< 6 >( d_ee, d_ee_2 );

  // d_ep_new = d_e; // d_e - d_ee;
  // d_ep_new -= d_ee;
  LvArray::tensorOps::copy< 6 >( d_ep_new, d_e );
  LvArray::tensorOps::subtract< 6 >( d_ep_new, d_ee );

  return returnFlag;
}

// Computes bisection between two points in transformed space
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::transformedBisection(real64 & z_0,
                                               real64 & r_0,
                                               const real64 & z_trial,
                                               const real64 & r_trial,
                                               const real64 & X,
                                               const real64 & Zeta,
									                             const real64 & coher,
									                             const real64 & a1,
									                             const real64 & a2,
									                             const real64 & a3,
									                             const real64 & a4,
                                               const real64 & r_to_rJ2
) const
{
  // Computes a bisection in transformed stress space between point sigma_0 (interior to the
  // yield surface) and sigma_trial (exterior to the yield surface).  Returns this new point,
  // which will be just outside the yield surface, overwriting the input arguments for
  // z_0 and r_0.

  // After the first iteration of the non-hardening return, the subsequent bisections will likely
  // converge with eta << 1.  It may be faster to put in some logic to try to start bisection
  // with tighter bounds, and only expand them to 0<eta<1 if the first eta_mid is too large.


  // (1) initialize bisection
  real64 eta_out = 1.0,
         eta_in  = 0.0,
         eta_mid,
         TOL = 1.0e-6,
         r_test,
         z_test;

  // (2) Test for convergence
  while (eta_out-eta_in > TOL){

    // (3) Transformed test point
    eta_mid = 0.5*(eta_out+eta_in);
    z_test = z_0 + eta_mid*(z_trial-z_0);
    r_test = r_0 + eta_mid*(r_trial-r_0);
    // (4) Check if test point is within the yield surface:
    if ( transformedYieldFunction( z_test,
                                   r_test,
                                   X,
                                   Zeta,
                                   coher,
                                   a1,
                                   a2,
                                   a3,
                                   a4,
                                   r_to_rJ2 ) !=1 )
    {
    	eta_in = eta_mid;
    }
    else
    {
    	eta_out = eta_mid;
    }
  }
  
  // (5) Converged, return {z_new,r_new}={z_test,r_test}
  z_0 = z_0 + eta_out*(z_trial-z_0); //z_0 = z_test;
  r_0 = r_0 + eta_out*(r_trial-r_0); //r_0 = r_test;
}

// computeTransformedYieldFunction from transformed inputs
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::transformedYieldFunction( const real64 & z,
                                                   const real64 & r,
                                                   const real64 & X,
                                                   const real64 & Zeta,
										                               const real64 & coher,
										                               const real64 & a1,
										                               const real64 & a2,
										                               const real64 & a3,
										                               const real64 & a4,
                                                   const real64 & r_to_rJ2
) const
{
  // Evaluate the yield criteria and return:
  //  -1: elastic
  //   0: on yield surface within tolerance
  //   1: plastic
	real64 sqrt_three = 1.732050807568877;
  
  // Untransformed values:
  real64 I1  = sqrt_three * z,
		     rJ2 = r_to_rJ2 * r;
  int    YIELD = computeYieldFunction( I1,
                                       rJ2,
                                       X,
                                       Zeta,
                                       coher,
                                       a1,
                                       a2,
                                       a3,
                                       a4 );
  return YIELD;
}

// computeYieldFunction from untransformed inputs
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int GeomechanicsUpdates::computeYieldFunction( const real64 & I1,
                                               const real64 & rJ2,
                                               const real64 & X,
                                               const real64 & Zeta,
									                             const real64 & GEOS_UNUSED_PARAM( coher ),
									                             const real64 & a1,
									                             const real64 & a2,
									                             const real64 & a3,
									                             const real64 & a4
) const    
{
	// Evaluate the yield criteria and return:
	//  -1: elastic
	//   0: on yield surface within tolerance (not used)
	//   1: plastic
	//
	//                        *** Developer Note ***
	// THIS FUNCTION IS DEEP WITHIN A NESTED LOOP AND IS CALLED THOUSANDS
	// OF TIMES PER TIMESTEP.  EVERYTHING IN THIS FUNCTION SHOULD BE
	// OPTIMIZED FOR SPEED.
	//
	int YIELD = -1;
	real64 I1mZ = I1 - Zeta;    // Shifted stress to evalue yield criteria

	// --------------------------------------------------------------------
	// *** SHEAR LIMIT FUNCTION (Ff) ***
	// --------------------------------------------------------------------
	// Read input parameters to specify strength model
	real64  Ff;
	Ff = a1 - a3*LvArray::math::exp(a2*I1mZ) - a4*I1mZ;

	// --------------------------------------------------------------------
	// *** Branch Point (Kappa) ***
	// --------------------------------------------------------------------
	real64  Kappa  = m_peakT1-m_cr*(m_peakT1-X); // Branch Point

	// --------------------------------------------------------------------
	// *** COMPOSITE YIELD FUNCTION ***
	// --------------------------------------------------------------------
	// Evaluate Composite Yield Function F(I1) = Ff(I1)*fc(I1) in each region.
	// The elseif statements have nested if statements, which is not equivalent
	// to them having a single elseif(A&&B&&C)
	if( I1mZ < X )
	{//---------------------------------------------------(I1<X)
		YIELD = 1;
	}
	else if(( I1mZ < Kappa )&&( I1mZ >= X )) 
  {// ---------------(X<I1<kappa)
		// p3 is the maximum achievable volumetric plastic strain in compresson
		// so if a value of 0 has been specified this indicates the user
		// wishes to run without porosity, and no cap function is used, i.e. fc=1

		// **Elliptical Cap Function: (fc)**
		// fc = LvArray::math::sqrt(1.0 - Pow((Kappa-I1mZ)/(Kappa-X)),2.0);
		// faster version: fc2 = fc^2
		real64 fc2 = 1.0 - ((Kappa-I1mZ)/(Kappa-X))*((Kappa-I1mZ)/(Kappa-X));
		if(rJ2*rJ2 > Ff*Ff*fc2 )
		{
			YIELD = 1;
		}
	}
	else if(( I1mZ <= m_peakT1 )&&( I1mZ >= Kappa ))
  { // -----(kappa<I1<PEAKI1)
		if(rJ2 > Ff) {
			YIELD = 1;
		}
	}
	else if( I1mZ > m_peakT1 )
	{// --------------------------------(peakI1<I1)
    YIELD = 1;
	};

  return YIELD;
} 

// Compute (dZeta/devp) Zeta and vol. plastic strain
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GeomechanicsUpdates::computedZetadevp( real64 const & fluid_pressure_initial,
                                              real64 const & Km,
                                              real64 const & Kf,
                                              real64 const & ev0,
                                              real64 const & phi_i,
                                              real64 Zeta,
                                              real64 evp
) const
{
  // Computes the partial derivative of the trace of the
  // isotropic backstress (Zeta) with respect to volumetric
  // plastic strain (evp).
  real64 dZetadevp = 0.0;           // Evolution rate of isotropic backstress

  if (evp <= ev0 && Kf != 0.0) { // .................................... Fluid effects are active
    real64 pfi = fluid_pressure_initial; // initial fluid pressure

    // This is an expensive calculation, but fasterexp() seemed to cause errors.
    dZetadevp = (3.0*LvArray::math::exp(evp)*Kf*Km)/(LvArray::math::exp(evp)*(Kf + Km)
                                      + LvArray::math::exp(Zeta/(3.0*Km))*Km*(-1.0 + phi_i)
                                      - LvArray::math::exp((3.0*pfi + Zeta)/(3.0*Kf))*Kf*phi_i);
  }
  return dZetadevp;
} 

// Compute (dZeta/devp) Zeta and vol. plastic strain
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeomechanicsUpdates::computeLimitParameters( real64 & a1,
		                                              real64 & a2,
		                                              real64 & a3,
		                                              real64 & a4,
		                                              const real64 & GEOS_UNUSED_PARAM( coher ),
                                                  const real64 & hardening
) const 
{ // Value of I1 at strength=0 (Perturbed by variability)
  // The shear limit surface is defined in terms of the a1,a2,a3,a4 parameters, but
  // the user inputs are the more intuitive set of FSLOPE. YSLOPE, STREN, and PEAKI1.

  // This routine computes the a_i parameters from the user inputs.  The code was
  // originally written by R.M. Brannon, with modifications by M.S. Swan.
  // harden peakI1 and stren together to not change slope:
	real64 stren_h = m_stren + hardening;
 	real64 peakT1_h;
  if (m_fSlope > 0.0)
  {
     peakT1_h = m_peakT1 + hardening/m_fSlope;
  } 
  else
  {
    peakT1_h = m_peakT1;
  }

  if (m_fSlope > 0.0 && m_peakT1 >= 0.0 && m_stren == 0.0 && m_ySlope == 0.0)
  {// ----------------------------------------------Linear Drucker-Prager
    a1 = m_peakT1 * m_fSlope;
    a2 = 0.0;
    a3 = 0.0;
    a4 = m_fSlope;
  }
  else if (m_fSlope == 0.0 && m_peakT1 == 0.0 && m_stren > 0.0 && m_ySlope == 0.0)
  { // ------------------------------------------------------- Von Mises
    a1 = stren_h;
    a2 = 0.0;
    a3 = 0.0;
    a4 = 0.0;
  }
  else if (m_fSlope > 0.0 && m_ySlope  == 0.0 && m_stren > 0.0 && m_peakT1 == 0.0)
  { // ------------------------------------------------------- 0 PEAKI1 to vonMises
    a1 = stren_h;
    a2 = m_fSlope / stren_h;
    a3 = stren_h;
    a4 = 0.0;
  }
  else if (m_fSlope > m_ySlope && m_ySlope > 0.0 && m_stren > m_ySlope*m_peakT1 && m_peakT1 >= 0.0)
  { // ------------------------------------------------------- Nonlinear Drucker-Prager
    a1 = stren_h;
    a2 = (m_fSlope-m_ySlope )/(stren_h-m_ySlope *peakT1_h);
    a3 = (stren_h-m_ySlope *peakT1_h)*LvArray::math::exp(-a2*peakT1_h);
    a4 = m_ySlope ;
  }
  else
  {
	  GEOS_LOG_RANK("Bad limit surface parameters.");
  }
}


/**
 * @class Geomechanics
 *
 * Geomechanics material model.
 */
class Geomechanics : public SolidBase
{
public:

  /// @typedef Alias for GeomechanicsUpdates
  using KernelWrapper = GeomechanicsUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Geomechanics( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Geomechanics() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Geomechanics";

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
    /// string/key for tangent elastic bulk modulus parameter 0 
    static constexpr char const * b0String() { return "b0"; }

    /// string/key for tangent elastic bulk modulus parameter 1 
    static constexpr char const * b1String() { return "b1"; }

    /// string/key for tangent elastic bulk modulus parameter 2 
    static constexpr char const * b2String() { return "b2"; }

    /// string/key for tangent elastic bulk modulus parameter 3 
    static constexpr char const * b3String() { return "b3"; }

    /// string/key for tangent elastic bulk modulus parameter 4 
    static constexpr char const * b4String() { return "b4"; }

    /// string/key for tangent elastic shear modulus parameter 0 
    static constexpr char const * g0String() { return "g0"; }

    /// string/key for tangent elastic shear modulus parameter 1 
    static constexpr char const * g1String() { return "g1"; }

    /// string/key for tangent elastic shear modulus parameter 2 
    static constexpr char const * g2String() { return "g2"; }

    /// string/key for tangent elastic shear modulus parameter 3 
    static constexpr char const * g3String() { return "g3"; }

    /// string/key for tangent elastic shear modulus parameter 4 
    static constexpr char const * g4String() { return "g4"; }

    /// string/key for crush curve parameter 0 
    static constexpr char const * p0String() { return "p0"; }

    /// string/key for crush curve parameter 1 
    static constexpr char const * p1String() { return "p1"; }

    /// string/key for crush curve parameter 2 
    static constexpr char const * p2String() { return "p2"; }

    /// string/key for crush curve parameter 3 
    static constexpr char const * p3String() { return "p3"; }

    /// string/key for crush curve parameter 4 
    static constexpr char const * p4String() { return "p4"; }

    /// string/key for cap shape parameter
    static constexpr char const * crString() { return "cr"; }

    /// string/key for fluid bulk modulus
    static constexpr char const * fluidBulkModulusString() { return "fluidBulkModulus"; }

    /// string/key for fluid initial pressure
    static constexpr char const * fluidInitialPressureString() { return "fluidInitialPressure"; }

    /// string/key for rate dependence parameter 1
    static constexpr char const * t1RateDependenceString() { return "t1RateDependence"; }

    /// string/key for rate dependence parameter 2
    static constexpr char const * t2RateDependenceString() { return "t2RateDependence"; }

    /// string/key for peak t1 shear limit parameter
    static constexpr char const * peakT1String() { return "peakT1"; }

    /// string/key for F slope shear limit parameter
    static constexpr char const * fSlopeString() { return "fSlope"; }
         
    /// string/key for stren shear limit parameter
    static constexpr char const * strenString() { return "stren"; }

    /// string/key for Y slope shear limit parameter
    static constexpr char const * ySlopeString() { return "ySlope"; }

    /// string/key for nonassociativity parameter
    static constexpr char const * betaString() { return "beta"; }

    /// string/key for fracture energy release rate
    static constexpr char const * fractureEnergyReleaseRateString() { return "fractureEnergyReleaseRate"; }

    /// string/key for creep flag
    static constexpr char const * creepString() { return "creep"; }

    /// string/key for creep C0 parameter
    static constexpr char const * creepC0String() { return "creepC0"; }

    /// string/key for creep C1 parameter
    static constexpr char const * creepC1String() { return "creepC1"; }

    /// string/key for creep A parameter
    static constexpr char const * creepAString() { return "creepA"; }

    /// string/key for creep B parameter
    static constexpr char const * creepBString() { return "creepB"; }

    /// string/key for creep C parameter
    static constexpr char const * creepCString() { return "creepC"; }

    /// string/key for strain-hardening N parameter
    static constexpr char const * strainHardeningNString() { return "strainHardeningN"; }

    /// string/key for strain-hardening K parameter
    static constexpr char const * strainHardeningKString() { return "strainHardeningK"; }

    //string/key for element/particle bulk modulus value
    static constexpr char const * bulkModulusString() { return "bulkModulus"; }

    //string/key for element/particle shear modulus value
    static constexpr char const * shearModulusString() { return "shearModulus"; }

    //string/key for element/particle velocityGradient value
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for quadrature point plasticStrain value 
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point porosity value 
    static constexpr char const * porosityString() { return "porosity"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }
  };

  /**
   * @brief Create a instantiation of the GeomechanicsUpdate class that refers to the data in this.
   * @return An instantiation of GeomechanicsUpdate.
   */
  GeomechanicsUpdates createKernelUpdates() const
  {
    return GeomechanicsUpdates( m_b0,
                                m_b1,
                                m_b2,
                                m_b3,
                                m_b4,
                                m_g0,
                                m_g1,
                                m_g2,
                                m_g3,
                                m_g4,
                                m_p0,
                                m_p1,
                                m_p2,
                                m_p3,
                                m_p4,
                                m_peakT1,
                                m_fSlope,
                                m_stren,
                                m_ySlope,
                                m_beta,
                                m_t1RateDependence,
                                m_t2RateDependence,
                                m_fractureEnergyReleaseRate,
                                m_cr,
                                m_fluidBulkModulus,
                                m_fluidInitialPressure,
                                m_creep,
                                m_creepC0,
                                m_creepC1,
                                m_creepA,
                                m_creepB,
                                m_creepC,
                                m_strainHardeningN,
                                m_strainHardeningK,
                                m_bulkModulus,
                                m_shearModulus,
                                m_velocityGradient,
                                m_plasticStrain,
                                m_porosity,
                                m_damage,
                                m_lengthScale,
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
                          m_b0,
                          m_b1,
                          m_b2,
                          m_b3,
                          m_b4,
                          m_g0,
                          m_g1,
                          m_g2,
                          m_g3,
                          m_g4,
                          m_p0,
                          m_p1,
                          m_p2,
                          m_p3,
                          m_p4,
                          m_peakT1,
                          m_fSlope,
                          m_stren,
                          m_ySlope,
                          m_beta,
                          m_t1RateDependence,
                          m_t2RateDependence,
                          m_fractureEnergyReleaseRate,
                          m_cr,
                          m_fluidBulkModulus,
                          m_fluidInitialPressure,
                          m_creep,
                          m_creepC0,
                          m_creepC1,
                          m_creepA,
                          m_creepB,
                          m_creepC,
                          m_strainHardeningN,
                          m_strainHardeningK,
                          m_bulkModulus,
                          m_shearModulus,
                          m_velocityGradient,
                          m_plasticStrain,
                          m_porosity,
                          m_damage,
                          m_lengthScale,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }

// CC: commmented out in case I need to use this for template of setters and getters for scalar values
//   /**
//    * @brief Getter for default transverse Young's modulus
//    * @return The value of the default transverse Young's modulus.
//    */
//   real64 getDefaultYoungModulusTransverse() const
//   {
//     return m_defaultYoungModulusTransverse;
//   }

//   /**
//    * @brief Setter for the default transverse Young's modulus.
//    * @param[in] input New value for the default transverse Young's modulus
//    */
//   void setDefaultYoungModulusTransverse( real64 const input )
//   {
//     m_defaultYoungModulusTransverse = input;
//   }

//   /**
//    * @brief Getter for default axial Young's modulus
//    * @return The value of the default axial Young's modulus.
//    */
//   real64 getDefaultYoungModulusAxial() const
//   {
//     return m_defaultYoungModulusAxial;
//   }

  /**
   * @brief Accessor for bulk modulus
   * @return A const reference to arrayView1d<real64> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const bulkModulus() { return m_bulkModulus; }

  /**
   * @brief Const accessor for bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         bulk modulus (at every element).
   */
  arrayView1d< real64 const > const bulkModulus() const { return m_bulkModulus; }

  /**
   * @brief Getter for bulk modulus.
   * @return reference to mutable bulk modulus.
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getBulkModulus() const override { return m_bulkModulus; }

  /**
   * @brief Accessor for shear modulus
   * @return A const reference to arrayView1d<real64> containing the shear
   *         modulus (at every element).
   */
  arrayView1d< real64 > const shearModulus() { return m_shearModulus; }

  /**
   * @brief Const accessor for shear modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         shear modulus (at every element).
   */
  arrayView1d< real64 const > const shearModulus() const { return m_shearModulus; }

  /**
   * @brief Getter for shear modulus.
   * @return reference to mutable shear modulus.
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getShearModulus() const override { return m_shearModulus; }

protected:
  virtual void postInputInitialization() override;

  // Tangent elastic bulk modulus parameters
  real64 m_b0;
  real64 m_b1;
  real64 m_b2;
  real64 m_b3;
  real64 m_b4;

  // Tangent elastic shear modulus parameters
  real64 m_g0;
  real64 m_g1;
  real64 m_g2;
  real64 m_g3;
  real64 m_g4;

  // Cruch curve parameters
  real64 m_p0;
  real64 m_p1;
  real64 m_p2;
  real64 m_p3;
  real64 m_p4;

  // Shear limit surface parameters
  real64 m_peakT1;
  real64 m_fSlope;
  real64 m_stren;
  real64 m_ySlope;

  // Nonassociativity parameter
  real64 m_beta;

  // Rate dependence paramters
  real64 m_t1RateDependence;
  real64 m_t2RateDependence;

  // Fracture energy release rate
  real64 m_fractureEnergyReleaseRate;

  // Cap shape parameter
  real64 m_cr;

  // Fluid parameters
  real64 m_fluidBulkModulus;
  real64 m_fluidInitialPressure;

  int m_creep;
  // deviatoric creep parameters
  real64 m_creepC0;
  real64 m_creepC1;
  // compaction creep parameters
  real64 m_creepA;
  real64 m_creepB;
  real64 m_creepC;

  // strain-hardening parameters
  real64 m_strainHardeningN;
  real64 m_strainHardeningK;

  /// The bulk modulus for each element/particle
  array1d< real64 > m_bulkModulus;

  /// The shear modulus for each element/particle
  array1d< real64 > m_shearModulus;
 
  ///State variable: The velocity gradient for each element/particle
  array3d< real64 > m_velocityGradient;

  ///State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  // State variable: The porosity for each element/particle
  array2d< real64 > m_porosity;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_GEOMECHANICS_HPP_ */
