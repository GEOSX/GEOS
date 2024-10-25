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
 * @file SolidUtilities.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDUTILITIES_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDUTILITIES_HPP_

namespace geos
{
namespace constitutive
{

struct SolidUtilities
{

  /**
   * @brief Perform a finite-difference check of the stiffness computation
   *
   * This method uses several stress evaluations and finite differencing to
   * approximate the 6x6 stiffness matrix, and then computes an error between
   * the coded stiffness method and the finite difference version.
   *
   * @note This method only works for models providing the smallStrainUpdate
   * method returning a 6x6 stiffness.
   *
   * @param solid the solid kernel wrapper
   * @param k the element number
   * @param q the quadrature index
   * @param timeIncrement the time increment
   * @param strainIncrement strain increment (on top of which a FD perturbation will be added)
   * @param print flag to decide if debug output is printed or not
   */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static bool
  checkSmallStrainStiffness( SOLID_TYPE const & solid,
                             localIndex k,
                             localIndex q,
                             real64 const & timeIncrement,
                             real64 const ( &strainIncrement )[6],
                             bool print = false )
  {
    real64 stiffness[6][6]{};     // coded stiffness
    real64 stiffnessFD[6][6]{};   // finite difference approximation
    real64 stress[6]{};           // original stress

    solid.smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );
    SolidUtilities::computeSmallStrainFiniteDifferenceStiffness( solid, k, q, timeIncrement, strainIncrement, stiffnessFD );

    // compute relative error between two versions

    real64 error = 0;
    real64 norm = 0;

    for( localIndex i=0; i<6; ++i )
    {
      for( localIndex j=0; j<6; ++j )
      {
        error += fabs( stiffnessFD[i][j]-stiffness[i][j] );
        norm += fabs( stiffnessFD[i][j] );
      }
    }
    error /= norm;

    // optional printing for debugging purposes

    if( print )
    {
      for( localIndex i=0; i<6; ++i )
      {
        for( localIndex j=0; j<6; ++j )
        {
          printf( "[%8.1e vs %8.1e] ", stiffnessFD[i][j], stiffness[i][j] );
        }
        printf( "\n" );
      }
    }

    return (error < 1e-3);
  }

  /**
   * @brief Perform a finite-difference stiffness computation
   *
   * This method uses stress evaluations and finite differencing to
   * approximate the 6x6 stiffness matrix.
   *
   * @note This method only works for models providing the smallStrainUpdate
   * method returning a 6x6 stiffness, as it will primarily be used to check
   * the hand coded tangent against a finite difference reference.
   * A similar method would need to be implemented to check compressed stiffness,
   * stress-only, or finite-strain interfaces.
   *
   * @param solid the solid kernel wrapper
   * @param k the element number
   * @param q the quadrature index
   * @param timeIncrement the time increment
   * @param strainIncrement strain increment (on top of which a FD perturbation will be added)
   * @param stiffnessFD finite different stiffness approximation
   */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  computeSmallStrainFiniteDifferenceStiffness( SOLID_TYPE const & solid,
                                               localIndex k,
                                               localIndex q,
                                               real64 const & timeIncrement,
                                               real64 const ( &strainIncrement )[6],
                                               real64 ( & stiffnessFD )[6][6] )
  {
    real64 stiffness[6][6]{};      // coded stiffness
    real64 stress[6]{};            // original stress
    real64 stressFD[6]{};          // perturbed stress
    real64 strainIncrementFD[6]{}; // perturbed strain
    real64 norm = 0;             // norm for scaling (note: method is fragile w.r.t. scaling)

    for( localIndex i=0; i<6; ++i )
    {
      strainIncrementFD[i] = strainIncrement[i];
      norm += fabs( strainIncrement[i] );
    }

    real64 eps = 1e-4*norm; // finite difference perturbation

    solid.smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );

    for( localIndex i=0; i<6; ++i )
    {
      strainIncrementFD[i] += eps;

      if( i>0 )
      {
        strainIncrementFD[i-1] -= eps;
      }

      solid.smallStrainUpdate( k, q, timeIncrement, strainIncrementFD, stressFD, stiffnessFD );

      for( localIndex j=0; j<6; ++j )
      {
        stiffnessFD[j][i] = (stressFD[j]-stress[j])/eps;
      }
    }
  }

  /**
   * @brief Hypo update (small strain, large rotation).
   *
   * This function uses a call to the small strain update, followed by
   * a rotation correction using the Hughes-Winget incrementally objective
   * algorithm.  One can imagine the material deforming in small strain, but in
   * a reference frame that rotates to track the body's local rotation.  This
   * provides a convenient way to extend small strain models to a finite rotation
   * regime.  From the assumption of small deformations, the Cauchy stress and
   * Kirchoff stress are approximately equal (det F ~ 1).
   *
   * We use a post-rotation of the new stress (as opposed to pre-rotation of the old stress)
   * as we don't have to unrotate the old stress in the event a rewind is required.
   * We do not post-rotate the stiffness tensor, which is an approximation, but
   * should be sufficient for small rotation increments.
   *
   * This method does not work for anisotropic properties or yield functions
   * (without some care) and a co-rotational formulation should be considered instead.
   *
   * @param solid the solid kernel wrapper
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] timeIncrement The time increment
   * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt)
   * @param[in] Rot The incremental rotation tensor
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New stiffness value
   */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  hypoUpdate( SOLID_TYPE const & solid,
              localIndex const k,
              localIndex const q,
              real64 const & timeIncrement,
              real64 const ( &Ddt )[6],
              real64 const ( &Rot )[3][3],
              real64 ( & stress )[6],
              real64 ( & stiffness )[6][6] )
  {
    solid.smallStrainUpdate( k, q, timeIncrement, Ddt, stress, stiffness );

    real64 temp[6]{};
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, Rot, stress );
    LvArray::tensorOps::copy< 6 >( stress, temp );
    solid.saveStress( k, q, stress );
  }

  /**
   * @brief Hypo update 2 (large total strain but small incremental strain, large rotation).
   *
   * Taken from http://www.sci.utah.edu/publications/Kam2011a/Kamojjala_ECTC2011.pdf
   *
   * The base class rotates the beginning-of-step stress and rate-of-deformation back
   * to the reference configuration using the beginning-of-step rotation (found via
   * polar decomposition of the beginning-of-step deformation gradient). It then calls
   * the small strain update to incrementally update the stress, followed by a rotation
   * of stress back to the end-of-step configuration (found using polar decomposition
   * of the end-of-step deformation gradient). This should be valid for most constitutive
   * models being explicitly integrated since the time steps are small enough that any given
   * step can be assumed to behave like a small-strain deformation with pre-stress.
   *
   * Note that if the derived class has tensorial state variables (beyond the
   * stress itself) care must be taken to rotate these as well.
   *
   * This method should work as-is for anisotropic properties and yield functions.
   *
   * @param solid the solid kernel wrapper
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] timeIncrement The time increment
   * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt) WITHOUT factors of 2 on the shear terms
   * @param[in] RotBeginning Beginning-of-step rotation tensor (obtained via polar decomposition)
   * @param[in] RotEnd End-of-step rotation tensor (obtained via polar decomposition)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New stiffness value
   */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  hypoUpdate2( SOLID_TYPE const & solid,
               localIndex const k,
               localIndex const q,
               real64 const timeIncrement,
               real64 ( & Ddt )[6],
               real64 const ( &RotBeginning )[3][3],
               real64 const ( &RotEnd )[3][3],
               real64 ( & stress )[6],
               real64 ( & stiffness )[6][6] )
  {
    // Prepare strain increment for rotation
    Ddt[3] *= 0.5;
    Ddt[4] *= 0.5;
    Ddt[5] *= 0.5;

    // Rotate m_oldStress and Ddt from beginning-of-step configuration to reference configuration.
    real64 temp[6] = { 0 };
    real64 RotBeginningTranpose[3][3] = { {0} };
    LvArray::tensorOps::transpose< 3, 3 >( RotBeginningTranpose, RotBeginning ); // We require the transpose since we're un-rotating
    LvArray::tensorOps::copy< 6 >( temp, solid.m_oldStress[ k ][ q ] );
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( solid.m_oldStress[ k ][ q ], RotBeginningTranpose, temp );
    LvArray::tensorOps::copy< 6 >( temp, Ddt );
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( Ddt, RotBeginningTranpose, temp );

    // Convert strain increment to Voigt notation by re-introducing factors of 2 on shear terms
    Ddt[3] *= 2;
    Ddt[4] *= 2;
    Ddt[5] *= 2;

    // Stress increment
    solid.smallStrainUpdate( k, q, timeIncrement, Ddt, stress, stiffness );

    // Rotate final stress to end-of-step (current) configuration
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, RotEnd, solid.m_newStress[ k ][ q ] );
    LvArray::tensorOps::copy< 6 >( stress, temp );
    solid.saveStress( k, q, stress );
  }

/**
 * @brief Hypo update, returning only stress
 *
 * @param solid the solid kernel wrapper
 * @param[in] k The element index.
 * @param[in] q The quadrature point index.
 * @param[in] timeIncrement The time increment
 * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt)
 * @param[in] Rot The incremental rotation tensor
 * @param[out] stress New stress value (Cauchy stress)
 */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  hypoUpdate_StressOnly( SOLID_TYPE const & solid,
                         localIndex const k,
                         localIndex const q,
                         real64 const & timeIncrement,
                         real64 const ( &Ddt )[6],
                         real64 const ( &Rot )[3][3],
                         real64 ( & stress )[6] )
  {
    solid.smallStrainUpdate_StressOnly( k, q, timeIncrement, Ddt, stress );

    real64 temp[6]{};
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, Rot, stress );
    LvArray::tensorOps::copy< 6 >( stress, temp );
    solid.saveStress( k, q, stress );
  }

/**
 * @brief Hypo update 2, returning only stress
 *
 * @param solid the solid kernel wrapper
 * @param[in] k The element index.
 * @param[in] q The quadrature point index.
 * @param[in] timeIncrement The time increment
 * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt) WITHOUT factors of 2 on the shear terms
 * @param[in] RotBeginning Beginning-of-step rotation tensor (obtained via polar decomposition)
 * @param[in] RotEnd End-of-step rotation tensor (obtained via polar decomposition)
 * @param[out] stress New stress value (Cauchy stress)
 * @param[out] stiffness New stiffness value
 */
  template< typename SOLID_TYPE >
  GEOS_HOST_DEVICE
  static void
  hypoUpdate2_StressOnly( SOLID_TYPE const & solid,
                          localIndex const k,
                          localIndex const q,
                          real64 const timeIncrement,
                          real64 ( & Ddt )[6],
                          real64 const ( &RotBeginning )[3][3],
                          real64 const ( &RotEnd )[3][3],
                          real64 ( & stress )[6] )
  {
    // Prepare strain increment for rotation
    Ddt[3] *= 0.5;
    Ddt[4] *= 0.5;
    Ddt[5] *= 0.5;

    // Rotate m_oldStress and Ddt from beginning-of-step configuration to reference configuration.
    real64 temp[6] = { 0 };
    real64 RotBeginningTranpose[3][3] = { {0} };
    LvArray::tensorOps::transpose< 3, 3 >( RotBeginningTranpose, RotBeginning ); // We require the transpose since we're un-rotating
    LvArray::tensorOps::copy< 6 >( temp, solid.m_oldStress[ k ][ q ] );
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( solid.m_oldStress[ k ][ q ], RotBeginningTranpose, temp );
    LvArray::tensorOps::copy< 6 >( temp, Ddt );
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( Ddt, RotBeginningTranpose, temp );

    // Convert strain increment to Voigt notation by re-introducing factors of 2 on shear terms
    Ddt[3] *= 2;
    Ddt[4] *= 2;
    Ddt[5] *= 2;

    // Stress increment
    solid.smallStrainUpdate_StressOnly( k, q, timeIncrement, Ddt, stress );

    // Rotate final stress to end-of-step (current) configuration
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, RotEnd, solid.m_newStress[ k ][ q ] );
    LvArray::tensorOps::copy< 6 >( stress, temp );
    solid.saveStress( k, q, stress );
  }


};

} // namespace constitutive

} // namespace geos

#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDUTILITIES_HPP_ */
