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
 * @file DamageSpectralUtilities.hpp
 * @brief Helper functions to perform spectral decomposition of stresses.
 *
 * A detailed description of the calculations performed here can be found in
 * Jiang, Wen et al. "Three-dimensional phase-field modeling of porosity dependent intergranular fracture in UO2"
 * Computational Materials Science 171 (2020): 109269
 *
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_DAMAGESPECTRALUTILITIES_HPP_
#define GEOS_CONSTITUTIVE_SOLID_DAMAGESPECTRALUTILITIES_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>
#include "LvArray/src/output.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"

namespace geos
{

//get the positive part only of tensor T using spectral split
GEOS_HOST_DEVICE inline
void PositivePartOfTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& positivePart)[6] )
{
  real64 positiveEigs[6]={};
  for( int i=0; i < 3; i++ )
  {
    positiveEigs[i] = fmax( 0.0, eigs[i] );
  }
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( positivePart, eigvecs, positiveEigs );
}

//get the negative part only of tensor T using spectral split
GEOS_HOST_DEVICE inline
void NegativePartOfTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& negativePart)[6] )
{
  real64 negativeEigs[6]={};
  for( int i=0; i < 3; i++ )
  {
    negativeEigs[i] = fmin( 0.0, eigs[i] );
  }
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( negativePart, eigvecs, negativeEigs );
}

//implements the : operator between two second-order tensors in voigt form
GEOS_HOST_DEVICE inline
real64 doubleContraction( real64 (& A)[6], real64 (& B)[6] )
{
  real64 ans = 0;
  for( int i=0; i < 6; i++ )
  {
    if( i < 3 )
    {
      ans = ans + A[i]*B[i];
    }
    else
    {
      ans = ans + 2*A[i]*B[i];
    }
  }
  return ans;
}


//heaviside function, return 1 for positive, 0 for negatives and 0.5 if argument is zero or very close
GEOS_HOST_DEVICE inline
real64 heaviside( real64 x )
{
  if( fabs( x ) < 1e-12 )
  {
    return 0.5;
  }
  if( x > 0 )
  {
    return 1;
  }
  if( x <= 0 )
  {
    return 0;
  }
  GEOS_ERROR( "This function was not supposed to reach this line" );
  return 1000000;
}

//computes a tensor that enters the calculation of the Jacobian of the Spectral Split - check reference paper for more details
GEOS_HOST_DEVICE inline
void QTensor( real64 const (&eigvector)[3], real64 (& Q)[6][6] )
{
  real64 M[6]={0};
  LvArray::tensorOps::symRij_eq_AiAj< 3 >( M, eigvector );
  for( int i = 0; i<6; i++ )
  {
    for( int j = 0; j<6; j++ )
    {
      Q[i][j] = M[ i ] * M[ j ];
    }
  }
}

//computes another tensor that enters the calculation of the Jacobian of the Spectral Split - check reference paper for more details
GEOS_HOST_DEVICE inline
void GTensor( real64 (& eigvec1)[3], real64 (& eigvec2)[3], real64 (& G)[6][6] )
{
  GEOS_UNUSED_VAR( eigvec1, eigvec2, G );
  real64 M1[6]={0};
  real64 M2[6]={0};
  LvArray::tensorOps::symRij_eq_AiAj< 3 >( M1, eigvec1 );
  LvArray::tensorOps::symRij_eq_AiAj< 3 >( M2, eigvec2 );

  G[0][0] = M1[0]*M2[0] + M1[0]*M2[0];
  G[0][1] = M1[5]*M2[5] + M1[5]*M2[5];
  G[0][2] = M1[4]*M2[4] + M1[4]*M2[4];
  G[0][3] = M1[5]*M2[4] + M1[4]*M2[5];
  G[0][4] = M1[0]*M2[4] + M1[4]*M2[0];
  G[0][5] = M1[0]*M2[5] + M1[5]*M2[0];
  G[1][0] = M1[5]*M2[5] + M1[5]*M2[5];
  G[1][1] = M1[1]*M2[1] + M1[1]*M2[1];
  G[1][2] = M1[3]*M2[3] + M1[3]*M2[3];
  G[1][3] = M1[1]*M2[3] + M1[3]*M2[1];
  G[1][4] = M1[5]*M2[3] + M1[3]*M2[5];
  G[1][5] = M1[5]*M2[1] + M1[1]*M2[5];
  G[2][0] = M1[4]*M2[4] + M1[4]*M2[4];
  G[2][1] = M1[3]*M2[3] + M1[3]*M2[3];
  G[2][2] = M1[2]*M2[2] + M1[2]*M2[2];
  G[2][3] = M1[3]*M2[2] + M1[2]*M2[3];
  G[2][4] = M1[4]*M2[2] + M1[2]*M2[4];
  G[2][5] = M1[4]*M2[3] + M1[3]*M2[4];
  G[3][0] = M1[5]*M2[4] + M1[5]*M2[4];
  G[3][1] = M1[1]*M2[3] + M1[1]*M2[3];
  G[3][2] = M1[3]*M2[2] + M1[3]*M2[2];
  G[3][3] = M1[1]*M2[2] + M1[3]*M2[3];
  G[3][4] = M1[5]*M2[2] + M1[3]*M2[4];
  G[3][5] = M1[5]*M2[3] + M1[1]*M2[4];
  G[4][0] = M1[0]*M2[4] + M1[0]*M2[4];
  G[4][1] = M1[5]*M2[3] + M1[5]*M2[3];
  G[4][2] = M1[4]*M2[2] + M1[4]*M2[2];
  G[4][3] = M1[5]*M2[2] + M1[4]*M2[3];
  G[4][4] = M1[0]*M2[2] + M1[4]*M2[4];
  G[4][5] = M1[0]*M2[3] + M1[5]*M2[4];
  G[5][0] = M1[0]*M2[5] + M1[0]*M2[5];
  G[5][1] = M1[5]*M2[1] + M1[5]*M2[1];
  G[5][2] = M1[4]*M2[3] + M1[4]*M2[3];
  G[5][3] = M1[5]*M2[3] + M1[4]*M2[1];
  G[5][4] = M1[0]*M2[3] + M1[4]*M2[5];
  G[5][5] = M1[0]*M2[1] + M1[5]*M2[5];
}

//this function takes the eigenvectors and eigenvalues of a tensor and builds the associated positive projector
/**
 *@brief This function takes the eigen-decomposition of a tensor and builds the 4th Positive Projector associated with it
 *@param[in] eigs array with the 3 eigenvalues of a tensor
 *@param[in] eigvecs 3x3 array with the 3 eigenvectors of a tensor (in rows)
 *@param[out] PositiveProjector empty array that will be populated with the voigt form of the Positive Projector
 *
 * Given a symmetric tensor T, the positive projector is defined as P+ = variation(T+)/variation(T). That is, if we
 * define the function f+ to be the positive spectral part of T, then, P+ is just the variational derivative of f+.
 * Note that we don't take the tensor T as a parameter, only its eigenvectors and eigenvalues.
 *
 */
GEOS_HOST_DEVICE inline
void PositiveProjectorTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& PositiveProjector)[6][6] )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  real64 tol = 1e-12;
  if( fabs( eigs[0] - eigs[1] ) < tol || fabs( eigs[0]-eigs[2] ) < tol || fabs( eigs[1]-eigs[2] ) < tol )
  {
    repeatedEigenvalues = true;
  }

  //init QVoigt
  real64 Qi[6][6] = {};
  //init GVoigt
  real64 Gsym[6][6] = {};
  real64 Gji[6][6] = {};

  //compute projector
  for( int i = 0; i < 3; i++ )
  {
    real64 ithEigenVector[3] = {};
    ithEigenVector[0]=eigvecs[0][i];
    ithEigenVector[1]=eigvecs[1][i];
    ithEigenVector[2]=eigvecs[2][i];
    //First Part
    //compute Qi
    QTensor( ithEigenVector, Qi );
    //do update
    LvArray::tensorOps::scale< 6, 6 >( Qi, heaviside( eigs[i] ));
    LvArray::tensorOps::add< 6, 6 >( PositiveProjector, Qi );
    if( !repeatedEigenvalues )
    {

      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3]={};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        //Do update
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (fmax( eigs[i], 0.0 ) - fmax( eigs[j], 0.0 ))/(2*(eigs[i]-eigs[j])));
        LvArray::tensorOps::add< 6, 6 >( PositiveProjector, Gsym );
      }
    }
    else
    {
      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3]={};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (heaviside( eigs[i] ) + heaviside( eigs[j] ))/4 );
        //do update
        LvArray::tensorOps::add< 6, 6 >( PositiveProjector, Gsym );
      }
    }
  }

}

//This is the negative projector, check the documentation of the positive projector for more details.
GEOS_HOST_DEVICE inline
void NegativeProjectorTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& NegativeProjector)[6][6] )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  real64 tol = 1e-12;
  if( fabs( eigs[0] - eigs[1] ) < tol || fabs( eigs[0]-eigs[2] ) < tol || fabs( eigs[1]-eigs[2] ) < tol )
  {
    repeatedEigenvalues = true;
  }

  //init QVoigt
  real64 Qi[6][6] = {};
  //init GVoigt
  real64 Gsym[6][6] = {};
  real64 Gji[6][6] = {};

  //compute projector
  for( int i = 0; i < 3; i++ )
  {
    real64 ithEigenVector[3] = {};
    ithEigenVector[0]=eigvecs[0][i];
    ithEigenVector[1]=eigvecs[1][i];
    ithEigenVector[2]=eigvecs[2][i];
    //First Part
    //compute Qi
    QTensor( ithEigenVector, Qi );
    //do update
    LvArray::tensorOps::scale< 6, 6 >( Qi, heaviside( -eigs[i] ));
    LvArray::tensorOps::add< 6, 6 >( NegativeProjector, Qi );
    if( !repeatedEigenvalues )
    {

      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3] = {};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        //Do update
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (fmin( eigs[i], 0.0 ) - fmin( eigs[j], 0.0 ))/(2*(eigs[i]-eigs[j])));
        LvArray::tensorOps::add< 6, 6 >( NegativeProjector, Gsym );
      }
    }
    else
    {
      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3] = {};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (heaviside( -eigs[i] ) + heaviside( -eigs[j] ))/4 );
        //do update
        LvArray::tensorOps::add< 6, 6 >( NegativeProjector, Gsym );
      }
    }
  }

}

//this function tests the getStiffness function from DamageSpectral.hpp
GEOS_HOST_DEVICE inline
void getStiffnessTest( real64 (& c)[6][6], real64 (& strain)[6], real64 damage )
{

  //Spectral Split
  real64 const damageFactor = (1-damage)*(1-damage);
  real64 const mu = 1;
  real64 const lambda = 1;
  //get strain tensor in voigt form
  real64 traceOfStrain = strain[0] + strain[1] + strain[2];
  //get eigenvalues and eigenvectors
  real64 eigenValues[3]={};
  real64 eigenVectors[3][3]={};
  LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );
  //construct 4th order IxI tensor
  real64 IxITensor[6][6]={};
  for( int i=0; i < 3; i++ )
  {
    for( int j=0; j < 3; j++ )
    {
      IxITensor[i][j] = 1.0;
    }
  }

  //construct positive part
  real64 cPositive[6][6]={};
  real64 positiveProjector[6][6]={};
  PositiveProjectorTensor( eigenValues, eigenVectors, positiveProjector );
  LvArray::tensorOps::scaledCopy< 6, 6 >( cPositive, IxITensor, lambda*heaviside( traceOfStrain ));
  LvArray::tensorOps::scale< 6, 6 >( positiveProjector, 2*mu );
  LvArray::tensorOps::add< 6, 6 >( cPositive, positiveProjector );

  //construct negative part
  real64 negativeProjector[6][6]={};
  NegativeProjectorTensor( eigenValues, eigenVectors, negativeProjector );
  LvArray::tensorOps::scaledCopy< 6, 6 >( c, IxITensor, lambda*heaviside( -traceOfStrain ));
  LvArray::tensorOps::scale< 6, 6 >( negativeProjector, 2*mu );
  LvArray::tensorOps::add< 6, 6 >( c, negativeProjector );
  //finish up
  LvArray::tensorOps::scale< 6, 6 >( cPositive, damageFactor );
  LvArray::tensorOps::add< 6, 6 >( c, cPositive );

}

//this function tests the GetStress function of DamageSpectral.hpp
GEOS_HOST_DEVICE inline
void getTestStress( real64 (& strain)[6], real64 (& stress)[6] )
{

  //Spectral split
  real64 const damageFactor = 0.25;
  real64 const mu = 1;
  real64 const lambda = 1;
  //get strain tensor in voigt form
  real64 traceOfStrain = strain[0] + strain[1] + strain[2];
  //get eigenvalues and eigenvectors
  real64 eigenValues[3] = {};
  real64 eigenVectors[3][3] = {};
  LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );
  //transpose eigenVectors matrix to match convention
  real64 temp[3][3] = {};
  LvArray::tensorOps::transpose< 3, 3 >( temp, eigenVectors );
  LvArray::tensorOps::copy< 3, 3 >( eigenVectors, temp );
  real64 tracePlus = fmax( traceOfStrain, 0.0 );
  real64 traceMinus = fmin( traceOfStrain, 0.0 );
  //build symmetric matrices of positive and negative eigenvalues
  real64 Itensor[6] = {};
  for( int i = 0; i < 3; i++ )
  {
    Itensor[i] = 1;
  }
  real64 positivePartOfStrain[6] = {};
  real64 negativePartOfStrain[6] = {};
  PositivePartOfTensor( eigenValues, eigenVectors, positivePartOfStrain );
  NegativePartOfTensor( eigenValues, eigenVectors, negativePartOfStrain );
  real64 positiveStress[6] = {};
  real64 negativeStress[6] = {};
  LvArray::tensorOps::scaledCopy< 6 >( positiveStress, Itensor, lambda*tracePlus );
  LvArray::tensorOps::scaledCopy< 6 >( negativeStress, Itensor, lambda*traceMinus );
  LvArray::tensorOps::scaledAdd< 6 >( positiveStress, positivePartOfStrain, 2*mu );
  LvArray::tensorOps::scaledAdd< 6 >( negativeStress, negativePartOfStrain, 2*mu );
  LvArray::tensorOps::copy< 6 >( stress, negativeStress );
  LvArray::tensorOps::scaledAdd< 6 >( stress, positiveStress, damageFactor );
}

} //namespacegeos

#endif /* GEOS_CONSTITUTIVE_SOLID_DAMAGESPECTRALUTILITIES_HPP_ */
