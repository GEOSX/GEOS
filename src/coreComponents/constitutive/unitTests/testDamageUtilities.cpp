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

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "../solid/DamageSpectralUtilities.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "LvArray/src/output.hpp"

#include <gtest/gtest.h>
using namespace geos;

TEST( DamageUtilities, PositivePartOfTensor )
{
  real64 eigs[3] = {1.7076, 3.3973, 6.8951};
  real64 eigvecs[3][3] = {{-0.8643, 0.0759, 0.4973}, {0.1706, -0.8857, 0.4317}, {0.4732, 0.4579, 0.7526}};
  real64 positivePart[6] = {0};
  PositivePartOfTensor( eigs, eigvecs, positivePart );
  std::cout << positivePart[0]<<" "<<positivePart[5]<<" "<<positivePart[4]<<std::endl;
  std::cout <<"  "<< positivePart[1]<<" "<<positivePart[3]<<std::endl;
  std::cout <<"    "<<positivePart[2]<<std::endl;
}

TEST( DamageUtilities, NegativePartOfTensor )
{
  real64 eigs[3] = {2.0, 1.0, -1.0};
  real64 eigvecs[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  real64 negativePart[6] = {0};
  NegativePartOfTensor( eigs, eigvecs, negativePart );
  std::cout << negativePart[0]<<" "<<negativePart[5]<<" "<<negativePart[4]<<std::endl;
  std::cout <<"  "<< negativePart[1]<<" "<<negativePart[3]<<std::endl;
  std::cout <<"    "<<negativePart[2]<<std::endl;
}

TEST( DamageUtilities, PositiveAndNegative )
{
  srand ( static_cast< unsigned >(time( 0 )));
  float r1 = -1 + static_cast< float >(rand())/static_cast< float >(RAND_MAX/2);
  float r2 = -1 + static_cast< float >(rand())/static_cast< float >(RAND_MAX/2);
  float r3 = -1 + static_cast< float >(rand())/static_cast< float >(RAND_MAX/2);
  real64 eigs[3] = {r1, r2, r3};
  real64 eigvecs[3][3] = {};
  for( int i=0; i<3; i++ )
  {
    r1 = -1 + static_cast< float >(rand())/static_cast< float >(RAND_MAX/2);
    r2 = -1 + static_cast< float >(rand())/static_cast< float >(RAND_MAX/2);
    r3 = -1 + static_cast< float >(rand())/static_cast< float >(RAND_MAX/2);
    eigvecs[i][0] = r1;
    eigvecs[i][1] = r2;
    eigvecs[i][2] = r3;
  }
  real64 negativePart[6] = {0};
  real64 positivePart[6] = {0};
  PositivePartOfTensor( eigs, eigvecs, positivePart );
  NegativePartOfTensor( eigs, eigvecs, negativePart );
  LvArray::tensorOps::add< 6 >( positivePart, negativePart );
  real64 eigenvaluesVoigt[6] = {0};
  for( int i=0; i < 3; i++ )
  {
    eigenvaluesVoigt[i] = eigs[i];
  }
  real64 originalMatrix[6] = {0};
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( originalMatrix, eigvecs, eigenvaluesVoigt );
  for( int i=0; i<6; i++ )
  {
    std::cout << "original: "<<originalMatrix[i]<<"/computed: "<<positivePart[i]<<std::endl;
  }
}

TEST( DamageUtilities, QTensor )
{
  using namespace LvArray;
  real64 vector[3] = {1, 2, 3};
  real64 Q[6][6] = {};
  QTensor( vector, Q );
  std::cout << "For vector [1,2,3], Q is: "<< std::endl;
  GEOS_LOG( Q );
}

TEST( DamageUtilities, GTensor )
{
  using namespace LvArray;
  real64 vector1[3] = {1, 2, 3};
  real64 vector2[3] = {-1, -2, -3};
  real64 G[6][6] = {};
  GTensor( vector1, vector2, G );
  std::cout << "For vector1 = [1,2,3] and vector2 = [-1,-2,-3], G is: "<< std::endl;
  GEOS_LOG( G );
}

TEST( DamageUtilities, PositiveProjectorTensor )
{
  using namespace LvArray;
  //distinct eigenvalues
  real64 eigs[3] = {2.0, 1.0, -1.0};
  real64 eigvecs[3][3] = {{1, 1, 0}, {-2, 1, 0}, {-1, 0, 1}};
  real64 PositiveProjector[6][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
  PositiveProjectorTensor( eigs, eigvecs, PositiveProjector );
  std::cout << "P+ is: "<< std::endl;
  GEOS_LOG( PositiveProjector );
  //repeated eigenvalues
  real64 eigs2[3] = {1.0, 1.0, -1.0};
  real64 eigvecs2[3][3] = {{1, 1, 0}, {-2, 1, 0}, {-1, 0, 1}};
  real64 PositiveProjector2[6][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
  PositiveProjectorTensor( eigs2, eigvecs2, PositiveProjector2 );
  std::cout << "P+ is: "<< std::endl;
  GEOS_LOG( PositiveProjector2 );
  //all eigenvalues are equal
  real64 eigs3[3] = {1.0, 1.0, 1.0};
  real64 eigvecs3[3][3] = {{1, 1, 0}, {-2, 1, 0}, {-1, 0, 1}};
  real64 PositiveProjector3[6][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
  PositiveProjectorTensor( eigs3, eigvecs3, PositiveProjector3 );
  std::cout << "P+ is: "<< std::endl;
  GEOS_LOG( PositiveProjector3 );
  //some eigenvalues are zero
  real64 eigs4[3] = {1.0, 1.0, 0.0};
  real64 eigvecs4[3][3] = {{1, 1, 0}, {-2, 1, 0}, {-1, 0, 1}};
  real64 PositiveProjector4[6][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
  PositiveProjectorTensor( eigs4, eigvecs4, PositiveProjector4 );
  std::cout << "P+ is: "<< std::endl;
  GEOS_LOG( PositiveProjector4 );
}

TEST( DamageUtilities, NegativeProjectorTensor )
{
  using namespace LvArray;
  real64 eigs[3] = {1.0, 1.0, -1.0};
  real64 eigvecs[3][3] = {{1, 1, 0}, {-2, 1, 0}, {-1, 0, 1}};
  real64 NegativeProjector[6][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
  NegativeProjectorTensor( eigs, eigvecs, NegativeProjector );
  std::cout << "P- is: "<< std::endl;
  GEOS_LOG( NegativeProjector );
}

TEST( DamageUtilities, GetStiffness )
{
  using namespace LvArray;
  real64 strain[6] = {3, 4, 5, 1, 2, 1};
  real64 d = 0.5;
  real64 c[6][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
  getStiffnessTest( c, strain, d );
  std::cout << "The Stiffness is: " << std::endl;
  GEOS_LOG( c );
}

TEST( DamageUtilities, GetStress )
{
  using namespace LvArray;
  real64 strain[6] = {3, 4, 5, 1, 2, 1};
  real64 stress[6] = {0, 0, 0, 0, 0, 0};
  getTestStress( strain, stress );
  std::cout << "The stress tensor is: " << std::endl;
  GEOS_LOG( stress );
}
