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
 * @file BB_Tetrahedron.hpp
 */

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_BBTETRAHEDRON_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_BBTETRAHEDRON_HPP_

#include "FiniteElementBase.hpp"
#include <utility>



namespace geos
{
namespace finiteElement
{

/**
 * This class contains the kernel accessible functions specific to the
 * Bernstein-Bézier (BB) modal any-order tetrahedron finite element with 
 * Gaussian quadrature rules. Available functions are tailored for
 * Discontinuous Galerkin (DG) applications.
 * In barycentric coordinates l1, l2, l3, l4, the function indexed by 
 * (i1, i2, i3, i4) corresponds to the function
 * (i1+i2+i3+i4+3)! / (i1! i2! i3! i4!) l1^i1 l2^i2 l3^i3 l4^i4 
 * and they integrate to one over the reference element defined
 * by 0<=l1,l2,l3,l3<=1, l1+l2+l3+l4=1
 */
template< int ORDER >
class BB_Tetrahedron final : public FiniteElementBase
{
public:

  /// The number of shape functions per element.
  constexpr static localIndex numNodes = ( ORDER + 1 ) * ( ORDER + 2 ) * ( ORDER + 3 ) / 6;

  /// The number of shape functions per face
  constexpr static localIndex numNodesPerFace = ( ORDER + 1 ) * ( ORDER + 2 ) / 2;

  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = numNodes;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~BB_Tetrahedron() override
  {}

  GEOS_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  /**
   * @brief Get the number of quadrature points.
   * @param stack Stack variables as filled by @ref setupStack.
   * @return The number of quadrature points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static localIndex getNumQuadraturePoints( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numQuadraturePoints;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual localIndex getNumSupportPoints() const override
  {
    return numNodes;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual localIndex getMaxSupportPoints() const override
  {
    return maxSupportPoints;
  }

  /**
   * @brief Get the number of support points.
   * @param stack Object that holds stack variables.
   * @return The number of support points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static localIndex getNumSupportPoints( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numNodes;
  }


  /**
   * @brief Returns the determinant of the Jacobian of the element 
   * @param[in] X The coordinates of the tetrahedron
   * @return the (absolute value of the) determinant of the Jacobian on the element
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 jacobianDeterminant( real64 const (&X)[4][3] )
                         
  {
    real64 m[3][3] = {};
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    return LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) ) / 2.0;
  }

  /**
   * @brief Calculate the determinant of the jacobian on the face opposite to the given vertex 
   * @param[in] face The index of the vertex opposite to the desired face
   * @return the (absolute value of the) determinant of the Jacobian on the face 
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 faceJacobianDeterminant( localIndex face,
                                         real64 const (&X)[4][3] )
  {
    int i1 = ( face + 1 ) % 4;
    int i2 = ( face + 2 ) % 4;
    int i3 = ( face + 3 ) % 4;
    real64 ab[3] = { X[ i2 ][ 0 ] - X[ i1 ][ 0 ],
                     X[ i2 ][ 1 ] - X[ i1 ][ 1 ],
                     X[ i2 ][ 2 ] - X[ i1 ][ 2 ] };
    real64 ac[3] = { X[ i3 ][ 0 ] - X[ i1 ][ 0 ],
                     X[ i3 ][ 1 ] - X[ i1 ][ 1 ],
                     X[ i3 ][ 2 ] - X[ i1 ][ 2 ] };
    real64 term1 = ab[1] * ac[2] - ab[2] * ac[1];
    real64 term2 = ab[2] * ac[0] - ab[0] * ac[2];
    real64 term3 = ab[0] * ac[1] - ab[1] * ac[0];
    return LvArray::math::sqrt( ( term1 * term1 + term2 * term2 + term3 * term3 ) / 3.0 );
  }

  /**
   * @brief Calculate shape functions values at a single point using De Casteljau's algorithm.
   * @param[in] lambda barycentric coordinates of the point in thetetrahedron
   * @param[out] ORDER The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const ( & lambda)[4],
                     real64 (& N)[numNodes] )
  {
    N[ 0 ] = 6.0;
    int prev;
    int c;
    int limits[ 4 ] = { 1, 1, 1, 1 };
    for( int np = 1; np <= ORDER; np++)
    {
      prev = np * ( np + 1 ) * ( np + 2 ) / 6 - 1; 
      c = ( np + 1 ) * ( np + 2 ) * ( np + 3 ) / 6 - 1; 
      for( int i = 0; i < 4; i++ )
      {
        int denominator = i == 0 ? np : 1;
        int offset = 0;
        int c1 = np - 1;
        int c2 = i + np - 2;
        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
        for( int j = 0; j < limits[ i ] ; j++ )
        {
          if( j - offset >=  repetitionCount )
          {
            denominator++;
            offset += repetitionCount;
            repetitionCount = repetitionCount * c1 / c2;
            c1--;
            c2--;
          }
          N[ c-- ] = N[ prev - j ] * lambda[ 3 - i ] * ( np + 3 ) / denominator;
        }
      }
      for( int i = 1; i < 4; i++ )
      {
        limits[ i ] += limits[ i - 1 ];
      }
    } 
  }

  /**
      * @brief Calculate shape functions values at a single point, given the coordinates of the tetrahedron vertices, using De Casteljau's algorithm.
      * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
      * @param[out] ORDER The shape function values.
      */
     GEOS_HOST_DEVICE
     GEOS_FORCE_INLINE
     static void calcN( real64 const (&X)[4][3],
                        real64 const (&coords)[3],
                        real64 (& N)[numNodes] )
     {
       real64 lambda[4] = {};
       real64 m[3][3] = {};
       for( int i = 0; i < 3; i++ )
       {
         for( int j = 0; j < 3; j++ )
         {
           m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
         }
       }
       real64 den = LvARray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) );
       for( int i = 0; i < 3; i++ )
       {
         for( int j = 0; j < 3; j++ )
         {
           m[ i ][ j ] = coords[ j ] - X[ 0 ][ j ];
         }
         lambda[ i + 1 ] = LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) ) / den;
         for( int j = 0; j < 3; j++ )
         {
           m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
         }
       }
       lambda[ 0 ] = 1.0 - lambda[ 1 ] - lambda[ 2 ] - lambda[ 3 ];
       return calcN( lambda, N );
     }

     /**
      * @brief Calculate the values and derivatives of shape functions with respect to barycentric coordinates at a single point using De Casteljau's algorithm.
      * @param[in] lambda barycentric coordinates of the point in thetetrahedron
      * @param[out] N The shape function values.
      * @param[out] gradN The derivatives of the shape functions with respect to the lambdas
      */
     GEOS_HOST_DEVICE
     GEOS_FORCE_INLINE
     static void calcNandGradN( real64 const ( & lambda)[4],
                                real64 const ( & N)[numNodes],
                                real64 (& gradN)[numNodes][ 4 ] )
     {                  
       gradN[ 0 ][ 0 ] = 0.0;
       gradN[ 0 ][ 1 ] = 0.0;
       gradN[ 0 ][ 2 ] = 0.0;
       gradN[ 0 ][ 3 ] = 0.0;
       N[ 0 ] = 6.0;
       int prev;
       int c;
       int limits[ 4 ] = { 1, 1, 1, 1 };
       for( int np = 1; np <= ORDER; np++)
       {
         prev = np * ( np + 1 ) * ( np + 2 ) / 6 - 1; 
         c = ( np + 1 ) * ( np + 2 ) * ( np + 3 ) / 6 - 1; 
         for( int i = 0; i < 4; i++ )
         {
           int denominator = i == 0 ? np : 1;
           int offset = 0;
           int c1 = np - 1;
        int c2 = i + np - 2;
        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
        for( int j = 0; j < limits[ i ] ; j++ )
        {
          if( j - offset >=  repetitionCount )
          {
            denominator++;
            offset += repetitionCount;
            repetitionCount = repetitionCount * c1 / c2;
            c1--;
            c2--;
          }
          gradN[ c ][ 0 ] = gradN[ prev - j ][ 0 ] * lambda[ 3 - i ] * ( np + 3 ) / denominator;
          gradN[ c ][ 1 ] = gradN[ prev - j ][ 1 ] * lambda[ 3 - i ] * ( np + 3 ) / denominator;
          gradN[ c ][ 2 ] = gradN[ prev - j ][ 2 ] * lambda[ 3 - i ] * ( np + 3 ) / denominator;
          gradN[ c ][ 3 ] = gradN[ prev - j ][ 3 ] * lambda[ 3 - i ] * ( np + 3 ) / denominator;
          gradN[ c ][ 3 - i ] += N[ prev - j ] * ( np + 3 ) / denominator;
          N[ c-- ] = N[ prev - j ] * lambda[ 3 - i ] * ( np + 3 ) / denominator;
        }
      }
      for( int i = 1; i < 4; i++ )
      {
        limits[ i ] += limits[ i - 1 ];
      }
    } 
  }


  /**
   * @brief Calculate the shape functions values and derivatives at a single point, given the coorginates of the tetrahedron vertices, using De Casteljau's algorithm.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] ORDER The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcNandGradN( real64 const (&X)[4][3],
                             real64 const (&coords)[3],
                             real64 (& N)[numNodes],
                             real64 (& gradN)[numNodes][ 4 ] )
  {
    real64 lambda[4] = {};
    real64 m[3][3] = {};
    real64 dNdLambda[numNodes][4];
    real64 dLambdadX[4][4];
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    real64 den = LvARray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) );
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = coords[ j ] - X[ 0 ][ j ];
      }
      lambda[ i + 1 ] = LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) ) / den;
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    lambda[ 0 ] = 1.0 - lambda[ 1 ] - lambda[ 2 ] - lambda[ 3 ];
    calcNandGradN( lambda, N, dNdLambda );
    for( int i = 0; i < numNodes; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        gradN[ i ][ j ] = ( ( ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) )* ( dNdLambda[ i ][ 1 ] - dNdLambda[ i ][ 0 ] ) +
                            ( ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) )* ( dNdLambda[ i ][ 2 ] - dNdLambda[ i ][ 0 ] ) +
                            ( ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) )* ( dNdLambda[ i ][ 3 ] - dNdLambda[ i ][ 0 ] )) / den 
      }
    }
  }

//  /**
//   * @brief Calculate shape functions values at a single point on a face using De Casteljau's algorithm.
//   * @param[in] lambda barycentric coordinates of the point in triangle
//   * @param[out] N The shape function values.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void calcN( real64 const ( & lambda)[3],
//                     real64 (& N)[numNodesPerFace] )
//  {
//    N[ 0 ] = 2.0;
//    int prev;
//    int c;
//    int limits[ 3 ] = { 1, 1, 1 };
//    for( int np = 1; np <= ORDER; np++)
//    {
//      prev = np * ( np + 1 ) / 2 - 1; 
//      c = ( np + 1 ) * ( np + 2 ) / 2 - 1; 
//      for( int i = 0; i < 3; i++ )
//      {
//        int denominator = i == 0 ? np : 1;
//        int offset = 0;
//        int c1 = np - 1;
//        int c2 = i + np - 2;
//        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
//        for( int j = 0; j < limits[ i ] ; j++ )
//        {
//          if( j - offset >=  repetitionCount )
//          {
//            denominator++;
//            offset += repetitionCount;
//            repetitionCount = repetitionCount * c1 / c2;
//            c1--;
//            c2--;
//          }
//          N[ c-- ] = N[ prev - j ] * lambda[2 - i] * ( np + 2 ) / denominator;
//        }
//      }
//      for( int i = 1; i < 3; i++ )
//      {
//        limits[ i ] += limits[ i - 1 ];
//      }
//    } 
//  }
//
//  /**
//   * @brief Calculate the derivatives of shape functions with respect to barycentric coordinates at a single point on a face using De Casteljau's algorithm.
//   * @param[in] lambda barycentric coordinates of the point in the triangle
//   * @param[out] ORDER The shape function values.
//   * @param[out] gradN The derivatives of the shape functions with respect to the lambdas
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void calcNandGradN( real64 const ( & lambda)[ 3 ],
//                             real64 const ( & N)[ numNodes ],
//                             real64 (& gradN)[ numNodes ][ 3 ] )
//  {                  
//    gradN[ 0 ][ 0 ] = 0.0;
//    gradN[ 0 ][ 1 ] = 0.0;
//    gradN[ 0 ][ 2 ] = 0.0;
//    N[ 0 ] = 2.0;
//    int prev;
//    int c;
//    int limits[ 3 ] = { 1, 1, 1 };
//    for( int np = 1; np <= ORDER; np++)
//    {
//      prev = np * ( np + 1 ) / 2 - 1; 
//      c = ( np + 1 ) * ( np + 2 ) / 2 - 1; 
//      for( int i = 0; i < 3; i++ )
//      {
//        int denominator = i == 0 ? np : 1;
//        int offset = 0;
//        int c1 = np - 1;
//        int c2 = i + np - 2;
//        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
//        for( int j = 0; j < limits[ i ] ; j++ )
//        {
//          if( j - offset >=  repetitionCount )
//          {
//            denominator++;
//            offset += repetitionCount;
//            repetitionCount = repetitionCount * c1 / c2;
//            c1--;
//            c2--;
//          }
//          gradN[ c ][ 0 ] = gradN[ prev - j ][ 0 ] * lambda[ 2 - i ] * ( np + 2 )/ denominator;
//          gradN[ c ][ 1 ] = gradN[ prev - j ][ 1 ] * lambda[ 2 - i ] * ( np + 2 )/ denominator;
//          gradN[ c ][ 2 ] = gradN[ prev - j ][ 2 ] * lambda[ 2 - i ] * ( np + 2 )/ denominator;
//          gradN[ c ][ 2 - i ] += N[ prev - j ] * ( np + 2 ) / denominator;
//          N[ c-- ] = N[ prev - j ] * lambda[ 2 - i ] * ( np + 2 ) / denominator;
//        }
//      }
//      for( int i = 1; i < 3; i++ )
//      {
//        limits[ i ] += limits[ i - 1 ];
//      }
//    } 
//  }

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the mesh support points.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           StackVariables const & GEOS_UNUSED_PARAM( stack ),
                           real64 ( &gradN )[numNodes][3] )
  {
    return calcGradN( q, X, gradN );
  }

  /**
   * @brief Computes a! / ( b! * c! ) with b + c >= a >= b, c
   * @param[in] a
   * @param[in] b
   * @param[in] c
   * @return a!/(b!*c!)
   */
   constexpr static real64 integralTerm(const int a, const int b, const int c)
   {
     real64 res = 1.0;
     int num = a;
     int den = c;
     for( int i = b; i > 0; i--)
     {
         res *= ( (real64) num ) /  i;
         num--;
     }
     for( int i = num; i > 0; i--)
     {
         res *= ( (real64) i ) /  den;
         den--;
     }
     for( int i = den; i > 0; i--)
     {
         res /= i;
     }
     
     return res;
  }

  /**
   * @brief Computes the superposition integral between Bernstein-Bézier functions indexed by 
   *  (i1, j1, k1, l1) and (i1, j2, k2, l2)
   * @param[in] i1
   * @param[in] j1
   * @param[in] k1
   * @param[in] l1
   * @param[in] i2
   * @param[in] j2
   * @param[in] k2
   * @param[in] l2
   * @return the superposition integral over the barycentric coordinates
   */
   constexpr static real64 computeSuperpositionIntegral( const int i1, const int j1, const int k1, const int l1, 
                                                         const int i2, const int j2, const int k2, const int l2 )
   {
     return (integralTerm(i1+i2, i1, i2)*
             integralTerm(j1+j2, j1, j2)*
             integralTerm(k1+k2, k1, k2)*
             integralTerm(l1+l2, l1, l2))/
             integralTerm(i1+j1+k1+l1+i2+j2+k2+l2+3, 
                          i1+j1+k1+l1+3, i2+j2+k2+l2+3);
   }

  /**
   * @brief Computes the superposition integral over a face between Bernstein-Bézier functions whose indices are given by 
   *  (i1, j1, k1, l1=0 ) and (i1, j2, k2, l2=0)
   * @param[in] i1
   * @param[in] j1
   * @param[in] k1
   * @param[in] i2
   * @param[in] j2
   * @param[in] k2
   * @return the superposition integral over the barycentric coordinates
   */
   constexpr static real64 computeFaceSuperpositionIntegral( const int i1, const int j1, const int k1, 
                                                             const int i2, const int j2, const int k2 )
   {
     return ((i1+k1+j1+3)*(i2+j2+k2+3)*
             integralTerm(i1+i2, i1, i2)*
             integralTerm(j1+j2, j1, j2)*
             integralTerm(k1+k2, k1, k2))/
             integralTerm(i1+j1+k1+i2+j2+k2+2, 
                          i1+j1+k1+2, i2+j2+k2+2);
   }

  /**
   * @brief Computes the local degree of freedom index given the shape function indices (i, j, k, l) for each vertex.
   *   Only i, j and k are needed since i + j + k + l = order
   * @tparam I The index with respect to the first vertex
   * @tparam J The index with respect to the second vertex
   * @tparam K The index with respect to the third vertex
   * @return The local degree of freedom index
   */
  template< int I, int J, int K >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static
  constexpr
  int
  dofIndex()
  {
    return ( ORDER - I ) * ( ORDER - I + 1 ) * ( ORDER - I + 2) / 6 +
           ( ORDER - I - J ) * ( ORDER - I - J + 1 ) / 2 +
           ( ORDER - I - J - K ); 
  }

  /**
   * @brief Computes the local degree of freedom index given the shape function indices for each vertex
   * @tparam C The dof index in the element
   * @tparam VTX the vertex with respect to 
   * @return i, j, k, l if VTX= 0, 1, 2 or 3 resepctively (with i + j + k + l = order)
   */
  template< int C, int VTX  >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static
  constexpr
  int
  indexToIJKL()
  {
    static_assert( VTX >= 0 && VTX < 4 ); 
    // compute the indices of c in the current element using tetrahedral and triangular roots
    constexpr int cc1 = C + 1;
    constexpr real64 tetr = cbrt( 3.0 * cc1 + sqrt( 9.0 * cc1 * cc1 - 1.0 / 27.0 ) )
                            + cbrt( 3.0 * cc1 - sqrt( 9.0 * cc1 * cc1 - 1.0 / 27.0 ) ) - 2;
    constexpr int i = ORDER - round( tetr * 10.0 ) / 10;
    constexpr int cc2 = C - ( ORDER - i ) * ( ORDER - i + 1 ) * ( ORDER - i + 2 ) / 6 + 1;
    constexpr real64 trir = ( sqrt( 8.0 * cc2 + 1.0 ) - 1.0 ) / 2.0 - 1;
    constexpr int j = ORDER - i - round( trir * 10.0 ) / 10;
    constexpr int k = ORDER - i - j - ( C - (ORDER - i ) * ( ORDER - i + 1 ) * ( ORDER - i + 2 ) / 6
                      - ( ORDER - i - j ) * ( ORDER - i - j + 1 ) / 2 );
    if constexpr ( VTX == 0 )
    {
      return i;
    }
    else if constexpr ( VTX == 1)
    {
      return j;
    } 
    else if constexpr ( VTX == 2)
    {
      return k;
    } 
    else if constexpr ( VTX == 3)
    {
      return ORDER - i - j - k;
    } 
    return -1;
  }


  /**
   * @brief Helper function for static for loop
   * @tparam FUNC the callback function
   * @tparam ...Is integer indices of the loop
   */
  template < typename FUNC, int... Is >
  static constexpr void loop( FUNC && func, std::integer_sequence< int, Is... > )
  {
      ( func( std::integral_constant< int, Is >{} ), ... );
  }
  
  /**
   * @brief Helper function for loop over barycentric coordinates
   * @tparam FUNC the callback function
   */
  template < typename FUNC >
  static constexpr void barycentricCoordinateLoop(FUNC && func) {
      loop( [&] ( auto const i ) {
        func( std::integral_constant< int, i >{} );
      }, std::make_integer_sequence< int, 4 >{} );
  }

  /**
   * @brief Helper function for loop over tet basis functions
   * @tparam FUNC the callback function
   */
  template < typename FUNC >
  static constexpr void basisLoop(FUNC && func) {
    loop( [&] ( auto const i )
    {
      constexpr int i1 = ORDER - i;
      loop( [&] ( auto const j )

        if constexpr ( j1 <= ORDER - i1 )
        {
          loop( [&] ( auto const k )
          {
            constexpr int k1 = ORDER - k;
            if constexpr ( k1 <= ORDER - i1 - j1 )
            {
              constexpr int l1 = ORDER - i1 - j1 - k1;
              constexpr int c1 = dofIndex< i1, j1, k1 >();
              func( std::integral_constant< int, c1 >{}, 
                    std::integral_constant< int, i1 >{},
                    std::integral_constant< int, j1 >{},
                    std::integral_constant< int, k1 >{},
                    std::integral_constant< int, l1 >{} );
            }
        }, std::make_integer_sequence< int, ORDER + 1 > {} );
      }
      }, std::make_integer_sequence< int, ORDER + 1 > {} );
    }, std::make_integer_sequence< int, ORDER + 1 > {} );
  }

  /**
   * @brief Helper function for loop over tet basis functions that have one index in a given set of indices.
   *   If multiple indices are in the given list, the callback is called multiple times. 
   * @tparam FUNC the callback function
   * @tparam Is the setindices
   */
  template < typename FUNC, int Is.. >
  static constexpr void conditionalBasisLoop(FUNC && func) {
    loop( [&] ( auto const i )
    {
      constexpr int i1 = ORDER - i;
      loop( [&] ( auto const j )

        if constexpr ( j1 <= ORDER - i1 )
        {
          loop( [&] ( auto const k )
          {
            constexpr int k1 = ORDER - k;
            if constexpr ( k1 <= ORDER - i1 - j1 )
            {
              constexpr int l1 = ORDER - i1 - j1 - k1;
              constexpr int c1 = dofIndex< i1, j1, k1 >();
              ( ( (i1 == Is) &&
               (void( func( 0,
                      std::integral_constant<int, i1>{},
                      std::integral_constant<int, c1>{}, 
                      std::integral_constant<int, j1>{},
                      std::integral_constant<int, k1>{},
                      std::integral_constant<int, l1>{}) ), 1 ) ) || ... );
              ( ( (j1 == Is) &&
               (void( func( 1,
                     std::integral_constant<int, j1>{},
                     std::integral_constant<int, c1>{}, 
                     std::integral_constant<int, k1>{},
                     std::integral_constant<int, l1>{},
                     std::integral_constant<int, i1>{}) ), 1 ) ) || ... );
              ( ( (k1 == Is) &&
               ( void( func( 2,
                     std::integral_constant<int, k1>{},
                     std::integral_constant<int, c1>{},
                     std::integral_constant<int, l1>{},
                     std::integral_constant<int, i1>{},
                     std::integral_constant<int, j1>{}) ), 1 ) ) || ... );
              ( ( (l1 == Is) &&
               ( void( func( 3,
                     std::integral_constant<int, l1>{},
                     std::integral_constant<int, c1>{}, 
                     std::integral_constant<int, i1>{},
                     std::integral_constant<int, j1>{},                 
                     std::integral_constant<int, k1>{}) ), 1 ) ) || ...);
            }
        }, std::make_integer_sequence< int, ORDER + 1 > {} );
      }
      }, std::make_integer_sequence< int, ORDER + 1 > {} );
    }, std::make_integer_sequence< int, ORDER + 1 > {} );
  }

  /**
   * @brief Computes the reference mass matrix, i.e., the superposition matrix of the shape functions
   *   in barycentric coordinates. The real-world mass matrix can be obtained by using the multiplying
   *   this matrix by the determinant of the Jacobian. 
   *   
   * @param[out] m The mass matrix 
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static
  void
  computeMassMatrix( real64 (& m)[numNodes][numNodes] )
  {
    basisLoop( [&] ( auto const c1, auto const i1; auto const j1, auto const k1, auto const l1)
    {
      basisLoop( [&] ( auto const c2, auto const i2; auto const j2, auto const k2, auto const l2)
      {
        constexpr real64 val = computeSuperpositionIntegral( i1, j1, k1, l1, i2, j2, k2, l2 ); 
        m[ c1 ][ c2 ] = val;
      } );
    } );
  }

  /**
   * @brief Computes the reference damping matrix, i.e., the superposition matrix of the shape functions
   *   in barycentric coordinates over faces. The real-world mass matrix can be obtained by using the multiplying
   *   this matrix by the determinant of the face Jacobian. 
   *   
   * @param[out] d The damping matrix 
   * @param[in] face1Damped Whether the first face contributes to the damping term 
   * @param[in] face2Damped Whether the second face contributes to the damping term 
   * @param[in] face3Damped Whether the third face contributes to the damping term 
   * @param[in] face4Damped Whether the fourth face contributes to the damping term 
   */
  template< typename DAMPING >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static
  void
  computeReferenceDampingMatrix( real64 (& d)[numNodes][numNodes], bool const face1Damped, bool const face2Damped, bool const face3Damped, bool const face4Damped)
  {
            
    for( int c1 = 0; c1 < numNodes; c1++ )
    {
      for( int c2 = 0; c2 < numNodes; c2++ )
      {
        d[ c1 ][ c2 ] = 0;
      }
    }
    conditionalBasisLoop< 0 >( [&] ( auto const f1, auto const, auto const c1; auto const i1, auto const j1, auto const k1)
    {
      conditionalBasisLoop< 0 >( [&] ( auto const f2, auto const, auto const c2; auto const i2, auto const j2, auto const k2)
      {
        if constexpr( f1 == f2 )
        {
          constexpr real64 val = computeFaceSuperpositionIntegral( j1, k1, l1, j2, k2, l2 );
          if ( ( f1 == 0 && face1Damped ) ||
               ( f1 == 1 && face2Damped ) ||
               ( f1 == 2 && face3Damped ) ||
               ( f1 == 3 && face4Damped ) )
          { 
            d[ c1 ][ c2 ] += val;
          }
        }
      } );
    } );
  }

  /**
   * @brief Computes the non-zero contributions inside the element of the
   *   mass matrix M, i.e., the superposition matrix of shape functions
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j (local d.o.f. inside the element) and M_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static
  constexpr
  void
  computeMassTerm( real64 const (&X)[4][3],
                   FUNC && func )
  {
    real64 detJ = jacobianDeterminant( X );
    basisLoop( [&] ( auto const c1, auto const i1, auto const j1, auto const k1, auto const l1 )
    {
      basisLoop( [&] ( auto const c2, auto const i2, auto const j2, auto const k2, auto const l2 )
      {
        constexpr real64 val = computeSuperpositionIntegral( i1, j1, k1, l1, i2, j2, k2, l2 );
        func( c1, c2, val * detJ );                
      } );
    } );
  }

  /**
   * @brief Computes the non-zero contributions inside the element of the
   *   stiffness matrix R, i.e., the superposition matrix of first derivatives of shape functions
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j (local d.o.f. inside the element) and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static
  constexpr
  void
  computeStiffnessTerm( real64 const (&X)[4][3],
                        FUNC && func )
  {
    real64 detJ = jacobianDeterminant( X );
    real64 dLambdadX[4][3] = {};
    for( int j = 0; j < 3; j ++ )
    {
      dLambdadX[1][j] = ( ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) ) / detJ;
      dLambdadX[2][j] = ( ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) ) / detJ;
      dLambdadX[3][j] = ( ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) ) / detJ;
      dLambdadX[0][j] = -dLambdadX[1][j] - dLambdadX[2][j] - dLambdadX[3][j];
    }
    basisLoop( [&] ( auto const c1, auto const i1, auto const j1, auto const k1, auto const l1 )
    {
      basisLoop( [&] ( auto const c2, auto const i2, auto const j2, auto const k2, auto const l2 )
      {
        barycentricCoordinateLoop( [&] ( auto const d1 ) 
        {
          barycentricCoordinateLoop( [&] ( auto const d2 ) 
          {
            constexpr real64 factor1 = ( i1 + j1 + k1 + l1 + 3 ) / ( i1 +j1 + k1 + l1 + 2);
            constexpr real64 factor2 = ( i2 + j2 + k2 + l2 + 3 ) / ( i2 +j2 + k2 + l2 + 2);
            constexpr int ii1 = ( d1 != 0 ) * i1 + ( d1 == 0 ) * ( i1 - 1 );
            constexpr int ij1 = ( d1 != 1 ) * j1 + ( d1 == 1 ) * ( j1 - 1 );
            constexpr int ik1 = ( d1 != 2 ) * k1 + ( d1 == 2 ) * ( k1 - 1 );
            constexpr int il1 = ( d1 != 3 ) * l1 + ( d1 == 3 ) * ( l1 - 1 );
            constexpr int ii2 = ( d2 != 0 ) * i2 + ( d2 == 0 ) * ( i2 - 1 );
            constexpr int ij2 = ( d2 != 1 ) * j2 + ( d2 == 1 ) * ( j2 - 1 );
            constexpr int ik2 = ( d2 != 2 ) * k2 + ( d2 == 2 ) * ( k2 - 1 );
            constexpr int il2 = ( d2 != 3 ) * l2 + ( d2 == 3 ) * ( l2 - 1 );
            if constexpr (ii1 >= 0 && ij1 >= 0 && ik1 >= 0 && il1 >= 0 &&
                          ii2 >= 0 && ij2 >= 0 && ik2 >= 0 && il2 >= 0)
            {
              constexpr real64 val = computeSuperpositionIntegral( ii1, ij1, ik1, il1, ii2, ij2, ik2, il2 ) * factor1 * factor2;
              func( c1, c2, val * (dLambdadX[d1][0]*dLambdadX[d2][0] + dLambdadX[d1][1]*dLambdadX[d2][1] + dLambdadX[d1][2]*dLambdadX[d2][2] );
            }
          } );
        } );
      } );
    } );
  }

  /**
   * @brief Computes the non-zero contributions inside the element of the surface terms, including the value of 
   *   the superposition integral of basis functions (used for the penalty and damping terms) and
   *   the superposition integral of the derivative of a function with the value of the other, used for the flux terms. 
   * @param X Array containing the coordinates of the support points.
   * @param funcP Callback function for non-zero penalty-type terms, accepting seven parameters: 
   *   c1, c2 (local d.o.f. inside the element), f1 (index of the face, i.e., index of the opposite vertex for this element), 
   *   i2, j2 and k2 (local indices for the second shape function) and value 
   * @param funcF Callback function for non-zero flux-type terms, accepting four parameters: 
   *   c2, c2 (local d.o.f. inside the element), f1 (index of the face, i.e., index of the opposite vertex for this element), and value 
   *   i2, j2 and k2 (local indices for the second shape function) and value 
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static
  constexpr
  void
  computeSurfaceTerms( real64 const (&X)[4][3],
                       FUNC && funcP,
                       FUNC && funcF )
  {
    real64 detJf = faceJacobianDeterminant( face, X );
    conditionalBasisLoop< 0, 1 >( [&] ( auto const f1, auto const d, auto const c1, auto const i1, auto const j1, auto const k1 )
    {
      conditionalBasisLoop< 0 >( [&] ( auto const f2, auto const, auto const c2, auto const i2, auto const j2, auto const k2 )
      {
        if constexpr ( f1 == f2 )
        {
          if constexpr( d == 0 )
          {
            constexpr real64 val = computeFaceSuperpositionIntegral( j1, k1, l1, j2, k2, l2 );
            funcP( c1, c2, f1, i2, j2, k2, val * detJf );
          }
          else if constexpr ( d == 1 )
          {
            constexpr real64 factor = ( i1 + j1 + k1 + 4 )/( i1 + j1 + k1 + 3 );
            constexpr real64 val = computeFaceSuperpositionIntegral( i1, j1, k1, i2, kj, k2 ) * factor; 
            funcF( c1, c2, f1, i2, j2, k2, val * detJf );
          }
        }
      } );
    } );
  }


//  /**
//   * @brief Calculate the integration weights for a quadrature point.
//   * @param q Index of the quadrature point.
//   * @param X Array containing the coordinates of the support points.
//   * @return The product of the quadrature rule weight and the determinate of
//   *   the parent/physical transformation matrix.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static real64 transformedQuadratureWeight( localIndex const q,
//                                             real64 const (&X)[numNodes][3] )
//  {
//    // TODO
//  }
//
//  /**
//   * @brief Calculates the isoparametric "Jacobian" transformation
//   *   matrix/mapping from the parent space to the physical space.
//   * @param q The quadrature point index
//   * @param X Array containing the coordinates of the mesh support points.
//   * @param J Array to store the Jacobian transformation.
//   * @return The determinant of the Jacobian transformation matrix.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static real64 invJacobianTransformation( int const q,
//                                           real64 const (&X)[8][3],
//                                           real64 ( & J )[3][3] )
//  {
//    // TODO
//  }
//
//
//  /**
//   * @brief Calculate the symmetric gradient of a vector valued support field
//   *   at a quadrature point using the stored inverse of the Jacobian
//   *   transformation matrix.
//   * @param q The quadrature point index
//   * @param invJ The inverse of the Jacobian transformation matrix.
//   * @param var The vector valued support field to apply the gradient
//   *   operator on.
//   * @param grad The symmetric gradient in Voigt notation.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void symmetricGradient( int const q,
//                                 real64 const (&invJ)[3][3],
//                                 real64 const (&var)[numNodes][3],
//                                 real64 ( &grad )[6] )
//  {
//    // TODO
//    
//  }
//
//  /**
//   * @brief Calculate the gradient of a vector valued support field at a point
//   *   using the stored basis function gradients for all support points.
//   * @param q The quadrature point index
//   * @param invJ The inverse of the Jacobian transformation matrix.
//   * @param var The vector valued support field to apply the gradient
//   *   operator on.
//   * @param grad The gradient.
//   *
//   * More precisely, the operator is defined as:
//   * \f[
//   * grad_{ij}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai}\right ),
//   * \f]
//   *
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void gradient( int const q,
//                        real64 const (&invJ)[3][3],
//                        real64 const (&var)[numNodes][3],
//                        real64 ( &grad )[3][3] )
//  {
//    // TODO
//  }
//
//  /**
//   * @brief Inner product of all basis function gradients and a rank-2
//   *   symmetric tensor evaluated at a quadrature point.
//   * @param q The 3d quadrature point index
//   * @param invJ The inverse of the Jacobian transformation matrix.
//   * @param var The rank-2 symmetric tensor at @p q.
//   * @param R The vector resulting from the tensor contraction.
//   *
//   * More precisely, the operator is defined as:
//   * \f[
//   * R_i = \sum_a^{nSupport} \left( \frac{\partial N_a}{\partial X_j} var_{ij} \right),
//   * \f]
//   * where \f$\frac{\partial N_a}{\partial X_j}\f$ is the basis function gradient,
//   *   \f$var_{ij}\f$ is the rank-2 symmetric tensor.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void plusGradNajAij( int const q,
//                              real64 const (&invJ)[3][3],
//                              real64 const (&var)[6],
//                              real64 ( &R )[numNodes][3] )
//  {
//    // TODO
//  }
//
//
//  /**
//   * @brief Calculates the isoparametric "Jacobian" transformation
//   *   matrix/mapping from the parent space to the physical space at a single point.
//   * @param coords The parent coordinates at which to evaluate the shape function value
//   * @param X Array containing the coordinates of the support points.
//   * @param J Array to store the Jacobian transformation.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void jacobianTransformation( real64 const (&coords)[3],
//                                      real64 const (&X)[numNodes][3],
//                                      real64 ( &J )[3][3] )
//  {
//    // TOOD
//  }
//
//
//  /**
//   * @brief computes the real-world coordinates of the support nodes
//   * @param[in] Xmesh Array containing the coordinates of the corners of the mesh element
//   * @param[out] X Array containing the coordinates of the support points.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void computeLocalCoords( real64 const (&Xmesh)[8][3],
//                                  real64 const (&X)[numNodes][3] )
//  {
//    // TODO 
//  }
//
//
//  
//  template< typename FUNC >
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void computeStiffnessxyTerm( localIndex const q,
//                                      real64 const (&X)[8][3],
//                                      FUNC && func )
//  { 
//    // TODO
//  }
//
//  /**
//   * @brief computes the non-zero contributions of the d.o.f. indexed by q to the
//   *   partial-stiffness matrix R, i.e., the superposition matrix of first derivatives in z only
//   *   of the shape functions. Warning, the matrix B is obtained by computeBzMatrix instead of usual one.
//   * @param q The quadrature point index
//   * @param X Array containing the coordinates of the support points.
//   * @param func Callback function accepting three parameters: i, j and R_ij
//   */
//  template< typename FUNC >
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void computeStiffnesszTerm( localIndex const q,
//                                     real64 const (&X)[8][3],
//                                     FUNC && func )
//  { 
//    // TODO
//  }
//
//
//
//  /**
//   * @brief Apply a Jacobian transformation matrix from the parent space to the
//   *   physical space on the parent shape function derivatives, producing the
//   *   shape function derivatives in the physical space.
//   * @param q The quadrature point index
//   * @param invJ The Jacobian transformation from parent->physical space.
//   * @param gradN Array to contain the shape function derivatives for all
//   *   support points at the coordinates of the quadrature point @p q.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void applyTransformationToParentGradients( int const q,
//                                                    real64 const ( &invJ )[3][3],
//                                                    real64 ( &gradN )[numNodes][3] )
//  {
//    // TODO
//  }
//
//  /**
//   * @brief Apply a Jacobian transformation matrix from the parent space to the
//   *   physical space on the parent shape function derivatives, producing the
//   *   shape function derivatives in the physical space at a single point.
//   * @param coords The parent coordinates at which to apply the transformation
//   * @param invJ The Jacobian transformation from parent->physical space.
//   * @param gradN Array to contain the shape function derivatives for all
//   *   support points at the coordinates of the quadrature point @p q.
//   */
//  GEOS_HOST_DEVICE
//  GEOS_FORCE_INLINE
//  static void applyTransformationToParentGradients( real64 const (&coords)[3],
//                                                    real64 const ( &invJ )[3][3],
//                                                    real64 ( &gradN )[numNodes][3] )
//  {
//    // TODO
//  }
//
};

/**
 * Fixed-degree classes
 */
using BB1_Tetrahedron = BB_Tetrahedron< 1 >;
using BB2_Tetrahedron = BB_Tetrahedron< 2 >;
using BB3_Tetrahedron = BB_Tetrahedron< 3 >;
using BB4_Tetrahedron = BB_Tetrahedron< 4 >;
using BB5_Tetrahedron = BB_Tetrahedron< 5 >;
using BB6_Tetrahedron = BB_Tetrahedron< 6 >;
using BB7_Tetrahedron = BB_Tetrahedron< 7 >;
using BB8_Tetrahedron = BB_Tetrahedron< 8 >;
using BB9_Tetrahedron = BB_Tetrahedron< 9 >;
using BB10_Tetrahedron = BB_Tetrahedron< 10 >;

}
}

#endif // GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_BBTETRAHEDRON_HPP_
