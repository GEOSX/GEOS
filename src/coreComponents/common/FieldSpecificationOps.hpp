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
 * @file FieldSpecificationOps.hpp
 */

#ifndef GEOS_COMMON_FIELDSPECIFICATIONOPS_HPP
#define GEOS_COMMON_FIELDSPECIFICATIONOPS_HPP

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

/**
 * @brief OpEqual Operator that sets a value
 */
struct OpEqual
{
  /**
   * @brief Pointwise set of a value
   * @tparam T type of the left-hand side
   * @tparam U type of the right-hand side
   * @param[in] lhs value to set
   * @param[in] rhs input value
   */
  template< typename T, typename U >
  GEOS_HOST_DEVICE static inline
  void apply( T & lhs, U const & rhs )
  {
    lhs = static_cast< T >( rhs );
  }
};

/**
 * @brief OpAdd Operator that adds a value
 */
struct OpAdd
{
  /**
   * @brief Pointwise update of a value
   * @tparam T type of the left-hand side
   * @tparam U type of the right-hand side
   * @param[in] lhs value to update
   * @param[in] rhs input value
   */
  template< typename T, typename U >
  GEOS_HOST_DEVICE static inline
  void apply( T & lhs, U const & rhs )
  {
    lhs += static_cast< T >( rhs );
  }
};

/**
 * @brief FieldSpecificationOp
 */
template< typename OP >
struct FieldSpecificationOp
{
  /// Alias for OP, the operator
  using OpType = OP;

  /**
   * @brief Pointwise application of a value to a field
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component not used.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index] (+)= value.
   */
  template< typename T >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView1d< T > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    GEOS_UNUSED_VAR( component );
    OP::template apply( field( index ), value );
  }

  /**
   * @brief Pointwise application of value to a field variable.
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index][component] (+)= value.
   */
  template< typename T >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView1d< T > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    OP::template apply( field( index )[component], value );
  }

  /**
   * @brief Read a field value
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component not used.
   * @param[out] value The value read from @p field.
   * @return type of the input field value.
   *
   * This function performs value (+)= field[index].
   */
  template< typename T >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  readFieldValue( arrayView1d< T const > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    GEOS_UNUSED_VAR( component );
    OP::template apply( value, field( index ) );
  }

  /**
   * @brief Read a field value
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component not used.
   * @param[out] value The value read from @p field.
   * @return type of the input field value.
   *
   * This function performs value (+)= field[index][component].
   */
  template< typename T >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  readFieldValue( arrayView1d< T const > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    OP::template apply( value, field( index )[component] );
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The index along second dimension of 2d array.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index][component] (+)= value.
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView2d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      OP::template apply( field( index, component ), value );
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        OP::template apply( field( index, a ), value );
      }
    }
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index][component] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView2d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        OP::template apply( field( index, a )[component], value );
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex c = 0; c < T::SIZE; ++c )
        {
          OP::template apply( field( index, a )[c], value );
        }
      }
    }
  }

  /**
   * @brief Read value from a 2d field.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component The index along second dimension of 2d array.
   * @param[out] value The value that is read from @p field.
   * @return type of the input field value.
   *
   * This function performs value (+)= field[index][component].
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  readFieldValue( arrayView2d< T const, USD > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    GEOS_ASSERT( component >= 0 );
    OP::template apply( value, field( index, component ) );
  }

  /**
   * @brief This function is not meaningful. It exists for generic purposes, but will result in an error if called.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array2d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component The index along second dimension of 2d array.
   * @param[out] value The value that is read from @p field.
   * @return type of the input field value.
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  readFieldValue( arrayView2d< T const, USD > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    GEOS_UNUSED_VAR( field );
    GEOS_UNUSED_VAR( index );
    GEOS_UNUSED_VAR( component );
    GEOS_UNUSED_VAR( value );
    GEOS_ERROR( "readFieldValue: unsupported operation" );
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array3d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array3d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The index along third dimension of 3d array.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView3d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        OP::template apply( field( index, a, component ), value );
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          OP::template apply( field( index, a, b ), value );
        }
      }
    }
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array3d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array3d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index][component] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView3d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          OP::template apply( field( index, a, b )[component], value );
        }
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          for( localIndex c = 0; c < T::size(); ++c )
          {
            OP::template apply( field( index, a, b )[c], value );
          }
        }
      }
    }
  }

  /**
   * @brief This function is not meaningful. It exists for generic purposes, but will result in an error if called.
   * @tparam T The type of the array3d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array3d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component The index along second dimension of 2d array.
   * @param[out] value The value that is read from @p field.
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline void
  readFieldValue( arrayView3d< T const, USD > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    GEOS_UNUSED_VAR( field );
    GEOS_UNUSED_VAR( index );
    GEOS_UNUSED_VAR( component );
    GEOS_UNUSED_VAR( value );
    GEOS_ERROR( "readFieldValue: unsupported operation" );
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array4d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array4d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The index along third dimension of 3d array.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< !traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView4d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          OP::template apply( field( index, a, b, component ), value );
        }
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          for( localIndex c = 0; c < field.size( 3 ); ++c )
          {
            OP::template apply( field( index, a, b, c ), value );
          }
        }
      }
    }
  }

  /**
   * @brief Pointwise application of a value to a field variable.
   * @tparam T The type of the array4d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array4d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value to apply to @p field.
   * @return type of the field value.
   *
   * This function performs field[index][component] (+)= value for all values of field[index].
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline typename std::enable_if< traits::is_tensorT< T >, void >::type
  SpecifyFieldValue( arrayView4d< T, USD > const & field,
                     localIndex const index,
                     integer const component,
                     real64 const value )
  {
    if( component >= 0 )
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          for( localIndex c = 0; c < field.size( 3 ); ++c )
          {
            OP::template apply( field( index, a, b, c )[component], value );
          }
        }
      }
    }
    else
    {
      for( localIndex a = 0; a < field.size( 1 ); ++a )
      {
        for( localIndex b = 0; b < field.size( 2 ); ++b )
        {
          for( localIndex c = 0; c < field.size( 3 ); ++c )
          {
            for( localIndex d = 0; d < T::size(); ++d )
            {
              OP::template apply( field( index, a, b, c )[d], value );
            }
          }
        }
      }
    }
  }

  /**
   * @brief This function is not meaningful. It exists for generic purposes, but will result in an error if called.
   * @tparam T The type of the array4d field variable specified in @p field.
   * @tparam USD the unit stride dimension of the array @p field.
   * @param[in] field The array4d field variable to read @p value from.
   * @param[in] index The index in field to read @p value from.
   * @param[in] component The index along second dimension of 2d array.
   * @param[out] value The value that is read from @p field.
   */
  template< typename T, int USD >
  GEOS_HOST_DEVICE
  static inline void
  readFieldValue( arrayView4d< T const, USD > const & field,
                  localIndex const index,
                  integer const component,
                  real64 & value )
  {
    GEOS_UNUSED_VAR( field );
    GEOS_UNUSED_VAR( index );
    GEOS_UNUSED_VAR( component );
    GEOS_UNUSED_VAR( value );
    GEOS_ERROR( "readFieldValue: unsupported operation" );
  }

};

/**
 * @struct FieldSpecificationEqual
 * this struct a collection of static functions which adhere to an assumed interface for overwriting
 * a value for a field.
 */
struct FieldSpecificationEqual : public FieldSpecificationOp< OpEqual >
{
  /// Alias for FieldSpecificationOp< OpEqual >
  using base_type = FieldSpecificationOp< OpEqual >;
  using base_type::SpecifyFieldValue;

  /**
   * @brief Function to apply a Dirichlet like boundary condition to a single dof in a system of
   *        equations.
   * @param[in] dof The degree of freedom that is to be set.
   * @param[in] dofRankOffset offset of dof indices on current rank
   * @param[in,out] matrix the local part of the system matrix
   * @param[out] rhs The rhs contribution resulting from the application of the BC.
   * @param[in] bcValue The target value of the Boundary Condition
   * @param[in] fieldValue The current value of the variable to be set.
   *
   * This function clears the matrix row for the specified \p dof, sets the diagonal to some
   * appropriately scaled value, and sets \p rhs to the negative product of the scaled value
   * of the diagonal and the difference between \p bcValue and \p fieldValue.
   *
   * @note This function assumes the user is doing a Newton-type nonlinear solve and will
   * negate the rhs vector upon assembly. Thus, it sets the value to negative of the desired
   * update for the field. For a linear problem, this may lead to unexpected results.
   */
  GEOS_HOST_DEVICE
  static inline void
  SpecifyFieldValue( globalIndex const dof,
                     globalIndex const dofRankOffset,
                     CRSMatrixView< real64, globalIndex const > const & matrix,
                     real64 & rhs,
                     real64 const bcValue,
                     real64 const fieldValue )
  {
    globalIndex const localRow = dof - dofRankOffset;
    if( localRow >= 0 && localRow < matrix.numRows() )
    {
      arraySlice1d< globalIndex const > const columns = matrix.getColumns( localRow );
      arraySlice1d< real64 > const entries = matrix.getEntries( localRow );
      localIndex const numEntries = matrix.numNonZeros( localRow );

      real64 diagonal = 0;
      real64 const minDiagonal = 1e-15;
      for( localIndex j = 0; j < numEntries; ++j )
      {
        if( columns[ j ] == dof )
        {
          // check that the entry is large enough to enforce the boundary condition
          if( entries[j] >= 0 && entries[j] < minDiagonal )
          {
            entries[ j ] = minDiagonal;
          }
          else if( entries[j] < 0 && entries[j] > -minDiagonal )
          {
            entries[ j ] = -minDiagonal;
          }
          diagonal = entries[j];
        }
        else
        { entries[ j ] = 0; }
      }

      rhs = -diagonal * (bcValue - fieldValue);
    }
    else
    {
      rhs = 0.0;
    }
  }

  /**
   * @brief Function to add some values of a vector.
   * @tparam POLICY the execution policy to use when setting values
   * @param rhs the target right-hand side vector
   * @param dof a list of global DOF indices to be set
   * @param dofRankOffset offset of dof indices on current rank
   * @param values a list of values corresponding to \p dof that will be added to \p rhs.
   */
  template< typename POLICY >
  static inline void prescribeRhsValues( arrayView1d< real64 > const & rhs,
                                         arrayView1d< globalIndex const > const & dof,
                                         globalIndex const dofRankOffset,
                                         arrayView1d< real64 const > const & values )
  {
    GEOS_ASSERT_EQ( dof.size(), values.size() );
    forAll< POLICY >( dof.size(), [rhs, dof, dofRankOffset, values] GEOS_HOST_DEVICE ( localIndex const a )
    {
      globalIndex const localRow = dof[ a ] - dofRankOffset;
      if( localRow >= 0 && localRow < rhs.size() )
      { rhs[ localRow ] = values[ a ]; }
    } );
  }
};

/**
 * @struct FieldSpecificationAdd
 * this struct a collection of static functions which adhere to an assumed interface for adding
 * a value for a field.
 */
struct FieldSpecificationAdd : public FieldSpecificationOp< OpAdd >
{
  /// Alias for FieldSpecificationOp< OpAdd >
  using base_type = FieldSpecificationOp< OpAdd >;
  using base_type::SpecifyFieldValue;

  /**
   * @brief Function to apply a value to a vector field for a single dof.
   * @param[in] dof The degree of freedom that is to be modified.
   * @param[in] dofRankOffset offset of dof indices on current rank
   * @param[in] matrix A ParalleMatrix object: the system matrix.
   * @param[out] rhs The rhs contribution to be modified
   * @param[in] bcValue The value to add to rhs
   * @param[in] fieldValue unused.
   *
   */
  GEOS_HOST_DEVICE
  static inline void
  SpecifyFieldValue( globalIndex const dof,
                     globalIndex const dofRankOffset,
                     CRSMatrixView< real64, globalIndex const > const & matrix,
                     real64 & rhs,
                     real64 const bcValue,
                     real64 const fieldValue )
  {
    GEOS_UNUSED_VAR( dof );
    GEOS_UNUSED_VAR( dofRankOffset );
    GEOS_UNUSED_VAR( matrix );
    GEOS_UNUSED_VAR( fieldValue );
    rhs += bcValue;
  }

  /**
   * @brief Function to add some values of a vector.
   * @tparam POLICY the execution policy to use when setting values
   * @param rhs the target right-hand side vector
   * @param dof a list of global DOF indices to be set
   * @param dofRankOffset offset of dof indices on current rank
   * @param values a list of values corresponding to \p dof that will be added to \p rhs.
   */
  template< typename POLICY >
  static inline void prescribeRhsValues( arrayView1d< real64 > const & rhs,
                                         arrayView1d< globalIndex const > const & dof,
                                         globalIndex const dofRankOffset,
                                         arrayView1d< real64 const > const & values )
  {
    GEOS_ASSERT_EQ( dof.size(), values.size() );
    forAll< POLICY >( dof.size(), [rhs, dof, dofRankOffset, values] GEOS_HOST_DEVICE ( localIndex const a )
    {
      globalIndex const localRow = dof[ a ] - dofRankOffset;
      if( localRow >= 0 && localRow < rhs.size() )
      { rhs[ localRow ] += values[ a ]; }
    } );
  }

};

} //namespace geos

#endif //GEOS_COMMON_FIELDSPECIFICATIONOPS_HPP
