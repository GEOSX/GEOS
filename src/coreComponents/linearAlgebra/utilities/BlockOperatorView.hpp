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
 * @file BlockOperatorView.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_

#include "codingUtilities/SFINAE_Macros.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockVector.hpp"
#include "common/traits.hpp"
#include "linearAlgebra/common/common.hpp"

namespace geos
{

/**
 * @brief Abstract view of a block operator.
 * @tparam VECTOR type of vector that sub-blocks of this view can operate on
 * @tparam OPERATOR type of operator that can operate on @p VECTOR
 *                  (can be base class or a more specialized derived class)
 *
 * This class does not deal with constructing or storing sub-blocks, only provides high-level access functions.
 * See derived classes BlockOperator and BlockOperatorWrapper for ways to construct a block operator.
 */
template< typename VECTOR, typename OPERATOR = LinearOperator< VECTOR > >
class BlockOperatorView : public LinearOperator< BlockVectorView< VECTOR > >
{

public:

  /// Alias for base type
  using Base = LinearOperator< BlockVectorView< VECTOR > >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /**
   * @name Constructors/destructors
   */
  ///@{

  /**
   * @brief Destructor.
   */
  virtual ~BlockOperatorView() override = default;

  /**
   * @brief Deleted copy assignment.
   * @return not callable
   */
  BlockOperatorView & operator=( BlockOperatorView const & ) = delete;

  /**
   * @brief Deleted move assignment.
   * @return not callable
   */
  BlockOperatorView & operator=( BlockOperatorView && ) = delete;

  ///@}

  /**
   * @name LinearOperator interface
   */
  ///@{

  /**
   * @brief Apply operator to a vector
   * @param src Input vector (x).
   * @param dst Output vector (b).
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override
  {
    for( localIndex i = 0; i < numBlockRows(); i++ )
    {
      dst.block( i ).zero();
      VECTOR temp( dst.block( i ) );
      for( localIndex j = 0; j < numBlockCols(); j++ )
      {
        if( m_operators( i, j ) != nullptr )
        {
          m_operators( i, j )->apply( src.block( j ), temp );
          dst.block( i ).axpy( 1.0, temp );
        }
      }
    }
  }

  /**
   * @brief Apply the transpose of block operator to a block vector.
   * @tparam OP dummy template parameter to enable SFINAE, do not provide
   * @param src source vector (rhs)
   * @param dst target vector (lhs)
   * @return nothing
   *
   * @note This method only exists if the underlying operator type has it.
   */
  template< typename OP = OPERATOR >
  std::enable_if_t< traits::VectorBasedTraits< Vector >::template HasMemberFunction_applyTranspose< OP > >
  applyTranspose( BlockVectorView< VECTOR > const & src,
                  BlockVectorView< VECTOR > & dst ) const
  {
    for( localIndex j = 0; j < numBlockCols(); j++ )
    {
      dst.block( j ).zero();
      VECTOR temp( dst.block( j ) );
      for( localIndex i = 0; i < numBlockRows(); i++ )
      {
        if( m_operators( j, i ) != nullptr )
        {
          m_operators( j, i )->applyTranspose( src.block( i ), temp );
          dst.block( j ).axpy( 1.0, temp );
        }
      }
    }
  }

  /**
   * @brief Get the number of global rows.
   * @return Number of global rows in the operator.
   */
  virtual globalIndex numGlobalRows() const override
  {
    return computeRowSize( []( OPERATOR const & block ) { return block.numGlobalRows(); } );
  }

  /**
   * @brief Get the number of global columns.
   * @return Number of global columns in the operator.
   */
  virtual globalIndex numGlobalCols() const override
  {
    return computeColSize( []( OPERATOR const & block ) { return block.numGlobalCols(); } );
  }

  /**
   * @brief Get the number of local rows.
   * @return Number of local rows in the operator.
   * @note Method only exists if the underlying operator type has it.
   */
  virtual localIndex numLocalRows() const override
  {
    return computeRowSize( []( OPERATOR const & block ) { return block.numLocalRows(); } );
  }

  /**
   * @brief Get the number of local columns.
   * @return Number of local columns in the operator.
   * @note Method only exists if the underlying operator type has it.
   */
  virtual localIndex numLocalCols() const override
  {
    return computeColSize( []( OPERATOR const & block ) { return block.numLocalCols(); } );
  }

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator passed in @p create...()
   *
   * @note when build without MPI, may return anything
   *       (MPI_Comm will be a mock type defined in MpiWrapper)
   */
  virtual MPI_Comm comm() const override
  {
    return m_operators( 0, 0 )->comm();
  }

  ///@}

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get number of block rows.
   * @return number of block rows
   */
  localIndex numBlockRows() const
  {
    return m_operators.size( 0 );
  }

  /**
   * @brief Get number of block columns.
   * @return number of block columns
   */
  localIndex numBlockCols() const
  {
    return m_operators.size( 1 );
  }

  /**
   * @brief Get the operator corresponding to a sub-block.
   * @param blockRowIndex block row index
   * @param blockColIndex block column index
   * @return a reference to the sub-block
   */
  OPERATOR const & block( localIndex const blockRowIndex, localIndex const blockColIndex ) const
  {
    GEOS_LAI_ASSERT( m_operators( blockRowIndex, blockColIndex ) != nullptr );
    return *m_operators( blockRowIndex, blockColIndex );
  }

  /**
   * @copydoc block( localIndex const, localIndex const ) const
   */
  OPERATOR & block( localIndex const blockRowIndex, localIndex const blockColIndex )
  {
    GEOS_LAI_ASSERT( m_operators( blockRowIndex, blockColIndex ) != nullptr );
    return *m_operators( blockRowIndex, blockColIndex );
  }

  /**
   * @copydoc block( localIndex const, localIndex const ) const
   */
  OPERATOR const & operator()( localIndex const blockRowIndex, localIndex const blockColIndex = 0 ) const
  {
    return block( blockRowIndex, blockColIndex );
  }

  /**
   * @copydoc block( localIndex const, localIndex const ) const
   */
  OPERATOR & operator()( localIndex const blockRowIndex, localIndex const blockColIndex = 0 )
  {
    return block( blockRowIndex, blockColIndex );
  }

  ///@}

protected:

  /**
   * @brief Create an operator with given number of block rows and columns.
   * @param nRows number of block rows
   * @param nCols number of block columns
   */
  BlockOperatorView( localIndex const nRows, localIndex const nCols )
    : m_operators( nRows, nCols )
  {
    GEOS_LAI_ASSERT_GT( nRows, 0 );
    GEOS_LAI_ASSERT_GT( nCols, 0 );
  }

  /**
   * @brief Copy constructor.
   * @param x the block vector to copy
   */
  BlockOperatorView( BlockOperatorView< VECTOR, OPERATOR > const & x ) = default;

  /**
   * @brief Move constructor.
   * @param x the block vector to move from
   */
  BlockOperatorView( BlockOperatorView< VECTOR, OPERATOR > && x ) = default;

  /**
   * @brief Set/replace a pointer to a block.
   * @param blockRowIndex row index of the block
   * @param blockColIndex column index of the block
   * @param op            the new pointer
   */
  void setPointer( localIndex const blockRowIndex, localIndex const blockColIndex, OPERATOR * op )
  {
    GEOS_LAI_ASSERT_GE( blockRowIndex, 0 );
    GEOS_LAI_ASSERT_GT( numBlockRows(), blockRowIndex );
    GEOS_LAI_ASSERT_GE( blockColIndex, 0 );
    GEOS_LAI_ASSERT_GT( numBlockCols(), blockColIndex );
    m_operators( blockRowIndex, blockColIndex ) = op;
  }

private:

  template< typename FUNC >
  auto computeRowSize( FUNC func ) const -> decltype( func( std::declval< OPERATOR const >() ) );

  template< typename FUNC >
  auto computeColSize( FUNC func ) const -> decltype( func( std::declval< OPERATOR const >() ) );

  /// Array of pointers to blocks
  array2d< OPERATOR * > m_operators;
};

template< typename VECTOR, typename OPERATOR >
template< typename FUNC >
auto BlockOperatorView< VECTOR, OPERATOR >::computeRowSize( FUNC func ) const -> decltype( func( std::declval< OPERATOR const >() ) )
{
  using sizeType = decltype( func( std::declval< OPERATOR const >() ) );
  sizeType rowSize = 0;
  for( localIndex i = 0; i < numBlockRows(); ++i )
  {
    for( localIndex j = 0; j < numBlockCols(); ++j )
    {
      if( m_operators( i, j ) != nullptr )
      {
        rowSize += func( block( i, j ) );
        break;
      }
    }
  }
  return rowSize;
}

template< typename VECTOR, typename OPERATOR >
template< typename FUNC >
auto BlockOperatorView< VECTOR, OPERATOR >::computeColSize( FUNC func ) const -> decltype( func( std::declval< OPERATOR const >() ) )
{
  using sizeType = decltype( func( std::declval< OPERATOR const >() ) );
  sizeType colSize = 0;
  for( localIndex j = 0; j < numBlockCols(); j++ )
  {
    for( localIndex i = 0; i < numBlockRows(); ++i )
    {
      if( m_operators( i, j ) != nullptr )
      {
        colSize += func( block( i, j ) );
        break;
      }
    }
  }
  return colSize;
}

}// end geos namespace


#endif /*GEOS_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_*/
