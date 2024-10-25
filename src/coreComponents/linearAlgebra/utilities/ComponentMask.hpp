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
 * @file ComponentMask.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_COMPONENTMASK_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_COMPONENTMASK_HPP_

#include "common/GeosxMacros.hpp"
#include "common/logger/Logger.hpp"

#include <cstdint>

namespace geos
{

namespace internal
{

/// Helps make static_assert's condition dependent on a non-type template parameter.
template< typename T >
constexpr bool always_false( T )
{
  return false;
}

template< int N >
struct ComponentMaskType
{
  static_assert( always_false( N ), "Unsupported template parameter N for ComponentMask. Use a value between 8 and 64." );
};

/**
 * @brief Macro for declaring a specialization of ComponentMaskTypeImpl template.
 */
#define GEOS_LAI_COMPONENTMASK_DECL( N ) \
  template<> \
  struct ComponentMaskType< N > \
  { \
    using type = std::uint ## N ## _t; \
  }

GEOS_LAI_COMPONENTMASK_DECL( 8 );
GEOS_LAI_COMPONENTMASK_DECL( 16 );
GEOS_LAI_COMPONENTMASK_DECL( 32 );
GEOS_LAI_COMPONENTMASK_DECL( 64 );

#undef GEOS_LAI_COMPONENTMASK_DECL

constexpr std::uint32_t roundToNextPowerOfTwo( std::uint32_t v )
{
  if( v <= 8 ) return 8;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

}

/**
 * @brief Utility class that represents a mask for included/excluded component of a mask.
 * @tparam MAX_COMP maximum number of components
 *
 * The main advantage over std::bitset<N> is being CUDA-friendly. In addition, it provides
 * convenience constructors and iteration over set bits.
 *
 * Class intended to be used by DofManager to specify component mappings when constructing
 * restriction/prolongation operators and during field-vector data copying.
 */
template< int MAX_COMP >
class ComponentMask
{
private:

  /// Number of bits in mask storage
  static constexpr int NUM_BITS = internal::roundToNextPowerOfTwo( MAX_COMP );

  /// Type used to represent the bit mask
  using mask_t = typename internal::ComponentMaskType< NUM_BITS >::type;

public:

  /**
   * @name Constructors.
   */
  ///@{

  /**
   * @brief Constructor.
   * @param numComp total number of components
   * @param includeAll if @p true, initializes the mask to include rather than exclude all components
   */
  GEOS_HOST_DEVICE
  ComponentMask( int const numComp, bool const includeAll )
    : m_mask( includeAll ? -1 : 0 )
  {
    setNumComp( numComp );
  }

  /**
   * @brief Constructor.
   * @param numComp total number of components
   */
  GEOS_HOST_DEVICE
  explicit ComponentMask( int const numComp = 0 )
    : ComponentMask( numComp, false )
  {}

  /**
   * @brief Range constructor, includes contiguous range of components
   * @param numComp total number of components
   * @param lo first component number to include
   * @param hi one-past last component number to include
   * @note [lo,hi) form a half-open range, indexing is zero-based.
   */
  GEOS_HOST_DEVICE
  ComponentMask( int const numComp, int const lo, int const hi )
    : ComponentMask( numComp )
  {
    for( int i = lo; i < hi; ++i )
    {
      set( i );
    }
  }

  /**
   * @brief Converting constructor from another mask with fewer max components.
   * @tparam M number of components in source mask
   * @param other source mask
   */
  template< int N, typename = std::enable_if_t< ( N < MAX_COMP ) > >
  GEOS_HOST_DEVICE
  ComponentMask( ComponentMask< N > const & other )
    : m_numComp( other.numComp() )
  {
    for( int i : other )
    {
      set( i );
    }
  }

  ///@}

  /**
   * @name Inspection methods.
   */
  ///@{

  /**
   * @brief Forward const iterator over the set components (bits).
   * @note Iterator is invalidated by any mutating operation on the ComponentMask object.
   */
  class Iterator
  {
public:

    using difference_type = int;                         ///< difference type
    using value_type = int;                              ///< value type
    using pointer = void;                                ///< pointer type (meaningless but makes it iterator_traits compatible)
    using reference = int;                               ///< reference type (no real references since this is a const iterator)
    using iterator_category = std::forward_iterator_tag; ///< iterator category

    /**
     * @brief Dereference operator.
     * @return index of currently pointed-to component
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    reference operator*() const
    {
      return m_index;
    }

    /**
     * @brief Prefix increment operator.
     * @return reference to @p this
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    Iterator & operator++()
    {
      GEOS_ASSERT_NE( m_mask, 0 );
      skipOne();
      findNextBit();
      return *this;
    }

    /**
     * @brief Postfix increment operator.
     * @return copy of @p this prior to increment
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    Iterator operator++( int ) &
    {
      Iterator old = *this;
      ++*this;
      return old;
    }

    /**
     * @brief Comparison operator.
     * @param other the iterator to compare to
     * @return @p true if iterators are equal (as determined by remaining mask)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    bool operator==( Iterator const & other ) const
    {
      return m_mask == other.m_mask;
    }

    /**
     * @brief Comparison operator.
     * @param other the iterator to compare to
     * @return @p true if iterators are not equal
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    bool operator!=( Iterator const & other ) const
    {
      return !(*this == other);
    }

private:

    GEOS_HOST_DEVICE
    inline
    void skipOne()
    {
      m_mask >>= 1;
      ++m_index;
    }

    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    void findNextBit()
    {
      if( m_mask )
      {
        // TODO: use hardware intrinsics (ffs/ctz) where possible
        while( !( m_mask & 1 ) )
        {
          skipOne();
        }
      }
    }

    GEOS_HOST_DEVICE
    Iterator( mask_t const mask, int const index )
      : m_mask( mask ),
      m_index( index )
    {
      findNextBit();
    }

    friend class ComponentMask;

    mask_t m_mask;
    int m_index;
  };

  /**
   * @brief @return iterator to the start of the component index range
   */
  GEOS_HOST_DEVICE
  Iterator begin() const
  {
    return Iterator( m_mask, 0 );
  }

  /**
   * @brief @return iterator past the end of the component index range
   */
  GEOS_HOST_DEVICE
  Iterator end() const
  {
    return Iterator( 0, -1 );
  }

  /**
   * @brief @return the number of set components.
   */
  GEOS_HOST_DEVICE
  int size() const
  {
    mask_t mask = m_mask;
    int c = 0;
    for(; mask != 0; ++c )
    {
      mask &= mask - 1;
    }
    return c;
  }

  /**
   * @brief @return number of components
   */
  GEOS_HOST_DEVICE
  int numComp() const
  {
    return m_numComp;
  }

  /**
   * @brief @return @p true if mask is empty (no components selected), @p false otherwise.
   */
  GEOS_HOST_DEVICE
  bool empty() const
  {
    return m_mask == 0;
  }

  /**
   * @brief Check if a particular component is selected.
   * @param i component index
   * @return @p true if i-th component is selected, @p false otherwise
   */
  GEOS_HOST_DEVICE
  bool operator[]( int const i ) const
  {
    GEOS_ASSERT_GT( m_numComp, i );
    return m_mask & (mask_t( 1 ) << i);
  }

  ///@}

  /**
   * @name Modification methods/operators.
   */
  ///@{

  /**
   * @brief Clear the mask, removing selected components.
   */
  GEOS_HOST_DEVICE
  void clear()
  {
    m_mask = 0;
  }

  /**
   * @brief Set new component limit and truncate all components above.
   * @param numComp new max components value
   */
  GEOS_HOST_DEVICE
  void setNumComp( int const numComp )
  {
    GEOS_ASSERT_GE( numComp, 0 );
    GEOS_ASSERT_GE( MAX_COMP, numComp );
    m_numComp = numComp;
    if( numComp < NUM_BITS )
    {
      m_mask &= ( mask_t( 1 ) << numComp ) - 1;
    }
  }

  /**
   * @brief Add a component.
   * @param i component index
   */
  GEOS_HOST_DEVICE
  void set( int const i )
  {
    GEOS_ASSERT_GE( i, 0 );
    GEOS_ASSERT_GT( m_numComp, i );
    m_mask |= mask_t( 1 ) << i;
  }

  /**
   * @brief Remove a component.
   * @param i component index
   */
  GEOS_HOST_DEVICE
  void unset( int const i )
  {
    GEOS_ASSERT_GE( i, 0 );
    GEOS_ASSERT_GT( m_numComp, i );
    m_mask &= ~(mask_t( 1 ) << i);
  }

  /**
   * @brief Invert component selection.
   */
  GEOS_HOST_DEVICE
  void invert()
  {
    m_mask = ~m_mask;
    setNumComp( m_numComp );
  }

  /**
   * @brief Add a component.
   * @param i component index
   * @return reference to this object
   */
  GEOS_HOST_DEVICE
  ComponentMask & operator+=( int const i )
  {
    set( i );
    return *this;
  }

  /**
   * @brief Remove a component.
   * @param i component index
   * @return reference to this object
   */
  GEOS_HOST_DEVICE
  ComponentMask & operator-=( int const i )
  {
    unset( i );
    return *this;
  }

  /**
   * @brief Negation operator (inverts component selection).
   * @return new mask that is the inversion of this object
   */
  GEOS_HOST_DEVICE
  ComponentMask operator~() const
  {
    ComponentMask mask( *this );
    mask.invert();
    return mask;
  }

  /**
   * @brief Make a new mask by adding a component.
   * @param mask the mask to append to
   * @param i component index
   * @return new component mask
   */
  GEOS_HOST_DEVICE
  friend ComponentMask operator+( ComponentMask const & mask, int const i )
  {
    ComponentMask new_mask = mask;
    new_mask += i;
    return new_mask;
  }

  /**
   * @brief Make a new mask by adding a component.
   * @param i component index
   * @param mask the mask to append to
   * @return new component mask
   */
  GEOS_HOST_DEVICE
  friend ComponentMask operator+( int const i, ComponentMask const & mask )
  {
    return mask + i;
  }

  /**
   * @brief Make a new mask by removing a component.
   * @param mask the mask to remove from
   * @param i component index
   * @return new component mask
   */
  GEOS_HOST_DEVICE
  friend ComponentMask operator-( ComponentMask const & mask, int const i )
  {
    ComponentMask new_mask = mask;
    new_mask -= i;
    return new_mask;
  }

  /**
   * @brief Make a new mask by removing a component.
   * @param i component index
   * @param mask the mask to remove from
   * @return new component mask
   */
  GEOS_HOST_DEVICE
  friend ComponentMask operator-( int const i, ComponentMask const & mask )
  {
    return mask - i;
  }

  ///@}

private:

  /// Bit representation of the component mask
  mask_t m_mask = 0;

  /// Runtime number of components (included or excluded) represented by this mask
  int m_numComp = MAX_COMP;
};

}

#endif //GEOS_LINEARALGEBRA_UTILITIES_COMPONENTMASK_HPP_
