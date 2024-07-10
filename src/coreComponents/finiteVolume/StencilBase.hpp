/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StencilBase.hpp
 */

#ifndef GEOS_FINITEVOLUME_STENCILBASE_HPP_
#define GEOS_FINITEVOLUME_STENCILBASE_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{

/**
 * @brief A collection of properties of a stencil type.
 * @tparam CONTAINER type of container used to store indices and weights
 * @tparam MAX_NUM_POINTS_IN_FLUX maximum number of points connected by a flux
 * @tparam MAX_STENCIL_SIZE maximum number of points in a stencil
 * @tparam MAX_NUM_CONNECTIONS maximum number of connections in a stencil
 */
template< template< typename ... > class CONTAINER,
          localIndex MAX_NUM_POINTS_IN_FLUX,
          localIndex MAX_STENCIL_SIZE,
          localIndex MAX_NUM_CONNECTIONS >
struct StencilTraits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = CONTAINER< localIndex >;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = LvArray::typeManipulation::NestedViewTypeConst< IndexContainerType >;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = CONTAINER< real64 >;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = LvArray::typeManipulation::NestedViewTypeConst< WeightContainerType >;

  /// The array view to type for the stencil weights
  using WeightContainerViewType = LvArray::typeManipulation::NestedViewType< WeightContainerType >;

  /// Maximum number of points the flux
  static constexpr localIndex maxNumPointsInFlux = MAX_NUM_POINTS_IN_FLUX;

  /// Maximum number of points in a stencil
  static constexpr localIndex maxStencilSize = MAX_STENCIL_SIZE;

  /// Maximum number of connections in a stencil
  static constexpr localIndex maxNumConnections = MAX_NUM_CONNECTIONS;
};

/**
 * @brief Describes properties of a standard two-point stencil.
 */
using TwoPointStencilTraits = StencilTraits< array2d, 2, 2, 1 >;

/**
 * @class StencilWrapperBase
 *
 * Class to provide access to the computation of stencil weights that may be
 * called from a kernel function.
 */
template< typename TRAITS >
class StencilWrapperBase : public TRAITS
{
public:

  /**
   * @brief Constructor
   * @param elementRegionIndices The container for the element region indices for each point in each stencil
   * @param elementSubRegionIndices The container for the element sub region indices for each point in each stencil
   * @param elementIndices The container for the element indices for each point in each stencil
   * @param weights The container for the weights for each point in each stencil
   */
  StencilWrapperBase( typename TRAITS::IndexContainerType const & elementRegionIndices,
                      typename TRAITS::IndexContainerType const & elementSubRegionIndices,
                      typename TRAITS::IndexContainerType const & elementIndices,
                      typename TRAITS::WeightContainerType const & weights ):
    m_elementRegionIndices( elementRegionIndices.toViewConst() ),
    m_elementSubRegionIndices( elementSubRegionIndices.toViewConst() ),
    m_elementIndices( elementIndices.toViewConst() ),
    m_weights( weights.toView() )
  {};

  /**
   * @brief Const access to the element regions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementRegionIndices() const { return m_elementRegionIndices; }

  /**
   * @brief Const access to the element subregions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementSubRegionIndices() const { return m_elementSubRegionIndices; }

  /**
   * @brief Const access to the element indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementIndices() const { return m_elementIndices; }

  /**
   * @brief Const access to the stencil weights.
   * @return A view to const
   */
  typename TRAITS::WeightContainerViewConstType
  getWeights() const { return m_weights; }

protected:

  /// The container for the element region indices for each point in each stencil
  typename TRAITS::IndexContainerViewConstType m_elementRegionIndices;

  /// The container for the element sub region indices for each point in each stencil
  typename TRAITS::IndexContainerViewConstType m_elementSubRegionIndices;

  /// The container for the element indices for each point in each stencil
  typename TRAITS::IndexContainerViewConstType m_elementIndices;

  /// The container for the weights for each point in each stencil
  typename TRAITS::WeightContainerViewType m_weights;
};


/**
 * @brief Provides management of the interior stencil points when using Two-Point flux approximation.
 * @tparam TRAITS the traits class describing properties of the stencil
 * @tparam LEAFCLASS derived type for CRTP
 */
template< typename TRAITS, typename LEAFCLASS >
class StencilBase : public TRAITS
{
public:

  /**
   * @brief Destructor.
   */
  virtual ~StencilBase() = default;

  /**
   * @brief Reserve the size of the stencil.
   * @param[in] size the size of the stencil to reserve
   */
  virtual void reserve( localIndex const size );

  /**
   * @brief Move the data arrays associated with the stencil to a specified
   *   memory space.
   * @param space The target memory space.
   *
   * @note The existence of this function indicates we need to redesign the
   * stencil classes.
   */
  virtual void move( LvArray::MemorySpace const space );


  /**
   * @brief Add an entry to the stencil.
   * @param[in] numPts The number of points in the stencil entry
   * @param[in] elementRegionIndices The element region indices for each point in the stencil entry
   * @param[in] elementSubRegionIndices The element sub-region indices for each point in the stencil entry
   * @param[in] elementIndices The element indices for each point in the stencil entry
   * @param[in] weights The weights each point in the stencil entry
   * @param[in] connectorIndex The index of the connector element that the stencil acts across
   */
  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) = 0;

  /**
   * @brief Zero weights for a stencil entry.
   * @param[in] connectorIndex The index of the connector element that the stencil acts across for which the weights are
   *                           to be zero.
   * @return True if a valid connectorIndex was found, and had its corresponding weights set to zero.
   */
  virtual bool zero( localIndex const connectorIndex );

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const = 0;

  /**
   * @brief Set the name used in data movement logging callbacks.
   * @param name the name prefix for the stencil's data arrays
   */
  void setName( string const & name );

  /**
   * @brief Const access to the element regions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementRegionIndices() const { return m_elementRegionIndices.toViewConst(); }

  /**
   * @brief Const access to the element subregions indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementSubRegionIndices() const { return m_elementSubRegionIndices.toViewConst(); }

  /**
   * @brief Const access to the element indices.
   * @return A view to const
   */
  typename TRAITS::IndexContainerViewConstType
  getElementIndices() const { return m_elementIndices.toViewConst(); }

  /**
   * @brief Const access to the stencil weights.
   * @return A view to const
   */
  typename TRAITS::WeightContainerViewConstType
  getWeights() const { return m_weights.toViewConst(); }

protected:

  /// The container for the element region indices for each point in each stencil
  typename TRAITS::IndexContainerType m_elementRegionIndices;

  /// The container for the element sub region indices for each point in each stencil
  typename TRAITS::IndexContainerType m_elementSubRegionIndices;

  /// The container for the element indices for each point in each stencil
  typename TRAITS::IndexContainerType m_elementIndices;

  /// The container for the weights for each point in each stencil
  typename TRAITS::WeightContainerType m_weights;

  /// The map that provides the stencil index given the index of the underlying connector object.
  unordered_map< localIndex, localIndex > m_connectorIndices;
};



template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 2 );
  m_elementSubRegionIndices.reserve( size * 2 );
  m_elementIndices.reserve( size * 2 );
  m_weights.reserve( size * 2 );
}


template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
bool StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::zero( localIndex const connectorIndex )
{
  return
    executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    for( localIndex i = 0; i < static_cast< LEAFCLASS * >(this)->stencilSize( connectorIndex ); ++i )
    {
      m_weights[connectionListIndex][i] = 0;
    }
  } );
}

template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::setName( string const & name )
{
  m_elementRegionIndices.setName( name + "/elementRegionIndices" );
  m_elementSubRegionIndices.setName( name + "/elementSubRegionIndices" );
  m_elementIndices.setName( name + "/elementIndices" );
  m_weights.setName( name + "/weights" );
}

template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::move( LvArray::MemorySpace const space )
{
  m_elementRegionIndices.move( space, true );
  m_elementSubRegionIndices.move( space, true );
  m_elementIndices.move( space, true );
  m_weights.move( space, true );
}

} /* namespace geos */

#endif /* GEOS_FINITEVOLUME_STENCILBASE_HPP_ */
