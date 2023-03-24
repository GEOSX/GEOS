/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_WELLBLOCK_HPP
#define GEOSX_WELLBLOCK_HPP

#include "mesh/generators/WellBlockABC.hpp"
#include "mesh/generators/InternalWellGenerator.hpp"


namespace geosx
{

/**
 * Implementation of the WellBlock responsible for modification/creation capabilities.
 */
class WellBlock : public WellBlockABC
{
public:
  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  WellBlock( string const & name, Group * const parent );

  /**
   * @name Getters / Setters
   */
  ///@{

  // getters for element data

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  globalIndex numElements() const override final { return m_numElems; }

  /**
   * @brief Set the global number of well elements.
   * @param numElems the global number of elements
   */
  void setNumElements( globalIndex numElems )  { m_numElems = numElems; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > getElemCoords() const override final { return m_elemCenterCoords; }

  /**
   * @brief Set the physical location of the centers of well elements.
   * @param elemCenterCoords  list of center locations of the well elements
   */
  void setElemCoords( arrayView2d< real64 const > elemCenterCoords )  { m_elemCenterCoords = elemCenterCoords; }

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  arrayView1d< globalIndex const > getNextElemIndex() const override final { return m_nextElemId; }

  /**
   * @brief Set the global indices mapping an element to the next.
   * @param nextElemId list providing the global index of the next element for each element
   */
  void setNextElemIndex( arrayView1d< globalIndex const > nextElemId )  { m_nextElemId = nextElemId; }

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  arrayView1d< arrayView1d< globalIndex const > const > getPrevElemIndices() const override final { return m_prevElemId.toNestedViewConst(); }


  /**
   * @brief Set the global indices mapping an element to the previous ones.
   * @param prevElemId list providing the global indices of the previous elements for each element
   */
  void setPrevElemIndices( arrayView1d< arrayView1d< globalIndex const > const > prevElemIndices );

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  arrayView2d< globalIndex const > getElemToNodesMap() const override final { return m_elemToNodesMap; }


  /**
   * @brief Set the global indices of the well nodes nodes connected to each element.
   * @param elemToNodesMap list providing the global index of the well nodes for each well element
   */
  void setElemToNodesMap( arrayView2d< globalIndex const > elemToNodesMap ) { m_elemToNodesMap = elemToNodesMap; }

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  arrayView1d< real64 const > getElemVolume() const override final { return m_elemVolume; }


  /**
   * @brief Set the volume of the well elements.
   * @param elemVolume list of volumes of the well elements
   */
  void setElemVolume( arrayView1d< real64 const > elemVolume ) { m_elemVolume = elemVolume; }

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  real64 getElementRadius() const override final { return m_radius; }


  /**
   * @brief Set the radius in the well.
   * @param radius the radius in the well
   */
  void setElementRadius( real64 radius ) { m_radius = radius; }

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  globalIndex numNodes() const override final { return m_numNodes; }


  /**
   * @brief Set the global number of well nodes.
   * @param numNodes the global number of nodes
   */
  void setNumNodes( globalIndex numNodes ) { m_numNodes = numNodes; }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView2d< real64 const > getNodeCoords() const override final { return m_nodeCoords; }

  /**
   * @brief Set the physical location of the centers of well elements.
   * @param nodeCoords list of center locations of the well elements
   */
  void setNodeCoords( arrayView2d< real64 const > nodeCoords ) { m_nodeCoords = nodeCoords; }



  // getters for perforation data

  /**
   * @brief Get the global number of perforations on this well.
   * @return the global number of elements
   */
  globalIndex numPerforations() const override final { return m_numPerforations; }


  /**
   * @brief Set the global number of perforations on this well.
   * @param numPerforations the global number of elements
   */
  void setNumPerforations( globalIndex numPerforations ) { m_numPerforations = numPerforations; }

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  arrayView2d< real64 const > getPerfCoords() const override final { return m_perfCoords; }


  /**
   * @brief Set the locations of the perforations.
   * @param  perfCoords list of locations of all the perforations on the well
   */
  void setPerfCoords( arrayView2d< real64 const > perfCoords ) { m_perfCoords = perfCoords; }

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  arrayView1d< real64 const > getPerfTransmissibility() const override final { return m_perfTransmissibility; }


  /**
   * @brief Set the well transmissibility at the perforations.
   * @param perfTransmissibility list of well transmissibility at all the perforations on the well
   */
  void setPerfTransmissibility( arrayView1d< real64 const > perfTransmissibility ) { m_perfTransmissibility = perfTransmissibility; }

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  arrayView1d< globalIndex const > getPerfElemIndex() const override final { return m_perfElemId; }

  /**
   * @brief Set the global indices of the well elements connected to each perforation.
   * @param perfElemId list providing the global index of the connected well element for each perforation
   */
  void setPerfElemIndex( arrayView1d< globalIndex const > perfElemId ) { m_perfElemId = perfElemId; }

  ///@}

  /// @endcond



private:


  // XML Input

  /// Radius area of the well (assumed to be valid for the entire well)
  real64 m_radius;

  // Geometry of the well (later passed to the WellElementSubRegion)

  // well element data

  /// Global number of well elements
  globalIndex m_numElems;

  /// Physical location of the center of the well element
  array2d< real64 > m_elemCenterCoords;

  /// Global index of the next well element
  array1d< globalIndex > m_nextElemId;

  /// Global indices of the prev well elements (maybe need multiple prevs for branching)
  array1d< array1d< globalIndex > > m_prevElemId;

  /// Connectivity between elements and nodes
  array2d< globalIndex > m_elemToNodesMap;

  /// Volume of well elements
  array1d< real64 > m_elemVolume;


  // well node data

  /// Number of nodes per well element
  globalIndex m_numNodesPerElem;

  /// Global number of well nodes
  globalIndex m_numNodes;

  /// Physical location of the nodes
  array2d< real64 > m_nodeCoords;

  // perforation data

  /// Global number of perforations
  globalIndex m_numPerforations;

  /// Absolute physical location of the perforation
  array2d< real64 > m_perfCoords;

  /// Well Peaceman index at the perforation
  array1d< real64 > m_perfTransmissibility;

  /// Global index of the well element
  array1d< globalIndex > m_perfElemId;



  // Auxiliary data

  // Number of physical dimensions
  int m_nDims;

  // Perforation data

  /// List of perforation names
  string_array m_perforationList;
};
}
#endif
