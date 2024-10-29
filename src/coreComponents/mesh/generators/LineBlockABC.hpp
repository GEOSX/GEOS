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

#ifndef GEOS_WELLBLOCKABC_HPP
#define GEOS_WELLBLOCKABC_HPP

#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"

#include <vector>

namespace geos
{

/**
 * Abstract base class defining the information provided by any well block implementation.
 *
 * It's noteworthy that the WellblockABC is immutable oriented.
 * The derived implementations need to have the modification/creation capabilities.
 */
class LineBlockABC : public dataRepository::Group
{
public:

  /**
   * @brief Struct to define the top and bottom node of a segment.
   */
  struct NodeLocation
  {
    static constexpr integer TOP    = 0; /**< Top of the well */
    static constexpr integer BOTTOM = 1; /**< Bottom of the well */
  };

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  LineBlockABC( string const & name,
                Group * const parent )
    :
    Group( name, parent )
  { }

  /**
   * @name Getters / Setters
   */
  ///@{

  // getters for element data

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  virtual globalIndex numElements() const = 0;

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  virtual arrayView2d< real64 const > getElemCoords() const = 0;

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  virtual arrayView1d< globalIndex const > getNextElemIndex() const = 0;

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  virtual arrayView1d< arrayView1d< globalIndex const > const > getPrevElemIndices() const = 0;

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  virtual arrayView2d< globalIndex const > getElemToNodesMap() const = 0;

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  virtual arrayView1d< real64 const > getElemVolume() const = 0;

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  virtual real64 getElementRadius() const = 0;

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  virtual globalIndex numNodes() const = 0;

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  virtual arrayView2d< real64 const > getNodeCoords() const = 0;



  // getters for perforation data

  /**
   * @brief Get the global number of perforations on this well.
   * @return the global number of elements
   */
  virtual globalIndex numPerforations() const = 0;

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  virtual arrayView2d< real64 const > getPerfCoords() const = 0;

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  virtual arrayView1d< real64 const > getPerfTransmissibility() const = 0;

  /**
   * @brief Get the well skin factor at the perforations.
   * @return list of well skin factor at all the perforations on the well
   */
  virtual arrayView1d< real64 const > getPerfSkinFactor() const = 0;

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  virtual arrayView1d< globalIndex const > getPerfElemIndex() const = 0;

  /**
   * @brief Get the well controls name
   * @return The well controls name
   */
  virtual string const & getWellControlsName() const = 0;

  /**
   * @brief Get the well generator name
   * @return The well generator name
   */
  virtual string const & getWellGeneratorName() const = 0;
};

}

#endif //GEOS_WELLBLOCKABC_HPP
