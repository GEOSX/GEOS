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

/*
 * @file PerforationData.hpp
 */

#ifndef GEOS_MESH_PERFORATIONDATA_HPP
#define GEOS_MESH_PERFORATIONDATA_HPP

#include "dataRepository/Group.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/ToElementRelation.hpp"
#include "mesh/generators/LineBlockABC.hpp"

namespace geos
{

class DomainPartition;
class MeshLevel;
class WellElementSubRegion;

/**
 * @class PerforationData
 * This class keeps track of all the local perforations on this rank
 */
class PerforationData : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for PerforationData Objects.
   * @param[in] name name of this instantiation of PerforationData in the repository
   * @param[in] parent parent group of this instantiation of PerforationData
   */
  explicit PerforationData( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~PerforationData() override;

  /**
   * @brief Deleted default constructor.
   */
  PerforationData() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  PerforationData( PerforationData const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  PerforationData( PerforationData && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a PerforationData object
   */
  PerforationData & operator=( PerforationData const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a PerforationData object
   */
  PerforationData & operator=( PerforationData && ) = delete;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "PerforationData"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Set the global number of perforations used for well initialization.
   * @param[in] nPerfs global number of perforations (obtained from LineBlockABC)
   */
  void setNumPerforationsGlobal( globalIndex nPerfs ) { m_numPerforationsGlobal = nPerfs; }


  /**
   * @brief Get the global number of perforations (used for well initialization).
   * @return global number of perforations
   */
  globalIndex getNumPerforationsGlobal() const { return m_numPerforationsGlobal; }


  /**
   * @brief Get perforation-to-mesh-element connectivity.
   * @return list of element region/subregion/index connected to each perforation
   */
  ToElementRelation< array1d< localIndex > > & getMeshElements() { return m_toMeshElements; }


  /**
   * @brief Provide an immutable accessor to a const perforation-to-mesh-element connectivity.
   * @return list of element region/subregion/index connected to each perforation
   */
  ToElementRelation< array1d< localIndex > > const & getMeshElements() const { return m_toMeshElements; }


  /**
   * @brief Get perforation-to-well-element connectivity.
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex > getWellElements() { return m_wellElementIndex; }


  /**
   * @brief Provide an immutable accessor to a const perforation-to-well-element connectivity.
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex const > getWellElements() const { return m_wellElementIndex; }

  /**
   * @brief Get perforation-to-reservoir-element connectivity.
   * @return list of global reservoir element index connected to each perforation
   */
  arrayView1d< globalIndex > getReservoirElementGlobalIndex() { return m_reservoirElementGlobalIndex; }



  /**
   * @brief Provide an immutable accessor to a const perforation-to-reservoir-element connectivity.
   * @return list of well element index connected to each perforation
   */
  arrayView1d< globalIndex const > getReservoirElementGlobalIndex() const { return m_reservoirElementGlobalIndex; }


  /**
   * @brief Get perforation locations.
   * @return list of perforation locations
   */
  arrayView2d< real64 > getLocation() { return m_location; }


  /**
   * @brief Provide an immutable accessor to a const perforation location arrayView.
   * @return list of perforation locations
   */
  arrayView2d< real64 const > getLocation() const { return m_location; }


  /**
   * @brief Provide an immutable accessor to a const perforation well indices array.
   * @return list of perforation well indices
   */
  arrayView1d< real64 const > getWellTransmissibility() const { return m_wellTransmissibility; }


  /**
   * @brief Get perforation well indices.
   * @return list of perforation well indices
   */
  arrayView1d< real64 > getWellTransmissibility() { return m_wellTransmissibility; }


  /**
   * @brief Provide an immutable accessor to a const perforation skin factor array.
   * @return list of perforation well skin factors
   */
  arrayView1d< real64 const > getWellSkinFactor() const { return m_wellSkinFactor; }


  /**
   * @brief Get perforation well skin factors.
   * @return list of perforation well skin factors
   */
  arrayView1d< real64 > getWellSkinFactor() { return m_wellSkinFactor; }


  ///@}

  /**
   * @name Well transmissibility computation
   */
  ///@{

  /**
   * @brief Compute the well transmissibility for each local perforation on this well.
   * @param[in] mesh target mesh level
   * @param[in] wellElemSubRegion  subRegion corresponding to this well
   * @param[in] perm the permeability in the reservoir
   */
  void computeWellTransmissibility( MeshLevel const & mesh,
                                    WellElementSubRegion const & wellElemSubRegion,
                                    array1d< array1d< arrayView3d< real64 const > > > const & perm );

  ///@}

  /**
   * @name Construction of the connectivity
   */
  ///@{

  /**
   * @brief Connect each perforation to a local wellbore element.
   * @param[in] lineBlock LineBlockABC containing the global well topology
   * @param[in] globalToLocalWellElementMap  global-to-local map of wellbore elements
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */
  void connectToWellElements( LineBlockABC const & lineBlock,
                              unordered_map< globalIndex, localIndex > const & globalToLocalWellElementMap,
                              globalIndex elemOffsetGlobal );

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the global number of perforations
    static constexpr char const * numPerforationsGlobalString() { return "numPerforationsGlobal"; }
  }
  /// ViewKey struct for the PerforationData class
  viewKeysPerforationData;

private:

  /**
   * @name Helpers for transmissibility computation
   */
  ///@{

  /**
   * @brief Compute the approximate dimensions of the reservoir element containing a perforation.
   *        This is done by computing a bounding box containing the element.
   * @param[in] mesh target mesh level
   * @param[in] er  index of the element region containing the reservoir element
   * @param[in] esr  index of the element subRegion containing the reservoir element
   * @param[in] ei  index of the reservoir element
   * @param[out] dx dimension of the element in the x-direction
   * @param[out] dy dimension of the element in the y-direction
   * @param[out] dz dimension of the element in the z-direction
   */
  void getReservoirElementDimensions( MeshLevel const & mesh,
                                      localIndex const er, localIndex const esr, localIndex const ei,
                                      real64 & dx, real64 & dy, real64 & dz ) const;

  ///@}


  /// Global number of perforations
  globalIndex m_numPerforationsGlobal;

  /// Indices of the mesh elements connected to perforations
  ToElementRelation< array1d< localIndex > > m_toMeshElements;

  /// Indices of the well elements to which perforations are attached
  array1d< localIndex > m_wellElementIndex;

  /// Global indices of reservoir cell containing perforation
  array1d< globalIndex > m_reservoirElementGlobalIndex;

  /// Location of the perforations
  array2d< real64 > m_location;

  /// Well transmissibility at the perforations
  array1d< real64 > m_wellTransmissibility;

  /// Well skin factor at the perforations
  array1d< real64 > m_wellSkinFactor;

};

} //namespace geos

#endif //GEOS_MESH_PERFORATIONDATA_HPP
