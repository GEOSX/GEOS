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

#ifndef GEOS_FINITEVOLUME_PROJECTIONEDFMHELPER_HPP_
#define GEOS_FINITEVOLUME_PROJECTIONEDFMHELPER_HPP_

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{

class MeshLevel;
class CellElementStencilTPFA;
class EmbeddedSurfaceToCellStencil;

/**
 * @struct CellDescriptor
 * @brief A structure containing a single cell (element) identifier triplet.
 */
struct CellDescriptor
{
  /// region index
  localIndex region;
  /// subregion index
  localIndex subRegion;
  /// cell index
  localIndex index;

  /**
   * @brief Constructor for the CellDescriptor struct
   * @param[r] region index
   * @param[sr] subregion index
   * @param[i] cell index
   */
  CellDescriptor( localIndex r, localIndex sr, localIndex i )
    : region( r ), subRegion( sr ), index( i )
  {}

  /**
   * @brief Comparison operator between two CellDescriptors.
   * @param[in] other the CellDescriptor to compare with
   * @return true if they represent the same mesh element
   */
  bool operator==( CellDescriptor const & other ) const
  {
    return( region==other.region && subRegion==other.subRegion && index==other.index );
  }
};

/**
 * @class ProjectionEDFMHelper
 * @brief A class that contains methods to modify cell and edfm stencils based on projection edfm.
 */
class ProjectionEDFMHelper
{
public:

  /**
   * @brief Constructor
   * @param mesh the meshLevel object
   * @param stencil the cellStencil
   * @param edfmStencil the edfmStencil
   * @param embeddedSurfaceRegionName name of the embeddedSurfaceRegion
   */
  ProjectionEDFMHelper( MeshLevel const & mesh,
                        CellElementStencilTPFA & stencil,
                        EmbeddedSurfaceToCellStencil & edfmStencil,
                        string const & embeddedSurfaceRegionName );

  /**
   * @brief add Fracture-matrix connections to the edfmStencil
   * and remove the appropriate connections from the cellStencil
   */
  void addNonNeighboringConnections() const;

private:

  /**
   * @brief select cell faces that will host non-neighboring fracture-matrix connections
   * @param[in] subRegionFaces faces of the subRegion
   * @param[in] hostCellID id of the hostCell (tre fractured one)
   * @param[in] fracElement the index of the fracture element
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   * @return a list of the faces that need to be disconnected.
   */
  std::vector< localIndex > selectFaces( FixedOneToManyRelation const & subRegionFaces,
                                         CellDescriptor const & hostCellID,
                                         localIndex const fracElement,
                                         EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  /**
   * @brief check the intersection a fracture element and an edge
   * @param[in] fracCenter the coordinates of the center of the fracture cell
   * @param[in] fracNomral the normal vector of the fracture element
   * @param[in] edgeIdx the edge index
   * @return whether the fracture element intersects the edge or not.
   */
  bool intersection( real64 const ( &fracCenter )[3],
                     arraySlice1d< real64 const > const & fracNormal,
                     localIndex const edgeIdx ) const;

  /**
   * @brief returns true is the face has only one neighbor
   * @param[in] faceIdx index of the face
   * @return whether the face is a boundary face or not.
   */
  bool isBoundaryFace( localIndex const faceIdx ) const;

  /**
   * @brief check if the center of a face is on the same side of the fracture as cell center
   * @param[in] faceIdx index of the face
   * @param[in] signedDistanceCellCenterToFrac the signed  distance between the cell center and the fracture
   * @param[in] fracCenter the coordinate of the frac center (center of the embSurf element)
   * @param[in] fracNormal the normal to the fracture segment
   * @return a bool that defines whether the center of a face is on the same side of the fracture as cell center
   */
  bool onLargerSide( localIndex const faceIdx,
                     real64 const signedDistanceCellCenterToFrac,
                     real64 const ( &fracCenter )[3],
                     arraySlice1d< real64 const > const & fracNormal ) const;

  /*
   * @brief compute the signed distance between the fracture and as cell center
   * @param[in] hostCellID id of the host cell
   * @param[in] fracNormal the normal vector of the fracture element
   * @param[in] fracCenter the coordinates of the center of the fracture cell
   * @param[out] cellCenterToFracCenter distance vector between the cell and the fracture.
   * return The signed distance between the fracture and the cell (based on the frac normal direction)
   */
  real64 getSignedDistanceCellCenterToFracPlane( CellDescriptor const & hostCellID,
                                                 arraySlice1d< real64 const > const & fracNormal,
                                                 real64 const (&fracCenter)[3],
                                                 real64 ( &cellCenterToFracCenter )[3] ) const;

  /*
   * @brief returns true if the signed distance from the neighbor center to the frac is of the same sign as signedDistanceCellCenterToFrac
   *(computed in the host cell)
   * @@param[in] faceIdx the face index
   * @param[in] signedDistanceCellCenterToFrac signed distance between the fracture and the cell (based on the frac normal direction)
   * @param[in] hostCellID id of the host cell (region, subregion and element index)
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   * return whether the neighboring cell has the fracture on the same side of the face or not.
   */
  bool neighborOnSameSide( localIndex const faceIdx,
                           real64 const signedDistanceCellCenterToFrac,
                           CellDescriptor const & hostCellID,
                           EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  /**
   * @brief given a face and its neighboring cell, return the id of the other cell
   * @param[in] faceIdx face index
   * @param[in] hostCellID id of the host cell
   * return the ide of the other cell throug a cell descriptor
   */
  CellDescriptor otherCell( localIndex const faceIdx, CellDescriptor const & hostCellID ) const;

  /**
   * @brief compute the absolute transmissibility for non-neighboring F-M connection
   * @param[in] neighborCell the neighboring cell
   * @param[in] fracElement the index of the embeddedSurface elmement (the fracture element)
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   * @param[in] faceIdx the face index
   * @param[out] trans The geometric transmissibility between fracture and matrix cells.
   */
  void computeFractureMatrixWeights( CellDescriptor const & neighborCell,
                                     localIndex const fracElement,
                                     EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                     localIndex const faceIdx,
                                     real64 ( &weights )[2] ) const;

  /**
   * @brief add non-neighboring F-M connection with given transmissibility tothe cell stencil
   * @param[in] fracElement fracture element index
   * @param[in] cell id of the cell
   * @param[in] transmissibility geometric transmissiblity
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   */
  void addNonNeighboringConnection( localIndex const fracElement,
                                    CellDescriptor const & cell,
                                    real64 const (&weights)[2],
                                    EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  // Private variables

  /// The element region manager
  ElementRegionManager const & m_elementManager;
  /// the nodes coordinates
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodesCoord;
  /// edge-to-nodes map
  arrayView2d< localIndex const > const m_edgeToNodes;
  /// face-to-regions map
  arrayView2d< localIndex const > const m_faceToRegions;
  /// face-to-subRegions map
  arrayView2d< localIndex const > const m_faceToSubRegions;
  /// face-to-element map
  arrayView2d< localIndex const > const m_faceToCells;
  /// face-to-edges map
  ArrayOfArraysView< localIndex const > m_faceToEdges;
  /// face-to-nodes map
  ArrayOfArraysView< localIndex const > const m_faceToNodes;
  /// view accessor to the cell centers
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const m_cellCenters;
  /// the cellStencil
  CellElementStencilTPFA & m_cellStencil;
  /// the pedfm stencil
  EmbeddedSurfaceToCellStencil & m_edfmStencil;
  ///
  string const m_embeddedSurfaceRegionName;

};

}  // end namespace geos


#endif // __PROJECTIONEDFMHELPER_H_
