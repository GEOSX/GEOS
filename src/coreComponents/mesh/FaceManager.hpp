/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FaceManager.hpp
 */

#ifndef GEOSX_MESH_FACEMANAGER_HPP_
#define GEOSX_MESH_FACEMANAGER_HPP_

#include "ToElementRelation.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class NodeManager;
class ElementRegionManager;
class CellElementSubRegion;


/**
 * @class FaceManager
 * @brief The FaceManager class provides an interface to ObjectManagerBase in order to manage face data.
 *
 * The FaceManager class manages the face data using face indexed or keyed containers.
 * This means that each field is stored in a way where each container entry
 * corresponds to a face.
 */
class FaceManager : public ObjectManagerBase
{
public:

  using NodeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;
  using ElemMapType = FixedToManyElementRelation;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Return the name of the face manager in the object catalog.
   * @return string that contains the catalog name to generate a new FaceManager object through the object catalog.
   */
  static const string CatalogName()
  { return "FaceManager"; }

  
  /** @brief Provide a virtual access to CatalogName().
   * @return string that contains the catalog name to generate a new FaceManager object through the object catalog.
   */
  virtual const string getCatalogName() const override
  { return FaceManager::CatalogName(); }
 ///@}


  
  /**
     @brief Get the default number of node per face in node list
     @return the default number of node per face in node list
   */
  static localIndex nodeMapExtraSpacePerFace()
  { return 4; }
  
  
  /**
     @brief Get the default number of edge per face in edge list
     @return the default number of edge per face in edge list
   */
  static localIndex edgeMapExtraSpacePerFace()
  { return 4; }




   /**
    *  @name Constructors/destructor
    */

  ///@{

  /**
   *  @brief Main Constructor for FaceManager
   *  @param[in] name the name of FaceManager
   *  @param[in] parent the parent Group object of FaceManager
   */
  FaceManager( string const & name, Group * const parent );
  
  
  /**
   @brief Destructor override from ObjectManager
   */
  virtual ~FaceManager() override;
  
  
  /**
   *@brief Deleted default constructor
   */
  FaceManager() = delete;
  
  
   /**
    *@brief Deleted default constructor
   */
 FaceManager( FaceManager const & ) = delete;
  
  
   /**
    *@brief Deleted copy constructor
   */
 FaceManager( FaceManager && ) = delete;

  ///@}

  /**
   @brief Extend base class resize method resizing  m_nodeList, m_edgeList member containers.
   @param[in] newsize new size of nodeList and edgeList containers
   */
  virtual void resize( localIndex const newsize ) override;

  
 /**
   @brief Build faces in filling face-to-node and face-to-element mappings.
   @param[in] nodeManager mesh node manager
   @param[in] elemManager element manager
 */
  void BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elemManager );

  
  /**
     @brief Compute faces center, area and normal.
     @param[in] nodeManager manager of mesh nodes defining the faces
   */
  void computeGeometry( NodeManager const * const nodeManager );

  
  /**
     @brief Return the number of nodes of the faces with the greatest number of nodes.
     @return the maximum number of nodes a face have
   */
  localIndex getMaxFaceNodes() const;

  
/**
   @brief Sort all faces nodes counterclockwisely in m_nodeList.
   @param[in] nodeManager mesh manager allowing access to face nodes coordinates
   @param[in] elemManager mesh manager allowing access to face elements
 */
  void SortAllFaceNodes( NodeManager const * const nodeManager,
                         ElementRegionManager const * const elemManager );

  
  /**
     @brief Reorder face nodes to be labeled counter-clockwise resulting in outgoing normal.
     @param[in] X array view of mesh nodes coordinates
     @param[in] elemCenter coordinate of the element center
     @param[inout] faceNodes reordered local label list of nodes
     @param[in] numFaceNodes number of nodes for the face
   */
  void SortFaceNodes( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                      R1Tensor const & elemCenter,
                      localIndex * const faceNodes,
                      localIndex const numFaceNodes );

  
  /**
     @brief Flag face and nodes'face with at least one element on the boundary.
     @param[in] nodeManager manager of mesh nodes
   */
  void SetDomainBoundaryObjects( NodeManager * const nodeManager );

  
  /**
     @brief Flag points on boundary or external to the Domain.
   */
  void SetIsExternal();

  
  /**
   * @name Packing methods
   */
  ///@{

  
  /**
   * @brief Create an array listing all excluded local face indices values.
   * @param [inout] exclusionList Sorted array with excluded local faces indices
   */
 virtual void ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  
 /**
  *@brief Calculate the size that a list would have if it were packed, but without actually packing it.
  *@param [in] packList the list of node indices that we wish to get the size of after packing
  *@return a localIndex value representing the size of packList if it were packed
  *@note This function does not perform any packing, it just evaluates and returns the possible packed size.
  */
  virtual localIndex PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Pack an array of node indices into a buffer.
   * @param [inout] buffer buffer to pack the node index data into
   * @param [in] packList the indices of nodes that should be packed
   * @return a localIndex value representing the size of the packed data
   */
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Unpack a buffer to an array of node indices.
   * @param [in] buffer buffer with the packed data
   * @param [inout] packList an array of localIndex values that we wish to unpack to
   * @param [in] overwriteUpMaps boolean: true to overwrite the previous Up maps
   * @param [in] overwriteDownMaps boolean: true to overwrite the previous Down maps
   * @return a localIndex value representing the size of the unpacked list
   */
  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  /**
   * @brief Call FixUpDownMaps for nodes-to-edges and nodes-to-faces maps.
   * @param [in] clearIfUnmapped boolean: true to remove if it is not mapped
   */
  void FixUpDownMaps( bool const clearIfUnmapped );
  ///@}

  /**
   * @brief Compress FaceManager face-to-node and face-to-edge containers so that the values of each array are contiguous with no extra capacity in between.
   * @note The method used here on each arrays (compress) does not free any memory.
   */
  void compressRelationMaps();

  
  /**
   * @brief Enforce child faces and parent faces to have opposite normals.
   * @param[in] targetIndices set of face indices for which the enforcement has to be done
   */
  virtual void enforceStateFieldConsistencyPostTopologyChange( std::set< localIndex > const & targetIndices ) override;
  
  
  /**
   * @brief Clean up the mappings from faces to elment index, region, subregion on a new (updated) list of faces, in order to keep only relevant mappings.
   * @param [in] receivedFaces the new list of target node indices
   * @param [in] elemRegionManager Element Region Manager
   */
  void depopulateUpMaps( std::set< localIndex > const & receivedFaces,
                         ElementRegionManager const & elemRegionManager );
  

  //void SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject );

  /**
   * @brief Extract a face-to-nodes map with global indexed for boundary faces.
   * @param[in] nodeManager mesh nodeManager
   * @param[out] faceToNodes face-to-node map
   */
  virtual void ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const nodeManager,
                                                   std::vector< std::vector< globalIndex > > & faceToNodes ) override;
  

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{
  /**
   *  @struct Containing added view access key to be bound with class data member
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto nodeListString              = "nodeList";
    static constexpr auto edgeListString              = "edgeList";
    static constexpr auto elementRegionListString     = "elemRegionList";
    static constexpr auto elementSubRegionListString  = "elemSubRegionList";
    static constexpr auto elementListString           = "elemList";
    constexpr static auto faceAreaString = "faceArea";
    constexpr static auto faceCenterString = "faceCenter";
    constexpr static auto faceNormalString = "faceNormal";

    dataRepository::ViewKey nodeList              = { nodeListString };
    dataRepository::ViewKey edgeList              = { edgeListString };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString };
    dataRepository::ViewKey elementList           = { elementListString };
  } viewKeys;

 /**
  * @struct Containing added group access key to be bound with class in group hierarchy
  */
  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;
  ///@}

 /**
    @name Accessors for FaceManager fixed data
  */

  ///@{
   /**
    * @brief Get the constant upper limit for numer of faces per Node.
    * @return constant expression of the maximal number of faces per node
   */
 constexpr int maxFacesPerNode() const { return 100; }

  
  /**
   * @brief Get an accessor to face area class member.
   * @return reference to the m_faceArea member
   */
 array1d< real64 > & faceArea()       { return m_faceArea; }
  
  
   /**
    * @brief Get an immutable accessor to face area class member.
    * @return reference to the m_faceArea member
   */
array1d< real64 > const & faceArea() const { return m_faceArea; }
  

  /**
   * @brief Get an accessor to the face center.
   * @return const reference to m_faceCenter class member
   */
  array1d< R1Tensor > & faceCenter()       { return m_faceCenter; }
  
  
   /**
    * @brief Get an immutable accessor to the face center.
    * @return const reference to m_faceCenter class member
   */
  array1d< R1Tensor > const & faceCenter() const { return m_faceCenter; }

  
  /**
   * @brief Get an accessor to the face normals.
   * @return const reference to m_faceNormal class member
   */
  array1d< R1Tensor > & faceNormal()       { return m_faceNormal; }
  
  
   /**
    * @brief Get an immutable accessor to the face normals.
    * @return const reference to m_faceNormal class member
   */
  array1d< R1Tensor > const & faceNormal() const { return m_faceNormal; }
  

  /**
   * @brief Get an accessor to face keyed node map.
   * @return non-const reference to m_nodeList class member
   */
  NodeMapType & nodeList()                    { return m_nodeList; }
  
  
   /**
    * @brief Get an immutable accessor to face keyed node map.
    * @return non-const reference to m_nodeList class member
   */
  NodeMapType const & nodeList() const { return m_nodeList; }
  

  /**
   * @brief  Get an accessor to face keyed edge map.
   * @return non-const reference to m_edgeList class member
   */
  EdgeMapType & edgeList()       { return m_edgeList; }
  
  
   /**
    * @brief  Get an immutable accessor to face keyed edge map.
    * @return non-const reference to m_edgeList class member
   */
  EdgeMapType const & edgeList() const { return m_edgeList; }
  

  /**
   * @brief Get an accessor to the faces-to-elements relation.
   * @return const reference to nodes-to-elements relation
   */
  array2d< localIndex > & elementRegionList()       { return m_toElements.m_toElementRegion; }
  
  
   /**
    * @brief Get an immutable accessor to the faces-to-elements relation.
    * @return const reference to nodes-to-elements relation
   */
  array2d< localIndex > const & elementRegionList() const { return m_toElements.m_toElementRegion; }
  
  
  /**
   * @brief Get an accessor to the faces-to-elements relation.
   * @return const reference to nodes-to-elements relation
   */

  array2d< localIndex > & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }
  
  
   /**
    * @brief Get an accessor to the faces-to-elements relation.
    * @return const reference to nodes-to-elements relation
   */
 array2d< localIndex > const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }
  

  /**
   * @brief Get an immutable accessor to the faces-to-elements list.
   * @return const reference to 2d face-to-element label array
   */
  array2d< localIndex > & elementList()       { return m_toElements.m_toElementIndex; }
  
  
   /**
    *@brief Get an immutable accessor to the faces-to-elements list.
    *@return const reference to 2d face-to-element label array
   */
 array2d< localIndex > const & elementList() const { return m_toElements.m_toElementIndex; }
  

  /**
   * @brief Get an accessor to the faces-to-elements relation.
   * @return const reference to nodes-to-elements relation
   */
  ElemMapType & toElementRelation()       { return m_toElements; }
  
  
   /**
    * @brief Get an immutable accessor to the faces-to-elements relation.
    * @return const reference to nodes-to-elements relation
   */
  ElemMapType const & toElementRelation() const { return m_toElements; }
  
  
  ///}@

private:
  
  
  /**
   * @brief Pack the upward and downward pointing maps.
   * @tparam DOPACK template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                packed and the function returns the size of the packing that would have occured if set to TRUE.
   * @param[inout] buffer the buffer to pack data into
   * @param[in] packList the indices of faces that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;
  

  /// face keyed map containing face-to-node relation
  NodeMapType m_nodeList;
  
  /// face keyed map containing face-to-edge relation
  EdgeMapType m_edgeList;
  
  /// face keyed map containing face-to-element relation
  ElemMapType m_toElements;

  /// map of global  to local  indices for nodes
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;
  
  /// map of global  to local  indices for edges
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  /// list of faces area
  array1d< real64 > m_faceArea;
  
  /// list of faces center
  array1d< R1Tensor > m_faceCenter;
  
  /// list of faces normal
  array1d< R1Tensor > m_faceNormal;

  /// constant expression of the maximum number of nodes per faces
  constexpr static int MAX_FACE_NODES = 9;

};

}
#endif /* FACEMANAGERT_H_ */
