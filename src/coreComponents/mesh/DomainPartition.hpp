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
 * @file DomainPartition.hpp
 */

#ifndef GEOS_MESH_DOMAINPARTITION_HPP_
#define GEOS_MESH_DOMAINPARTITION_HPP_

#include "common/MpiWrapper.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

namespace geos
{

class SiloFile;
// namespace dataRepository
// {
// namespace keys
// {
// /// @return PartitionManager string key
// string const partitionManager( "partitionManager" );
// }
// }

class ObjectManagerBase;
class PartitionBase;

/**
 * @brief Partition of the decomposed physical domain. It also manages the connexion information to its neighbors.
 *
 * Two types of neighors are being managed. One MPI cartesian communicator that DomainPartition shall build.
 * Or through `metis`, but it is not the responsibility of DomainPartition to build the decomposition in that case.
 */
class DomainPartition : public dataRepository::Group
{
public:
  /**
   * @brief Constructor.
   * @param[in] name Name of this object manager
   * @param[in] parent Parent Group
   */
  DomainPartition( string const & name,
                   Group * const parent );

  /**
   * @brief Destructor.
   */
  ~DomainPartition() override;

  /**
   * @name Default, copy and assignment constructors are deleted.
   */
  ///@{
  /// @cond DO_NOT_DOCUMENT
  DomainPartition() = delete;

  DomainPartition( DomainPartition const & ) = delete;

  DomainPartition( DomainPartition && ) = delete;

  DomainPartition & operator=( DomainPartition const & ) = delete;

  DomainPartition & operator=( DomainPartition && ) = delete;
  /// @endcond
  ///@}

  void initializationOrder( string_array & order ) override final;

  /**
   * @name MPI functionality
   */
  ///@{
  /**
   * @brief Constructs the communications between this DomainPartition and its neighbors.
   * @param use_nonblocking If true complete the communications of each phase in the order they are received.
   */
  void setupCommunications( bool use_nonblocking );

  /**
   * @brief Constructs the global information of this DomainPartition, needed to set up ghosting
   */
  void setupBaseLevelMeshGlobalInfo();

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * This recursive function builds the neighbors recursively by increasing
   * the dimension index of the current DomainPartition until all the dimensions (3) are done.
   * The relevant values for initiating this functions are therefore @p ibim = 0
   * and a non-initialized vector @p ncoords of length 3.
   *
   * This functions should have been implemented `private`
   * and an additional functions to initiate the recursion could have been implemented.
   */
  void addNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );
  ///@}


  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeysStruct
   */
  struct groupKeysStruct
  {
    /// @return String key to the Group holding the MeshBodies
    static constexpr char const * meshBodiesString() { return "MeshBodies"; }
    /// @return String key to the Group holding the ConstitutiveManager
    static constexpr char const * constitutiveManagerString() { return "Constitutive"; }

    /// View key to the Group holding the MeshBodies
    dataRepository::GroupKey meshBodies = { meshBodiesString() };
    /// View key to the Group holding the ConstitutiveManager
    dataRepository::GroupKey constitutiveManager = { constitutiveManagerString() };
    /// View key to the Group holding the CommunicationManager
    dataRepository::GroupKey communicationManager = { "communicationManager" };
    /// View key to the Group holding the partitionManager
    dataRepository::GroupKey partitionManager = { "partitionManager" };
  }
  /// groupKey struct for the DomainPartition class
  groupKeys;

  /**
   * @brief Get the constitutive manager, const version.
   * @return Pointer to a const instance of a ConstitutiveManager.
   */
  constitutive::ConstitutiveManager const & getConstitutiveManager() const
  { return this->getGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager ); }

  /**
   * @brief Get the constitutive manager.
   * @return Pointer to an instance of a ConstitutiveManager.
   */
  constitutive::ConstitutiveManager & getConstitutiveManager()
  { return this->getGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager ); }

  /**
   * @brief @return Return a reference to const NumericalMethodsManager from ProblemManager
   */
  NumericalMethodsManager const & getNumericalMethodManager() const
  { return this->getParent().getGroup< NumericalMethodsManager >( "NumericalMethods" ); }

  /**
   * @brief @return Return a reference to NumericalMethodsManager from ProblemManager
   */
  NumericalMethodsManager & getNumericalMethodManager()
  { return this->getParent().getGroup< NumericalMethodsManager >( "NumericalMethods" ); }

  /**
   * @brief Get the mesh bodies, const version.
   * @return Reference to a const instance of a Group that contains MeshBody instances.
   */
  Group const & getMeshBodies() const
  { return this->getGroup( groupKeys.meshBodies ); }

  /**
   * @brief Get the mesh bodies.
   * @return Reference to a instance of a Group that contains MeshBody instances.
   */
  Group & getMeshBodies()
  { return this->getGroup( groupKeys.meshBodies ); }

  /**
   * @brief Check if a MeshBody is present given a name.
   * @tparam KEY_TYPE The type of the key used to look up the MeshBody.
   * @param key The key to the MeshBody.
   * @return True is the MeshBody exists in the domain.
   */
  template< typename KEY_TYPE >
  bool hasMeshBody( KEY_TYPE const & key ) const
  { return getMeshBodies().hasGroup< MeshBody >( key ); }

  /**
   * @brief Get a MeshBody by name, const version.
   * @tparam KEY_TYPE The type of the key used to look up the MeshBody.
   * @param key The key to the MeshBody.
   * @return Reference to a const MeshBody instance matching @p key.
   */
  template< typename KEY_TYPE >
// TODO uncomment these to disallow using a non string or char const * key
//            std::enable_if_t< std::is_same< T, string >::value ||
//                              std::is_same< T, const char * >::value, bool > = false >
  MeshBody const & getMeshBody( KEY_TYPE const & key ) const
  { return getMeshBodies().getGroup< MeshBody >( key ); }

  /**
   * @brief Get a MeshBody by name.
   * @tparam KEY_TYPE The type of the key used to look up the MeshBody.
   * @param key The key to the MeshBody.
   * @return Reference to a const MeshBody instance matching @p key.
   */
  template< typename KEY_TYPE >
//  std::enable_if_t< std::is_same< T, string >::value ||
//                    std::is_same< T, const char * >::value, bool > = false >
  MeshBody & getMeshBody( KEY_TYPE const & key )
  { return getMeshBodies().getGroup< MeshBody >( key ); }


  /**
   * @brief Apply the given functor to all meshBodies.
   * @tparam FUNCTION the type of functor to call
   * @param[in] function  the functor to call
   */
  template< typename FUNCTION >
  void forMeshBodies( FUNCTION && function ) const
  {
    getMeshBodies().forSubGroups< MeshBody >( std::forward< FUNCTION >( function ) );
  }

  /**
   * @copydoc forMeshBodies(FUNCTION &&) const
   */
  template< typename FUNCTION >
  void forMeshBodies( FUNCTION && function )
  {
    getMeshBodies().forSubGroups< MeshBody >( std::forward< FUNCTION >( function ) );
  }

  /**
   * @copydoc forMeshBodies(FUNCTION &&) const
   */
  template< typename FUNCTION >
  void forMeshBodiesIndex( FUNCTION && function ) const
  {
    getMeshBodies().forSubGroupsIndex< MeshBody >( std::forward< FUNCTION >( function ) );
  }

  /**
   * @copydoc forMeshBodies(FUNCTION &&) const
   */
  template< typename FUNCTION >
  void forMeshBodiesIndex( FUNCTION && function )
  {
    getMeshBodies().forSubGroupsIndex< MeshBody >( std::forward< FUNCTION >( function ) );
  }

  /**
   * @brief Get the neighbor communicators. @see DomainPartition#m_neighbors.
   * @return Container of communicators.
   */
  std::vector< NeighborCommunicator > & getNeighbors()
  { return m_neighbors; }

  /**
   * @brief Get the neighbor communicators, const version. @see DomainPartition#m_neighbors.
   * @return Container of communicators.
   */
  std::vector< NeighborCommunicator > const & getNeighbors() const
  { return m_neighbors; };

private:

  /**
   * @brief Contains all the communicators from this DomainPartition to its neighbors.
   */
  std::vector< NeighborCommunicator > m_neighbors; 

};

} /* namespace geos */

#endif /* GEOS_MESH_DOMAINPARTITION_HPP_ */
