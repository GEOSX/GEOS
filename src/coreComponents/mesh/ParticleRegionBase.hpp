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

#ifndef GEOS_MESH_PARTICLEREGIONBASE_HPP
#define GEOS_MESH_PARTICLEREGIONBASE_HPP

#include "ParticleSubRegion.hpp"
#include "mesh/ObjectManagerBase.hpp"

namespace geos
{

/**
 * @class ParticleRegionBase
 * @brief The ParticleRegionBase is the base class to manage the data stored at the particle level.
 *
 * The ParticleRegion base is the base class for the ParticleRegion class. It may be depreciated at
 * some point since no other classes are currently derived from ParticleRegionBase.
 */
class ParticleRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Deleted default constructor.
   */
  ParticleRegionBase() = delete;

  /**
   * @brief Main constructor.
   * @param name the name of the particle region
   * @param parent the pointer to the parent group
   */
  ParticleRegionBase( string const & name, Group * const parent );


  /**
   * @brief Copy constructor.
   * @param init the particle region to be copied
   */
  ParticleRegionBase( const ParticleRegionBase & init );

  /**
   * @brief Default destructor.
   */
  virtual ~ParticleRegionBase() override;

  ///@}


  /**
   * @brief verify that the meshBody name specified exists in the meshBodies group.
   * If there is only one meshBody it returns the name of the only existing mesh body.
   *
   * @param[in] meshBodies the meshBodies group.
   * @param[in] meshBodyBlockName name of the meshbody to be verified
   * @return the name of the only mesh body present if there is only one
   * OR the name specified by @p  meshBodyName if the meshBody is found.
   */
  static string verifyMeshBodyName( Group const & meshBodies,
                                    string const & meshBodyBlockName );

  /**
   * @name Generation of the mesh region
   */
  ///@{

  /**
   * @brief Generate mesh.
   * @param particleBlocks particle blocks where the mesh is generated
   */
  virtual void generateMesh( Group & particleBlocks )
  {
    GEOS_UNUSED_VAR( particleBlocks );
    GEOS_ERROR( "ParticleRegionBase::GenerateMesh() should be overriden if called." );
  }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @copydoc getSubRegions() const
   */
  subGroupMap & getSubRegions()
  {
    return getGroup( viewKeyStruct::particleSubRegions() ).getSubGroups();
  }

  /**
   * @brief Get a collection of the subregions.
   * @return a collection of the subregions
   */
  subGroupMap const & getSubRegions() const
  {
    return getGroup( viewKeyStruct::particleSubRegions() ).getSubGroups();
  }


  /**
   * @brief Get a reference to a subregion.
   * @tparam SUBREGIONTYPE the type that will be used to attempt casting the subregion
   * @tparam KEY_TYPE The type of the key used to lookup the subregion.
   * @param key The key to the subregion.
   * @return A reference to the subregion
   */
  template< typename SUBREGIONTYPE=ParticleSubRegionBase, typename KEY_TYPE=void >
  SUBREGIONTYPE const & getSubRegion( KEY_TYPE const & key ) const
  {
    return this->getGroup( viewKeyStruct::particleSubRegions() ).getGroup< SUBREGIONTYPE >( key );
  }

  /**
   * @copydoc getSubRegion( KEY_TYPE const & key ) const
   */
  template< typename SUBREGIONTYPE=ParticleSubRegionBase, typename KEY_TYPE=void >
  SUBREGIONTYPE & getSubRegion( KEY_TYPE const & key )
  {
    return this->getGroup( viewKeyStruct::particleSubRegions() ).getGroup< SUBREGIONTYPE >( key );
  }

  /**
   * @brief Get the number of subregions in the region.
   * @return the number of subregions  in the region
   */
  localIndex numSubRegions() const
  {
    return this->getGroup( viewKeyStruct::particleSubRegions() ).getSubGroups().size();
  }

  /**
   * @brief Get the number of particles in the region
   *        for specific subregion types provided as template arguments.
   * @tparam SUBREGIONTYPE  the first type that will be used in the attempted casting of the subregion
   * @tparam SUBREGIONTYPES a variadic list of types that will be used in the attempted casting of the subregion
   * @return the number of particles contained in the particle region
   * @note This function requires that the subRegion types specified
   *       in the variadic template argument can be casted to ParticleSubRegionBase
   */
  template< typename SUBREGIONTYPE = ParticleSubRegionBase, typename ... SUBREGIONTYPES >
  localIndex getNumberOfParticles() const
  {
    localIndex numParticle = 0;
    this->forParticleSubRegions< SUBREGIONTYPE, SUBREGIONTYPES... >( [&]( Group const & group ) -> void
    {
      numParticle += group.size();
    } );
    return numParticle;
  }

  /**
   * @copydoc getMaterialList() const
   */
  string_array & getMaterialList() {return m_materialList;}

  /**
   * @brief Get the material list in the particle region.
   * @return the material list
   */
  string_array const & getMaterialList() const {return m_materialList;}

  /**
   * @brief Get the name of the constitutive in the particle region.
   * @tparam CONSTITUTIVE_TYPE the type of the constitutive model
   * @return the string_array of the constitutive names
   */
  template< typename CONSTITUTIVE_TYPE >
  string_array getConstitutiveNames() const;


  ///@}


  /**
   * @name Functor-based loops over particle subregions
   */
  ///@{


  /**
   * @brief Apply a lambda to all particle subregions.
   * @param lambda the functor to be applied
   */
  template< typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda ) const
  {
    forParticleSubRegions< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @copydoc forParticleSubRegions( LAMBDA && lambda ) const
   */
  template< typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda )
  {
    forParticleSubRegions< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief Apply LAMBDA to the subregions with the specific subregion types
   *        listed in the template.
   * @param lambda the functor to be applied
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda ) const
  {
    this->getGroup( viewKeyStruct::particleSubRegions() ).forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @copydoc forParticleSubRegions( LAMBDA && lambda ) const
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda )
  {
    this->getGroup( viewKeyStruct::particleSubRegions() ).forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief Apply LAMBDA to the subregions, loop using subregion indices.
   * @tparam LAMBDA type of functor to call
   * @param lambda the functor to be applied to all subregions
   */
  template< typename LAMBDA >
  void forParticleSubRegionsIndex( LAMBDA && lambda ) const
  {
    forParticleSubRegionsIndex< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @copydoc forParticleSubRegionsIndex( LAMBDA && lambda ) const
   */
  template< typename LAMBDA >
  void forParticleSubRegionsIndex( LAMBDA && lambda )
  {
    forParticleSubRegionsIndex< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief Apply LAMBDA to the subregions with the specific subregion types
   *        listed in the template, loop using subregion indices.
   * @tparam SUBREGIONTYPE the first subregion type
   * @tparam SUBREGIONTYPES a variadic list of subregion types
   * @tparam LAMBDA type of functor to call
   * @param lambda the functor to be applied
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forParticleSubRegionsIndex( LAMBDA && lambda ) const
  {
    for( localIndex esr=0; esr<this->numSubRegions(); ++esr )
    {
      ParticleSubRegionBase const & subRegion = this->getSubRegion( esr );
      applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      } );
    }
  }

  /**
   * @copydoc forParticleSubRegionsIndex( LAMBDA && lambda ) const
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forParticleSubRegionsIndex( LAMBDA && lambda )
  {
    for( localIndex esr=0; esr<this->numSubRegions(); ++esr )
    {
      ParticleSubRegionBase & subRegion = this->getSubRegion( esr );
      applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      } );
    }
  }

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the material list
    static constexpr char const * materialListString() { return "materialList"; }
    /// @return String key for the material list
    static constexpr char const * meshBodyString() { return "meshBody"; }
    /// @return String key for the particle subregions
    static constexpr char const * particleSubRegions() { return "particleSubRegions"; }
  };

private:

  ParticleRegionBase & operator=( const ParticleRegionBase & rhs );

  /// List of materials for the particle region
  string_array m_materialList;

  /// Name of the mesh body that contains this region
  string m_meshBody;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


/**
 * @brief Get the names of all constitutive models of a specific type
 *
 * @tparam CONSTITUTIVE_TYPE type of constitutive model
 * return string array with the names of the constitutive models
 */
template< typename CONSTITUTIVE_TYPE >
string_array ParticleRegionBase::getConstitutiveNames() const
{
  string_array rval;
  for( string const & matName : m_materialList )
  {
    if( this->getSubRegion( 0 ).getConstitutiveModels().hasGroup< CONSTITUTIVE_TYPE >( matName ) )
    {
      rval.emplace_back( matName );
    }
  }
  return rval;
}

}



#endif /* GEOS_MESH_PARTICLEREGIONBASE_HPP */
