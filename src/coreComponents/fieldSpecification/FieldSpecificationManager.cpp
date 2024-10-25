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

#include "FieldSpecificationManager.hpp"

#include "common/format/StringUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geos
{

FieldSpecificationManager * FieldSpecificationManager::m_instance = nullptr;

using namespace dataRepository;
using namespace constitutive;
FieldSpecificationManager::FieldSpecificationManager( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  GEOS_ERROR_IF( m_instance != nullptr, "Only one FieldSpecificationManager can exist at a time." );
  m_instance = this;

}

FieldSpecificationManager::~FieldSpecificationManager()
{
  GEOS_ERROR_IF( m_instance != this, "m_instance != this should not be possible." );
  m_instance = nullptr;
}


FieldSpecificationManager & FieldSpecificationManager::getInstance()
{
  GEOS_ERROR_IF( m_instance == nullptr,
                 "FieldSpecificationManager has not been constructed, or is already been destructed." );
  return *m_instance;
}

Group * FieldSpecificationManager::createChild( string const & childKey, string const & childName )
{
  std::unique_ptr< FieldSpecificationBase > bc = FieldSpecificationBase::CatalogInterface::factory( childKey, getDataContext(),
                                                                                                    childName, this );
  return &this->registerGroup( childName, std::move( bc ) );
}


void FieldSpecificationManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BoundaryConditionBase here
  for( auto & catalogIter: FieldSpecificationBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

void FieldSpecificationManager::validateBoundaryConditions( MeshLevel & mesh ) const
{
  // loop over all the FieldSpecification of the XML file
  this->forSubGroups< FieldSpecificationBase >( [&] ( FieldSpecificationBase const & fs )
  {
    localIndex isFieldNameFound = 0;
    // map from set name to a flag (1 if targetSet empty, 0 otherwise)
    map< string, localIndex > isTargetSetEmpty;
    // map from set name to a flag (1 if targetSet has been created, 0 otherwise)
    map< string, localIndex > isTargetSetCreated;

    // Step 1: collect all the set names in a map (this is made necessary by the "apply" loop pattern

    array1d< string > const & setNames = fs.getSetNames();
    for( localIndex i = 0; i < setNames.size(); ++i )
    {
      isTargetSetEmpty[setNames[i]] = 1;
      isTargetSetCreated[setNames[i]] = 0;
    }

    // We have to make sure that the meshLevel is in the target of the boundary conditions
    // This is important for multi-level simulations, such as high-order wave propagation
    MeshObjectPath const & objectPath = fs.getMeshObjectPaths();
    if( !objectPath.containsMeshLevel( mesh ) )
    {
      return;
    }

    // Step 2: apply the boundary condition

    fs.apply< dataRepository::Group >( mesh,
                                       [&]( FieldSpecificationBase const &,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const & targetSet,
                                            Group & targetGroup,
                                            string const fieldName )
    {
      dataRepository::InputFlags const flag = fs.getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).getInputFlag();

      // 2.a) If we enter this loop, we know that the set has been created
      //      Fracture/fault sets are created later and the "apply" call silently ignores them
      isTargetSetCreated.at( setName ) = 1;

      // 2.b) If the fieldName is registered on this target, we record it
      //      Unfortunately, we need two exceptions:
      //       - FieldSpecification that do not target a field, like Aquifer, Traction, Equilibrium, etc. For these, the check is not
      // necessary (the user cannot mess up)
      //       - Face boundary conditions that target cell-based quantities, like the face BC of the flow solvers
      if( targetGroup.hasWrapper( fieldName ) ||
          flag == InputFlags::FALSE || // no need to check input if the input flag is false (Aquifer, Traction, Equilibrium do not target a
                                       // field)
          targetGroup.getName() == MeshLevel::groupStructKeys::faceManagerString() ) // the field names of the face BCs are not always
                                                                                     // registered on
                                                                                     // the faceManager...
      {
        isFieldNameFound = 1;
      }

      // 2.c) If the target set is not empty, we record it
      //      Fracture/fault sets are sometimes at initialization and the "apply" call silently ignores them
      if( targetSet.size() > 0 )
      {
        isTargetSetEmpty.at( setName ) = 0;
      }
    } );

    // Step 3: MPI synchronization

    isFieldNameFound = MpiWrapper::max( isFieldNameFound );

    for( std::pair< string const, localIndex > & mapEntry : isTargetSetEmpty )
    {
      mapEntry.second = MpiWrapper::min( mapEntry.second );
    }
    bool areAllSetsEmpty = true;
    for( std::pair< string const, localIndex > & mapEntry : isTargetSetEmpty )
    {
      if( mapEntry.second == 0 ) // target set is not empty
      {
        areAllSetsEmpty = false;
        break;
      }
    }

    for( std::pair< string const, localIndex > & mapEntry : isTargetSetCreated )
    {
      mapEntry.second = MpiWrapper::max( mapEntry.second );
    }
    bool areAllSetsMissing = true;
    for( std::pair< string const, localIndex > & mapEntry : isTargetSetCreated )
    {
      if( mapEntry.second == 1 ) // target set has been created
      {
        areAllSetsMissing = false;
      }
    }

    // Step 4: issue an error or a warning if the field was not found

    // if all sets are missing, we stop the simulation.
    if( areAllSetsMissing )
    {
      // loop again over the map to collect the set names
      array1d< string > missingSetNames;
      for( auto const & mapEntry : isTargetSetCreated )
      {
        missingSetNames.emplace_back( mapEntry.first );
      }
      GEOS_THROW( GEOS_FMT( "\n{}: there is/are no set(s) named `{}` under the {} `{}`.\n",
                            fs.getWrapperDataContext( FieldSpecificationBase::viewKeyStruct::objectPathString() ),
                            fmt::join( missingSetNames, ", " ),
                            FieldSpecificationBase::viewKeyStruct::objectPathString(), fs.getObjectPath() ),
                  InputError );
    }

    // if a target set is empty, we issue a warning
    // ideally we would just stop the simulation, but the SurfaceGenerator relies on this behavior
    for( auto const & mapEntry : isTargetSetEmpty )
    {
      GEOS_LOG_RANK_0_IF( ( mapEntry.second == 1 ), // target set is empty
                          GEOS_FMT( "\nWarning!\n{}: this FieldSpecification targets (an) empty set(s)"
                                    "\nIf the simulation does not involve the SurfaceGenerator, check the content of the set `{}` in `{}`. \n",
                                    fs.getDataContext(), mapEntry.first, fs.getObjectPath() ) );
    }

    if( isFieldNameFound == 0 )
    {
      char const fieldNameNotFoundMessage[] =
        "\n{}: there is no {} named `{}` under the {} `{}`.\n";
      string const errorMsg =
        GEOS_FMT( fieldNameNotFoundMessage,
                  fs.getWrapperDataContext( FieldSpecificationBase::viewKeyStruct::fieldNameString() ),
                  FieldSpecificationBase::viewKeyStruct::fieldNameString(),
                  fs.getFieldName(), FieldSpecificationBase::viewKeyStruct::objectPathString(), fs.getObjectPath() );
      if( areAllSetsEmpty )
      {
        GEOS_LOG_RANK_0( errorMsg );
      }
      else
      {
        GEOS_THROW( errorMsg, InputError );
      }
    }
  } );
}

void FieldSpecificationManager::applyInitialConditions( MeshLevel & mesh ) const
{
  this->forSubGroups< FieldSpecificationBase >( [&] ( FieldSpecificationBase const & fs )
  {
    if( fs.initialCondition() )
    {
      fs.apply< dataRepository::Group >( mesh,
                                         [&]( FieldSpecificationBase const & bc,
                                              string const &,
                                              SortedArrayView< localIndex const > const & targetSet,
                                              Group & targetGroup,
                                              string const fieldName )
      {
        bc.applyFieldValue< FieldSpecificationEqual >( targetSet, 0.0, targetGroup, fieldName );
      } );
    }
  } );
}

} /* namespace geos */
