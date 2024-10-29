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
 * @file FieldSpecificationManager.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONMANAGER_HPP_
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONMANAGER_HPP_

#include "FieldSpecificationBase.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{
namespace dataRepository
{
namespace keys
{
/**
 * @brief The key for BoundaryConditionManager
 * @return the key
 */
string const boundaryConditionManager( "BoundaryConditionManager" );
}
}

/**
 * @class FieldSpecificationManager
 * This class contains the field objects and provides an interface for administering
 * to specify them The class is a singleton.
 */
class FieldSpecificationManager : public dataRepository::Group
{
public:

  /**
   * @brief private constructor for the singleton BoundaryConditionManager.
   * @param name The name of the BoundaryConditionManager in the data repository.
   * @param parent The parent of BoundaryConditionManager in the data repository.
   */
  FieldSpecificationManager( string const & name, dataRepository::Group * const parent );

  virtual ~FieldSpecificationManager() override;

  /**
   * @brief @return A pointer to the FieldSpecificationManager
   */
  static FieldSpecificationManager & getInstance();

  /**
   * @brief Create a new FieldSpecificationBase object as a child of this group.
   * @param childKey the catalog key of the new FieldSpecificationBase derived type to create
   * @param childName the name of the new FieldSpecificationBase object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /**
   * @brief Function to apply a value directly to a field variable.
   * @tparam POLICY the policy for kernels launched inside this function.
   * @param time The time at which the value will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param mesh The MeshLevel object.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   *
   * Interface function that applies a boundary condition directly to a field variable by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the FieldSpecificationBase object,
   * and calls FieldSpecificationBase::applyFieldValue().
   *
   */
  template< typename POLICY=parallelHostPolicy >
  void applyFieldValue( real64 const time,
                        MeshLevel & mesh,
                        string const & fieldName ) const
  {
    GEOS_MARK_FUNCTION;

    applyFieldValue< POLICY >( time, mesh, fieldName,
                               [&]( FieldSpecificationBase const &,
                                    SortedArrayView< localIndex const > const & ){} );
  }

  /**
   * @brief Function to apply a value directly to a field variable and applies a lambda
   *        for any post operations that are needed.
   * @tparam POLICY The execution policy for kernels launched in this function.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param mesh The MeshLevel object.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   * @param lambda A lambda function that defines any operations that should be performed
   *               after application of value to the field.
   *
   * Interface function that applies a value directly to a field by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the FieldSpecificationBase object,
   * and calls FieldSpecificationBase::applyFieldValue(), and calls the lambda function
   * to apply any operations required for completing the application of the value to the field in addition to
   * setting the target field.
   */
  template< typename POLICY=parallelHostPolicy, typename LAMBDA=void >
  void applyFieldValue( real64 const time,
                        MeshLevel & mesh,
                        string const & fieldName,
                        LAMBDA && lambda ) const;

  /**
   * @brief Function to apply a value directly to a field variable and applies a lambda
   *        for any post operations that are needed.
   * @tparam POLICY The execution policy for kernels launched in this function.
   * @tparam PRELAMBDA The type of the lambda function to be called before the applyFieldValue
   * @tparam POSTLAMBDA The type of the lambda function to be called after the applyFieldValue
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param mesh The MeshLevel objectt.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   * @param preLambda A lambda function that defines any operations that should be performed
   *                  before application of value to the field.
   * @param postLambda A lambda function that defines any operations that should be performed
   *                   after application of value to the field.
   *
   * Interface function that applies a value directly to a field by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Calls the
   * preLambda function to apply any operations required before the application of the value to the
   * field in addition to setting the target field. Searches for the string specified in fieldPath
   * as a substring in the objectPath specified in the input file, checks if fieldName is equal to
   * the fieldName specified in the input file, and check if the time parameter falls within the
   * beginTime and endTime of the FieldSpecificationBase object, and calls
   * FieldSpecificationBase::applyFieldValue(), and calls the postLambda function to apply any
   * operations required for completing the application of the value to the field in addition to
   * setting the target field.
   */
  template< typename POLICY=parallelHostPolicy, typename PRELAMBDA=void, typename POSTLAMBDA=void >
  void applyFieldValue( real64 const time,
                        MeshLevel & mesh,
                        string const & fieldName,
                        PRELAMBDA && preLambda,
                        POSTLAMBDA && postLambda ) const;


  /**
   * @brief function to apply initial conditions
   * @param mesh the MeshLevel object
   */
  void applyInitialConditions( MeshLevel & mesh ) const;

  /**
   * @brief function to validate the application of boundary conditions
   * @param mesh the MeshLevel object
   */
  void validateBoundaryConditions( MeshLevel & mesh ) const;


  /**
   * @brief This function is the main driver for the field applications
   * @tparam OBJECT_TYPE the type of object that the application targets
   * @tparam BCTYPE the type of boundary condition object that is being applied
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param mesh The MeshLevel object.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the specification.
   * @param lambda A lambda function that defines the application of the field.
   *
   * This function loops through all available fields, checks to see if they
   * should be applied, and applies them. More specifically, this function simply checks
   * values of fieldPath,fieldName, against each FieldSpecificationBase object contained in the
   * FieldSpecificationManager and decides on whether or not to call the user defined lambda.
   */
  template< typename OBJECT_TYPE=dataRepository::Group,
            typename BCTYPE = FieldSpecificationBase,
            typename LAMBDA >
  void apply( real64 const time,
              MeshLevel & mesh,
              string const & fieldName,
              LAMBDA && lambda ) const
  {
    GEOS_MARK_FUNCTION;

    string const meshBodyName = mesh.getParent().getParent().getName();
    string const meshLevelName = mesh.getName();

    // loop over all FieldSpecificationBase objects
    this->forSubGroups< BCTYPE >( [&] ( BCTYPE const & fs )
    {
      integer const isInitialCondition = fs.initialCondition();
      if( ( isInitialCondition && fieldName=="") || // this only use case for this line is in the unit test for field specification
          ( !isInitialCondition && time >= fs.getStartTime() && time < fs.getEndTime() && fieldName == fs.getFieldName() ) )
      {
        fs.template apply< OBJECT_TYPE, BCTYPE, LAMBDA >( mesh, std::forward< LAMBDA >( lambda ) );
      }
    } );
  }

private:
  static FieldSpecificationManager * m_instance;

};

template< typename POLICY, typename LAMBDA >
void
FieldSpecificationManager::
  applyFieldValue( real64 const time,
                   MeshLevel & mesh,
                   string const & fieldName,
                   LAMBDA && lambda ) const
{
  GEOS_MARK_FUNCTION;

  apply( time, mesh, fieldName,
         [&]( FieldSpecificationBase const & fs,
              string const &,
              SortedArrayView< localIndex const > const & targetSet,
              Group & targetGroup,
              string const & targetField )
  {
    fs.applyFieldValue< FieldSpecificationEqual, POLICY >( targetSet, time, targetGroup, targetField );
    lambda( fs, targetSet );
  } );
}

template< typename POLICY, typename PRELAMBDA, typename POSTLAMBDA >
void
FieldSpecificationManager::
  applyFieldValue( real64 const time,
                   MeshLevel & mesh,
                   string const & fieldName,
                   PRELAMBDA && preLambda,
                   POSTLAMBDA && postLambda ) const
{
  GEOS_MARK_FUNCTION;

  apply( time, mesh, fieldName,
         [&]( FieldSpecificationBase const & fs,
              string const &,
              SortedArrayView< localIndex const > const & targetSet,
              Group & targetGroup,
              string const & targetField )
  {
    preLambda( fs, targetSet );
    fs.applyFieldValue< FieldSpecificationEqual, POLICY >( targetSet, time, targetGroup, targetField );
    postLambda( fs, targetSet );
  } );
}

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONMANAGER_HPP_ */
