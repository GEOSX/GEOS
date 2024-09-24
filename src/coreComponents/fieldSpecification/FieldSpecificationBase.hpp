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

/**
 * @file FieldSpecificationBase.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP

#include "common/DataTypes.hpp"
#include "common/TypeDispatch.hpp"
#include "codingUtilities/traits.hpp"
#include "codingUtilities/Utilities.hpp"
#include "dataRepository/Group.hpp"
#include "functions/FunctionBase.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/MeshObjectPath.hpp"
#include "functions/FunctionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{
class Function;


/**
 * @class FieldSpecificationBase
 * A class to hold values for and administer a single boundary condition
 */
class FieldSpecificationBase : public dataRepository::Group
{
public:

  /**
   * @defgroup alias and functions to defined statically initialized catalog
   * @{
   */

  /**
   * alias to define the catalog type for this base type
   */
  using CatalogInterface = dataRepository::CatalogInterface< FieldSpecificationBase,
                                                             string const &,
                                                             dataRepository::Group * const >;

  /**
   * @brief static function to return static catalog.
   * @return the static catalog to create derived types through the static factory methods.
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "FieldSpecification"; }

  /**
   * @brief return the catalog name
   * @return the catalog name
   */
  virtual const string getCatalogName() const
  {
    return FieldSpecificationBase::catalogName();
  }

  /**
   * @}
   */


  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationBase in the data repository
   * @param parent the parent group of this group.
   */
  FieldSpecificationBase( string const & name, dataRepository::Group * parent );

  /**
   * destructor
   */
  virtual ~FieldSpecificationBase() override;


  /// Deleted copy constructor
  FieldSpecificationBase( FieldSpecificationBase const & ) = delete;

  /// Defaulted move constructor
  FieldSpecificationBase( FieldSpecificationBase && ) = default;

  /// deleted copy assignment
  FieldSpecificationBase & operator=( FieldSpecificationBase const & ) = delete;

  /// deleted move assignement
  FieldSpecificationBase & operator=( FieldSpecificationBase && ) = delete;

  /**
   * @brief Apply this field specification to the discretization
   *
   * @tparam OBJECT_TYPE The type of discretization/mesh object that the
   *   specification is being applied to.
   * @tparam BC_TYPE The type of BC being applied
   * @tparam LAMBDA
   * @param mesh The MeshLevel that the specification is applied to
   * @param lambda The being executed
   */
  template< typename OBJECT_TYPE,
            typename BC_TYPE = FieldSpecificationBase,
            typename LAMBDA >
  void apply( MeshLevel & mesh,
              LAMBDA && lambda ) const
  {
    MeshObjectPath const & meshObjectPaths = this->getMeshObjectPaths();
    meshObjectPaths.forObjectsInPath< OBJECT_TYPE >( mesh,
                                                     [&] ( OBJECT_TYPE & object )
    {
      {
        dataRepository::Group const & setGroup = object.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );
        string_array setNames = this->getSetNames();
        for( auto & setName : setNames )
        {
          if( setGroup.hasWrapper( setName ) )
          {
            SortedArrayView< localIndex const > const & targetSet = setGroup.getReference< SortedArray< localIndex > >( setName );
            lambda( dynamic_cast< BC_TYPE const & >(*this), setName, targetSet, object, getFieldName() );
          }
        }
      }
    } );
  }

  /**
   * @tparam FIELD_OP type that contains static functions to apply the value to the field
   * @param[in] field the field to apply the value to.
   * @param[in] targetSet the set of indices which the value will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the value.
   * @param[in] dataGroup the Group that contains the field to apply the value to.
   *
   * This function applies the value to a field variable.
   */
  template< typename FIELD_OP, typename POLICY, typename T, int N, int USD >
  void applyFieldValueKernel( ArrayView< T, N, USD > const & field,
                              SortedArrayView< localIndex const > const & targetSet,
                              real64 const time,
                              Group & dataGroup ) const;

  /**
   * @tparam FIELD_OP type that contains static functions to apply the value to the field
   * @param[in] targetSet the set of indices which the value will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the value.
   * @param[in] dataGroup the Group that contains the field to apply the value to.
   * @param[in] fieldname the name of the field to apply the value to.
   *
   * This function applies the value to a field variable. This function is typically
   * called from within the lambda to a call to FieldSpecificationManager::applyFieldValue().
   */
  template< typename FIELD_OP, typename POLICY=parallelHostPolicy >
  void applyFieldValue( SortedArrayView< localIndex const > const & targetSet,
                        real64 const time,
                        dataRepository::Group & dataGroup,
                        string const & fieldname ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over the target set.
   * @tparam T Data type of the field.
   * @tparam NDIM Number of dimensions in the field array.
   * @tparam USD Unit stride dimension of the field array.
   * @param targetSet The set of indices which the boundary condition will be applied.
   * @param time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param dataGroup The Group that contains the field to apply the boundary condition to.
   * @param dofMap The map from the local index of the primary field to the global degree of freedom number.
   * @param dofRankOffset Offset of dof indices on current rank.
   * @param matrix Local part of the system matrix.
   * @param rhs Local part of the system rhs vector.
   * @param fieldView Array view of the field data.
   *
   * @note This function is rarely used directly. More often it is called by other ApplyBoundaryCondition functions.
   */
  template< typename FIELD_OP, typename POLICY, typename T, int NDIM, int USD >
  void applyBoundaryConditionToSystemKernel( SortedArrayView< localIndex const > const & targetSet,
                                             real64 const time,
                                             dataRepository::Group const & dataGroup,
                                             arrayView1d< globalIndex const > const & dofMap,
                                             globalIndex const dofRankOffset,
                                             CRSMatrixView< real64, globalIndex const > const & matrix,
                                             arrayView1d< real64 > const & rhs,
                                             ArrayView< T const, NDIM, USD > const & fieldView ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *            application of the boundary condition.
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] fieldName The name of the field to apply the boundary condition to.
   * @param[in] dofMapName The name of the map from the local index of the primary field to the
   *                       global degree of freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[in,out] matrix Local part of the system matrix.
   * @param[in,out] rhs Local part of the system rhs vector.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename POLICY >
  void applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                       real64 const time,
                                       dataRepository::Group const & dataGroup,
                                       string const & fieldName,
                                       string const & dofMapName,
                                       globalIndex const dofRankOffset,
                                       CRSMatrixView< real64, globalIndex const > const & matrix,
                                       arrayView1d< real64 > const & rhs ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[inout] matrix Local part of the system matrix.
   * @param[inout] rhs Local part of the system rhs vector.
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                   provided by the FIELD_OP templated type.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename POLICY, typename LAMBDA >
  void
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dt time step size which is applied as a factor to bc values
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[inout] matrix Local part of the system matrix.
   * @param[inout] rhs Local part of the system rhs vector.
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                   provided by the FIELD_OP templated type.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename POLICY, typename LAMBDA >
  void
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const;

  /**
   * @brief Compute the contributions that will be added/enforced to the right-hand side, and collect the corresponding dof numbers
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dt time step size which is applied as a factor to bc values
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[inout] matrix Local part of the system matrix.
   * @param[inout] dof array storing the degrees of freedom of the rhsContribution, to know where in the rhs they will be added/enforced
   * @param[inout] rhsContribution array storing the values that will be added/enforced to the right-hand side
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                   provided by the FIELD_OP templated type.
   *
   * Note that this function only computes the rhs contributions, but does not apply them to the right-hand side.
   * The application of these rhs contributions is done in applyBoundaryConditionToSystem.
   *
   * Why did we have to extract the computation of the rhs contributions from applyBoundaryConditionToSystem?
   * Because applyBoundaryConditionToSystem is not very well suited to apply the rhsContributions to the equation layout used in the
   * compositional solvers.
   * Therefore, the compositional solvers do not call applyBoundaryConditionToSystem, but instead call computeRhsContribution directly, and
   * apply these rhs contributions "manually" according to the equation layout used in the solver
   */
  template< typename FIELD_OP, typename POLICY, typename LAMBDA >
  void
  computeRhsContribution( SortedArrayView< localIndex const > const & targetSet,
                          real64 const time,
                          real64 const dt,
                          dataRepository::Group const & dataGroup,
                          arrayView1d< globalIndex const > const & dofMap,
                          globalIndex const dofRankOffset,
                          CRSMatrixView< real64, globalIndex const > const & matrix,
                          arrayView1d< globalIndex > const & dof,
                          arrayView1d< real64 > const & rhsContribution,
                          LAMBDA && lambda ) const;


  /**
   * @brief Function to zero matrix rows to apply boundary conditions
   * @tparam POLICY the execution policy to use when zeroing rows
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[inout] matrix the local system matrix
   *
   * This function zeroes the rows of the matrix that correspond to boundary conditions.
   */
  template< typename POLICY >
  void zeroSystemRowsForBoundaryCondition( SortedArrayView< localIndex const > const & targetSet,
                                           arrayView1d< globalIndex const > const & dofMap,
                                           CRSMatrixView< real64, globalIndex const > const & matrix ) const;

  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    /// @return The key for setName
    constexpr static char const * setNamesString() { return "setNames"; }
    /// @return The key for constitutivePath
    constexpr static char const * constitutivePathString() { return "constitutivePath"; }
    /// @return The key for objectPath
    constexpr static char const * objectPathString() { return "objectPath"; }
    /// @return The key for fieldName
    constexpr static char const * fieldNameString() { return "fieldName"; }
    /// @return The key for dataType
    constexpr static char const * dataTypeString() { return "dataType"; }
    /// @return The key for component
    constexpr static char const * componentString() { return "component"; }
    /// @return The key for direction
    constexpr static char const * directionString() { return "direction"; }
    /// @return The key for bcApplicationTableName
    constexpr static char const * bcApplicationTableNameString() { return "bcApplicationTableName"; }
    /// @return The key for scale
    constexpr static char const * scaleString() { return "scale"; }
    /// @return The key for functionName
    constexpr static char const * functionNameString() { return "functionName"; }
    /// @return The key for initialCondition
    constexpr static char const * initialConditionString() { return "initialCondition"; }
    /// @return The key for beginTime
    constexpr static char const * beginTimeString() { return "beginTime"; }
    /// @return The key for endTime
    constexpr static char const * endTimeString() { return "endTime"; }
  };

  /**
   * Accessor
   * @return const reference to m_function
   */
  string const & getFunctionName() const
  {
    return m_functionName;
  }

  /**
   * Accessor
   * @return const reference to m_objectPath
   */
  virtual const string & getObjectPath() const
  {
    return m_objectPath;
  }

  /**
   * Accessor
   * @return const reference to m_fieldName
   */
  virtual const string & getFieldName() const
  {
    return m_fieldName;
  }

  /**
   * Accessing the considered component.
   * @return The component axis or a special value.
   */
  virtual int getComponent() const
  {
    return m_component;
  }

  /**
   * Accessor
   * @return const reference to m_direction
   */
  virtual R1Tensor const & getDirection() const
  {
    GEOS_UNUSED_VAR( time );
    return m_direction;
  }

  /**
   * Accessor
   * @return const m_beginTime
   */
  real64 getStartTime() const
  {
    return m_beginTime;
  }

  /**
   * Accessor
   * @return const m_endTime
   */
  real64 getEndTime() const
  {
    return m_endTime;
  }

  /**
   * Accessor
   * @return const reference to m_setNames
   */
  string_array const & getSetNames() const
  {
    return m_setNames;
  }

  /**
   * Accessor
   * @return const m_initialCondition
   */
  int initialCondition() const
  {
    return m_initialCondition;
  }

  /**
   * Accessor
   * @return const m_scale
   */
  real64 getScale() const
  {
    return m_scale;
  }

  /**
   * Mutator
   * @param[in] fieldName The name of the field
   */
  void setFieldName( string const & fieldName )
  {
    m_fieldName = fieldName;
  }

  /**
   * Mutator
   * @param[in] objectPath The path for the object
   */
  void setObjectPath( string const & objectPath )
  {
    m_objectPath = objectPath;
  }

  /**
   * Mutator
   * @param[in] scale Scaling factor
   */
  void setScale( real64 const & scale )
  {
    m_scale = scale;
  }

  /**
   * Mutator
   * @param[in] isInitialCondition Logical value to indicate if it is an initial condition
   */
  void initialCondition( bool isInitialCondition )
  {
    m_initialCondition = isInitialCondition;
  }

  /**
   * Mutator
   * @param[in] setName The name of the set
   */
  void addSetName( string const & setName )
  {
    m_setNames.emplace_back( setName );
  }

  /**
   * @brief Set the Mesh Object Path object
   *
   * @param meshBodies The group containing all the MeshBody objects
   */
  void setMeshObjectPath( Group const & meshBodies );

  /**
   * @brief Get the Mesh Object Paths object
   *
   * @return reference to const m_meshObjectPaths
   */
  MeshObjectPath const & getMeshObjectPaths() const
  {
    return *(m_meshObjectPaths.get());
  }


protected:


private:


  /// the names of the sets that the boundary condition is applied to
  string_array m_setNames;

  /// the path to the object which contains the fields that the boundary condition is applied to
  string m_objectPath;

  std::unique_ptr< MeshObjectPath > m_meshObjectPaths;

  /// the name of the field the boundary condition is applied to or a key string to use for
  /// determining whether or not to apply the boundary condition.
  string m_fieldName;


  /// The component the boundary condition acts on. Not used if field is a scalar.
  int m_component;

  /// The direction the boundary condition acts in.
  R1Tensor m_direction;

  /// Whether or not the boundary condition is an initial condition.
  int m_initialCondition;

  /// The name of the function used to generate values for application.
  string m_functionName;

  /// The scale factor to use on the value of the boundary condition.
  real64 m_scale;

  /// Time after which the bc is allowed to be applied
  real64 m_beginTime;

  /// Time after which the bc will no longer be applied.
  real64 m_endTime;

  /// The name of a function used to turn on and off the boundary condition.
  string m_bcApplicationFunctionName;

};


template< typename FIELD_OP, typename POLICY, typename T, int N, int USD >
void FieldSpecificationBase::applyFieldValueKernel( ArrayView< T, N, USD > const & field,
                                                    SortedArrayView< localIndex const > const & targetSet,
                                                    real64 const time,
                                                    Group & dataGroup ) const
{
  integer const component = getComponent();
  FunctionManager & functionManager = FunctionManager::getInstance();

  if( m_functionName.empty() )
  {
    real64 const value = m_scale;
    forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      FIELD_OP::SpecifyFieldValue( field, a, component, value );
    } );
  }
  else
  {
    FunctionBase const & function = [&]() -> FunctionBase const &
    {
      try
      {
        return functionManager.getGroup< FunctionBase >( m_functionName );
      }
      catch( std::exception const & e )
      {
        throw InputError( e, GEOS_FMT( "Error while reading {}:\n",
                                       getWrapperDataContext( viewKeyStruct::functionNameString() ) ) );
      }
    }();

    if( function.isFunctionOfTime()==2 )
    {
      real64 const value = m_scale * function.evaluate( &time );
      forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        FIELD_OP::SpecifyFieldValue( field, a, component, value );
      } );
    }
    else
    {
      real64_array result( static_cast< localIndex >( targetSet.size() ) );
      function.evaluate( dataGroup, time, targetSet, result );
      arrayView1d< real64 const > const & resultView = result.toViewConst();
      real64 const scale = m_scale;
      forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        FIELD_OP::SpecifyFieldValue( field, a, component, scale * resultView[i] );
      } );
    }
  }
}


template< typename FIELD_OP, typename POLICY >
void FieldSpecificationBase::applyFieldValue( SortedArrayView< localIndex const > const & targetSet,
                                              real64 const time,
                                              dataRepository::Group & dataGroup,
                                              string const & fieldName ) const
{
  dataRepository::WrapperBase & wrapper = dataGroup.getWrapperBase( fieldName );

  // // This function is used in setting boundary/initial conditions on simulation fields.
  // // This is meaningful for 1/2/3D real arrays and sometimes 1D integer (indicator) arrays.
  using FieldTypes = types::ListofTypeList< types::Join< types::ArrayTypes< types::RealTypes, types::DimsUpTo< 3 > >,
                                                         types::ArrayTypes< types::TypeList< integer >, types::DimsSingle< 1 > > > >;


  types::dispatch( FieldTypes{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    auto & wrapperT = dataRepository::Wrapper< ArrayType >::cast( wrapper );
    applyFieldValueKernel< FIELD_OP, POLICY >( wrapperT.reference().toView(), targetSet, time, dataGroup );
  }, wrapper );
}

template< typename FIELD_OP, typename POLICY, typename T, int NDIM, int USD >
void FieldSpecificationBase::applyBoundaryConditionToSystemKernel( SortedArrayView< localIndex const > const & targetSet,
                                                                   real64 const time,
                                                                   dataRepository::Group const & dataGroup,
                                                                   arrayView1d< globalIndex const > const & dofMap,
                                                                   globalIndex const dofRankOffset,
                                                                   CRSMatrixView< real64, globalIndex const > const & matrix,
                                                                   arrayView1d< real64 > const & rhs,
                                                                   ArrayView< T const, NDIM, USD > const & fieldView ) const
{
  integer const component = getComponent();
  this->applyBoundaryConditionToSystem< FIELD_OP, POLICY >( targetSet, time, dataGroup, dofMap, dofRankOffset, matrix, rhs,
                                                            [fieldView, component] GEOS_HOST_DEVICE ( localIndex const a )
  {
    real64 value = 0.0;
    FieldSpecificationEqual::readFieldValue( fieldView, a, component, value );
    return value;
  } );
}

template< typename FIELD_OP, typename POLICY >
void FieldSpecificationBase::applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                                             real64 const time,
                                                             dataRepository::Group const & dataGroup,
                                                             string const & fieldName,
                                                             string const & dofMapName,
                                                             globalIndex const dofRankOffset,
                                                             CRSMatrixView< real64, globalIndex const > const & matrix,
                                                             arrayView1d< real64 > const & rhs ) const
{
  dataRepository::WrapperBase const & wrapper = dataGroup.getWrapperBase( fieldName );
  arrayView1d< globalIndex const > const & dofMap = dataGroup.getReference< array1d< globalIndex > >( dofMapName );

  // We're reading values from a field, which is only well-defined for dims 1 and 2
  using FieldTypes = types::ListofTypeList< types::ArrayTypes< types::RealTypes, types::DimsUpTo< 2 > > >;
  types::dispatch( FieldTypes{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    auto const & wrapperT = dataRepository::Wrapper< ArrayType >::cast( wrapper );
    applyBoundaryConditionToSystemKernel< FIELD_OP, POLICY >( targetSet,
                                                              time,
                                                              dataGroup,
                                                              dofMap,
                                                              dofRankOffset,
                                                              matrix,
                                                              rhs,
                                                              wrapperT.reference() );
  }, wrapper );
}

template< typename FIELD_OP, typename POLICY, typename LAMBDA >
void
FieldSpecificationBase::
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const
{
  return applyBoundaryConditionToSystem< FIELD_OP, POLICY >( targetSet,
                                                             time,
                                                             1.0,
                                                             dataGroup,
                                                             dofMap,
                                                             dofRankOffset,
                                                             matrix,
                                                             rhs,
                                                             std::forward< LAMBDA >( lambda ) );
}

template< typename FIELD_OP, typename POLICY, typename LAMBDA >
void
FieldSpecificationBase::
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const
{
  array1d< globalIndex > dofArray( targetSet.size() );
  arrayView1d< globalIndex > const & dof = dofArray.toView();

  array1d< real64 > rhsContributionArray( targetSet.size() );
  arrayView1d< real64 > const & rhsContribution = rhsContributionArray.toView();

  computeRhsContribution< FIELD_OP, POLICY, LAMBDA >( targetSet,
                                                      time,
                                                      dt,
                                                      dataGroup,
                                                      dofMap,
                                                      dofRankOffset,
                                                      matrix,
                                                      dof,
                                                      rhsContribution,
                                                      std::forward< LAMBDA >( lambda ) );

  FIELD_OP::template prescribeRhsValues< POLICY >( rhs, dof, dofRankOffset, rhsContribution );
}

template< typename FIELD_OP, typename POLICY, typename LAMBDA >
void
FieldSpecificationBase::
  computeRhsContribution( SortedArrayView< localIndex const > const & targetSet,
                          real64 const time,
                          real64 const dt,
                          dataRepository::Group const & dataGroup,
                          arrayView1d< globalIndex const > const & dofMap,
                          globalIndex const dofRankOffset,
                          CRSMatrixView< real64, globalIndex const > const & matrix,
                          arrayView1d< globalIndex > const & dof,
                          arrayView1d< real64 > const & rhsContribution,
                          LAMBDA && lambda ) const
{
  integer const component = ( getComponent() >=0 ) ? getComponent() : 0;
  string const & functionName = getReference< string >( viewKeyStruct::functionNameString() );
  FunctionManager & functionManager = FunctionManager::getInstance();

  // Compute the value of the rhs terms, and collect the dof numbers
  // The rhs terms will be assembled in applyBoundaryConditionToSystem (or in the solver for CompositionalMultiphaseBase)

  if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
  {
    real64 value = m_scale * dt;
    if( !functionName.empty() )
    {
      FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
      value *= function.evaluate( &time );
    }

    forAll< POLICY >( targetSet.size(),
                      [targetSet, dof, dofMap, dofRankOffset, component, matrix, rhsContribution, value, lambda] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      dof[ i ] = dofMap[ a ] + component;
      FIELD_OP::SpecifyFieldValue( dof[ i ],
                                   dofRankOffset,
                                   matrix,
                                   rhsContribution[ i ],
                                   value,
                                   lambda( a ) );
    } );
  }
  else
  {
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );

    real64_array resultsArray( targetSet.size() );
    function.evaluate( dataGroup, time, targetSet, resultsArray );
    arrayView1d< real64 const > const & results = resultsArray.toViewConst();
    real64 const value = m_scale * dt;

    forAll< POLICY >( targetSet.size(),
                      [targetSet, dof, dofMap, dofRankOffset, component, matrix, rhsContribution, results, value, lambda] GEOS_HOST_DEVICE (
                        localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      dof[ i ] = dofMap[ a ] + component;
      FIELD_OP::SpecifyFieldValue( dof[ i ],
                                   dofRankOffset,
                                   matrix,
                                   rhsContribution[ i ],
                                   value * results[ i ],
                                   lambda( a ) );
    } );
  }
}


template< typename POLICY >
void FieldSpecificationBase::zeroSystemRowsForBoundaryCondition( SortedArrayView< localIndex const > const & targetSet,
                                                                 arrayView1d< globalIndex const > const & dofMap,
                                                                 CRSMatrixView< real64, globalIndex const > const & matrix ) const

{
  integer const component = ( getComponent() >=0 ) ? getComponent() : 0;
  forAll< POLICY >( targetSet.size(), [targetSet, dofMap, matrix, component] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    globalIndex const dof = dofMap[ a ] + component;

    arraySlice1d< real64 > const entries = matrix.getEntries( dof );
    localIndex const numEntries = matrix.numNonZeros( dof );

    for( localIndex j = 0; j < numEntries; ++j )
    {
      entries[ j ] = 0;
    }
  } );
}

}

#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP
