/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FunctionBase.hpp
 */

#ifndef GEOS_FUNCTIONS_FUNCTIONBASE_HPP_
#define GEOS_FUNCTIONS_FUNCTIONBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/TypeDispatch.hpp"
#include "dataRepository/Group.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{
/**
 * @brief The key for inputVarNames
 * @return the key
 */
string const inputVarNames( "inputVarNames" );
}
}

/**
 * @class FunctionBase
 *
 * An object for interfacing with arbitrary N-dimensional functions.
 */
class FunctionBase : public dataRepository::Group
{
public:

  /// Maximum total number of independent variables (including components of multidimensional variables)
  static constexpr int MAX_VARS = 4;

  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  FunctionBase( const string & name,
                dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  virtual ~FunctionBase() override = default;

  /**
   * @brief Function initialization
   */
  virtual void initializeFunction() = 0;

  /**
   * @brief Test to see if the function is a 1D function of time
   * @return integer value:
   *         - 0 is the function does not have time as parameter
   *         - 1 is the function has time as one of the parameters
   *         - 2 is the function has time as only parameter
   */
  integer isFunctionOfTime() const;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void evaluate( dataRepository::Group const & group,
                         real64 const time,
                         SortedArrayView< localIndex const > const & set,
                         arrayView1d< real64 > const & result ) const = 0;

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function evaluation
   */
  virtual real64 evaluate( real64 const * const input ) const = 0;

  /// Alias for the catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< FunctionBase, string const &, Group * const >;

  /**
   * @brief return the catalog entry for the function
   * @return the catalog entry
   */
  static CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  /**
   * @brief This generates statistics by applying a function to an object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @return An array holding the min, average, max values of the results
   */
  real64_array evaluateStats( dataRepository::Group const & group,
                              real64 const time,
                              SortedArray< localIndex > const & set ) const;

  /**
   * @brief Set the input variable names
   * @param inputVarNames A list of input variable names
   */
  void setInputVarNames( string_array inputVarNames ) { m_inputVarNames = std::move( inputVarNames ); }


protected:
  /// names for the input variables
  string_array m_inputVarNames;

  /**
   * @brief Method to apply an function with an arbitrary type of output
   * @tparam LEAF the return type
   * @param[in] group a pointer to the object holding the function arguments
   * @param[in] time current time
   * @param[in] set the subset of nodes to apply the function to
   * @param[out] result the results
   */
  template< typename LEAF, typename POLICY = serialPolicy >
  void evaluateT( dataRepository::Group const & group,
                  real64 const time,
                  SortedArrayView< localIndex const > const & set,
                  arrayView1d< real64 > const & result ) const;

  virtual void postInputInitialization() override { initializeFunction(); }

};

template< typename LEAF, typename POLICY >
void FunctionBase::evaluateT( dataRepository::Group const & group,
                              real64 const time,
                              SortedArrayView< localIndex const > const & set,
                              arrayView1d< real64 > const & result ) const
{
  real64 const * inputPtrs[MAX_VARS]{};
  localIndex varSize[MAX_VARS]{};
  localIndex varStride[MAX_VARS][2]{};

  integer const numVars = LvArray::integerConversion< integer >( m_inputVarNames.size() );
  localIndex totalVarSize = 0;
  for( integer varIndex = 0; varIndex < numVars; ++varIndex )
  {
    string const & varName = m_inputVarNames[varIndex];

    if( varName == "time" )
    {
      inputPtrs[varIndex] = &time;
      varSize[varIndex] = 1;
    }
    else
    {
      dataRepository::WrapperBase const & wrapper = group.getWrapperBase( varName );
      varSize[varIndex] = wrapper.numArrayComp();

      using Types = types::ListofTypeList< types::ArrayTypes< types::TypeList< real64 >, types::DimsUpTo< 2 > > >;
      types::dispatch( Types{}, [&]( auto tupleOfTypes )
      {
        using ArrayType = camp::first< decltype( tupleOfTypes ) >;
        auto const view = dataRepository::Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();
        view.move( hostMemorySpace, false );
        for( int dim = 0; dim < ArrayType::NDIM; ++dim )
        {
          varStride[varIndex][dim] = view.strides()[dim];
        }
        inputPtrs[varIndex] = view.data();
      }, wrapper );
    }
    totalVarSize += varSize[varIndex];
  }

  // Make sure the inputs do not exceed the maximum length
  GEOS_ERROR_IF_GT_MSG( totalVarSize, MAX_VARS,
                        getDataContext() << ": Function input size exceeded" );

  // Make sure the result / set size match
  GEOS_ERROR_IF_NE_MSG( result.size(), set.size(),
                        getDataContext() << ": To apply a function to a set, the size of the result and set must match" );

  forAll< POLICY >( set.size(), [=]( localIndex const i )
  {
    localIndex const index = set[i];
    real64 input[MAX_VARS]{};
    int offset = 0;
    for( integer varIndex = 0; varIndex < numVars; ++varIndex )
    {
      for( localIndex compIndex = 0; compIndex < varSize[varIndex]; ++compIndex )
      {
        input[offset++] = inputPtrs[varIndex][index * varStride[varIndex][0] + compIndex * varStride[varIndex][1]];
      }
    }
    result[i] = static_cast< LEAF const * >( this )->evaluate( input );
  } );
}
} /* namespace geos */

#endif /* GEOS_FUNCTIONS_FUNCTIONBASE_HPP_ */
