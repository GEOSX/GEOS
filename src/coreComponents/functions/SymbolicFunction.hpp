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

#ifndef GEOS_FUNCTIONS_SYMBOLICFUNCTION_HPP_
#define GEOS_FUNCTIONS_SYMBOLICFUNCTION_HPP_

#include "FunctionBase.hpp"

#include <mathpresso/mathpresso.h>

namespace geos
{

/**
 * @class SymbolicFunction
 *
 * An interface for an arbitrary symbolic function
 */
class SymbolicFunction : public FunctionBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  SymbolicFunction( const string & name,
                    dataRepository::Group * const parent );

  /**
   * @brief The destructor
   */
  virtual ~SymbolicFunction() override;

  /**
   * @brief The catalog name interface
   * @return name of the TableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "SymbolicFunction"; }

  /**
   * @brief Initialize the table function
   */
  virtual void initializeFunction() override;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  inline void evaluate( dataRepository::Group const & group,
                        real64 const time,
                        SortedArrayView< localIndex const > const & set,
                        arrayView1d< real64 > const & result ) const override final
  {
    FunctionBase::evaluateT< SymbolicFunction >( group, time, set, result );
  }

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function result
   */
  inline real64 evaluate( real64 const * const input ) const override final
  {
    return parserExpression.evaluate( reinterpret_cast< void * >( const_cast< real64 * >(input) ) );
  }


  /**
   * @brief Set the symbolic variable names
   * @param variableNames An array of variable names used in the expression
   */
  void setSymbolicVariableNames( string_array variableNames ) { m_variableNames = std::move( variableNames ); }

  /**
   * @brief Set the symbolic expression
   * @param expression A string containing the symbolic expression
   */
  void setSymbolicExpression( string expression ) { m_expression = std::move( expression ); }



private:
  // Symbolic math driver objects
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;

  /// Symbolic expression variable names
  string_array m_variableNames;

  /// Symbolic expression
  string m_expression;
};


} /* namespace geos */

#endif /* GEOS_FUNCTIONS_SYMBOLICFUNCTION_HPP_ */
