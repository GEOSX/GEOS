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
 * @file MultivariableTableFunction.hpp
 */

#ifndef GEOS_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_
#define GEOS_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

#include "codingUtilities/EnumStrings.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

/**
 * @class MultivariableTableFunction
 *
 * An interface class for multivariable table function (function with multiple inputs and outputs) with uniform discretization
 * Prepares input data for MultilinearInterpolatorStaticKernel, which performes actual interpolation
 */

class MultivariableTableFunction : public FunctionBase
{
public:

  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  MultivariableTableFunction( const string & name,
                              Group * const parent );

  /**
   * @brief The catalog name interface
   * @return name of the MultivariableTableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "MultivariableTableFunction"; }

  /**
   * @brief Set table coordinates
   *
   * @param[in] numDims number of table dimensions (number of inputs)
   * @param[in] numOps number of functions to interpolate (number of outputs)
   * @param[in] axisMinimums minimum coordinate for each axis
   * @param[in] axisMaximums maximum coordinate for each axis
   * @param[in] axisPoints number of discretization points between minimum and maximum for each axis
   */
  void setTableCoordinates( integer const numDims,
                            integer const numOps,
                            real64_array const & axisMinimums,
                            real64_array const & axisMaximums,
                            integer_array const & axisPoints );

  /**
   * @brief Set the table values
   * @param values An array of table values in C order (the most rapidly changing index is the last)
   */
  void setTableValues( real64_array const values );



  /**
   * @brief Initialize the table function after setting table coordinates and values
   */
  virtual void initializeFunction() override;

  /**
   * @brief Initialize the table function using data from file
   * @param[in] filename The name of the file to read.
   */
  void initializeFunctionFromFile( string const & filename );


  /**
   * @brief Method to evaluate a function on a target object (not supported)
   * @param[in] group a pointer to the object holding the function arguments
   * @param[in] time current time
   * @param[in] set the subset of nodes to apply the function to
   * @param[out] result an array to hold the results of the function
   */
  void evaluate( dataRepository::Group const & group,
                 real64 const time,
                 SortedArrayView< localIndex const > const & set,
                 arrayView1d< real64 > const & result ) const override final
  {
    GEOS_UNUSED_VAR( group );
    GEOS_UNUSED_VAR( time );
    GEOS_UNUSED_VAR( set );
    GEOS_UNUSED_VAR( result );

    GEOS_ERROR( "This method is not supported by MultivariableTableFunction" );

  };

  /**
   * @brief Method to evaluate a function (not supported)
   * @param[in] input a scalar input
   * @return the function evaluation
   */
  virtual real64 evaluate( real64 const * const input ) const override final
  {
    GEOS_UNUSED_VAR( input );

    GEOS_ERROR( "This method is not supported by MultivariableTableFunction" );
    return 0;
  };

  /**
   * @brief Get the table axes minimums
   * @return a reference to an array of table axes minimums
   */
  arrayView1d< real64 const > getAxisMinimums() const { return m_axisMinimums.toViewConst(); }

  /**
   * @brief Get the table axes maximums
   * @return a reference to an array of table axes maximums
   */
  arrayView1d< real64 const > getAxisMaximums() const { return m_axisMaximums.toViewConst(); }

  /**
   * @brief Get the table axes discretization points numbers
   * @return a reference to an array of table axes discretization points numbers
   */
  arrayView1d< integer const > getAxisPoints() const { return m_axisPoints.toViewConst(); }


  /**
   * @brief Get the table axes step sizes
   * @return a reference to an array of table axes step sizes
   */
  arrayView1d< real64 const > getAxisSteps() const { return m_axisSteps.toViewConst(); }

  /**
   * @brief Get the table axes step sizes inversions
   * @return a reference to an array of table axes step sizes inversions
   */
  arrayView1d< real64 const > getAxisStepInvs() const { return m_axisStepInvs.toViewConst(); }

  /**
   * @brief Get the table axes hypercube index multiplicators
   * @return a reference to an array of table axes hypercube index multiplicators
   */
  arrayView1d< __uint128_t const > getAxisHypercubeMults() const { return m_axisHypercubeMults.toViewConst(); }

  /**
   * @brief Get the table values stored per-hypercube
   * @return a reference to an array of table values stored per-hypercube
   */
  arrayView1d< real64 const > getHypercubeData() const { return m_hypercubeData.toViewConst(); }

  /**
   * @brief Get the number of table dimensions
   * @return the number of table dimensions
   */
  integer  numDims() const {return m_numDims;};

  /**
   * @brief Get the number of operators (functions to be interpolated)
   * @return the number of operators (functions to be interpolated)
   */
  integer  numOps() const {return m_numOps;};

private:

  /**
   * @brief Get indexes of all vertices of a hypercube
   *
   * @param[in] hypercubeIndex index of a hypercube
   * @param[out] hypercubePoints array of indexes of hypercube vertexes
   */
  void getHypercubePoints( globalIndex const hypercubeIndex, globalIndex_array & hypercubePoints ) const;

  /// Number of table dimensions (inputs)
  integer m_numDims;

  /// Number of operators (functions to interpolate - ouputs)
  integer m_numOps;

  /// Number of vertices in each hypercube
  integer m_numVerts;

  /// Array [numDims] of axis minimum values
  real64_array m_axisMinimums;

  /// Array [numDims] of axis maximum values
  real64_array m_axisMaximums;

  /// Array [numDims] of axis discretization points numbers
  integer_array m_axisPoints;

  // inputs : data derived from table coordinates data

  ///  Array [numDims] of axis interval lengths (axes are discretized uniformly)
  real64_array m_axisSteps;

  ///  Array [numDims] of inversions of axis interval lengths (axes are discretized uniformly)
  real64_array m_axisStepInvs;

  ///  Array [numDims] of point index mult factors for each axis
  array1d< __uint128_t > m_axisPointMults;

  ///  Array [numDims] of hypercube index mult factors for each axis
  array1d< __uint128_t > m_axisHypercubeMults;

  // inputs: operator sample data

  ///  Main table data stored per point: all operator values for each table coordinate
  real64_array m_pointData;

  ///  Main table data stored per hypercube: all values required for interpolation withing give hypercube are stored contiguously
  real64_array m_hypercubeData;
};


} /* namespace geos */

#endif /* GEOS_FUNCTIONS_MULTIVARIABLETABLEFUNCTION_HPP_ */
