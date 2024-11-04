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
 * @file PythonFunction.hpp
 */

#ifndef GEOS_FUNCTIONS_PYTHONFUNCTION_HPP_
#define GEOS_FUNCTIONS_PYTHONFUNCTION_HPP_

#include <Python.h>
#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "LvArray/src/python/PyFunc.hpp"

#if defined(GEOS_USE_PYGEOSX)
  #include "python/PyPythonFunctionType.hpp"
#endif

namespace std
{
template<>
struct hash< __uint128_t >
{
  size_t operator()( const __uint128_t & x ) const noexcept
  {
    size_t h1 = std::hash< uint64_t >{} (static_cast< uint64_t >(x));
    size_t h2 = std::hash< uint64_t >{} (static_cast< uint64_t >(x >> 64));
    return h1 ^ (h2 * 0x9e3779b97f4a7c15 + 0x7f4a7c15);    // Use a large prime multiplier and a random offset
  }
};
};

namespace geos
{
/**
 * @class PythonFunction
 * @brief A wrapper class to interface Python functions for use in C++ computations.
 *
 * @tparam IN_ARRAY type of input array of coordinates
 * @tparam OUT_ARRAY type of output array of values
 * @tparam OUT_2D_ARRAY type of output array of derivatives
 * @tparam INDEX_T datatype used for indexing in multidimensional space
 */
template< typename INDEX_T = __uint128_t >
class PythonFunction : public dataRepository::Group
{
public:
  /// Number of table dimensions (inputs)
  integer numDims;

  /// Number of operators (interpolated functions, outputs)
  integer numOps;

  /// Number of hypercube vertices
  integer numVerts;

  /// Datatype used for indexing
  typedef INDEX_T longIndex;

  typedef LvArray::python::PythonFunction< array1d< real64 > const &, array1d< real64 > & > PyEvaluateFunction;
  typedef LvArray::python::PythonFunction< array1d< real64 > const &, array2d< real64 > & > PyDerivativeFunction;

  /// Python function reference to be called in evaluate
  // PyObject* py_evaluate_func = nullptr;
  mutable std::unique_ptr< PyEvaluateFunction > py_evaluate_func;

  /// Python function reference to be called for evaluation of derivatives
  // PyObject* py_derivative_func = nullptr;
  mutable std::unique_ptr< PyDerivativeFunction > py_derivative_func;

  /**
   * @brief Construct a new wrapper for python function
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  PythonFunction( const string & name, Group * const parent ):
    dataRepository::Group( name, parent ),
    py_evaluate_func( nullptr ),
    py_derivative_func( nullptr )
  {
    if( !Py_IsInitialized())
      Py_Initialize();
  }
  /**
   * @brief The catalog name interface
   * @return name of the PythonFunction in the FunctionBase catalog
   */
  static string catalogName() { return "PythonFunction"; }
  /**
   * @brief Set Python functions for values and (optionally) for their derivatives
   *
   * @param[in] pyFunc pointer to the function
   * @param[in] pyDerivativeFunc pointer to the function evaluating derivatives
   */
  void setEvaluateFunction( PyObject * pyFunc, PyObject * pyDerivativeFunc = nullptr )
  {
    // Set the evaluation function (required)
    if( PyCallable_Check( pyFunc ))
      py_evaluate_func = std::make_unique< PyEvaluateFunction >( pyFunc );
    else
      GEOS_THROW( GEOS_FMT( "{}: Provided Python function is not callable.", getName()), InputError );

    // Set the derivative function (optional)
    if( pyDerivativeFunc )
    {
      if( PyCallable_Check( pyDerivativeFunc ))
        py_derivative_func = std::make_unique< PyDerivativeFunction >( pyDerivativeFunc );
      else
        GEOS_THROW( GEOS_FMT( "{}: Provided Python derivative function is not callable.", getName()), InputError );
    }
  }
  /**
   * @brief Set dimensions of input and output mapping spaces
   *
   * @param[in] numberOfDimesions dimension of input states
   * @param[in] numberOfOperators dimension of output operators
   */
  void setDimensions( integer numberOfDimesions,
                      integer numberOfOperators )
  {
    numDims = numberOfDimesions;
    numOps = numberOfOperators;
    numVerts = 1 << numDims;
  }
  /**
   * @brief Set the properties of parametrization
   *
   * @param[in] numberOfDimesions dimension of input states
   * @param[in] numberOfOperators dimension of output operators
   * @param[in] axisMinimums array with state axes' minimums
   * @param[in] axisMaximums array with state axes' maximums
   * @param[in] axisPoints array defining number of points along each state axis
   */
  void setAxes( integer numberOfDimesions,
                integer numberOfOperators,
                real64_array axisMinimums,
                real64_array axisMaximums,
                integer_array axisPoints )
  {
    numDims = numberOfDimesions;
    numOps = numberOfOperators;
    numVerts = 1 << numDims;

    m_axisMinimums = axisMinimums;
    m_axisMaximums = axisMaximums;
    m_axisPoints = axisPoints;

    m_axisSteps.resize( numDims );
    m_axisStepInvs.resize( numDims );
    m_axisHypercubeMults.resize( numDims );

    for( integer dim = 0; dim < numDims; dim++ )
    {
      m_axisSteps[dim] = (m_axisMaximums[dim] - m_axisMinimums[dim]) / (m_axisPoints[dim] - 1);
      m_axisStepInvs[dim] = 1 / m_axisSteps[dim];
    }

    m_axisHypercubeMults[numDims - 1] = 1;
    for( integer dim = numDims - 2; dim >= 0; --dim )
    {
      m_axisHypercubeMults[dim] = m_axisHypercubeMults[dim + 1] * (m_axisPoints[dim + 1] - 1);
    }
  }
  /**
   * @brief interpolate all operators values at a given point
   *
   * @param[in] state function arguments
   * @param[out] values function values
   */
  GEOS_HOST_DEVICE
  inline
  void
  evaluate( array1d< real64 > const & state,
            array1d< real64 > & values ) const
  {
    if( py_evaluate_func )
      (*py_evaluate_func)( state, values );
    else
      GEOS_THROW( GEOS_FMT( "{}: Evaluation function not set.", getName()), InputError );
  }

  /**
   * @brief interpolate all operators values at a given point
   *
   * @param[in] state function arguments
   * @param[out] values function values
   * @param[out] derivatives derivatives of function values w.r.t. arguments
   */
  GEOS_HOST_DEVICE
  inline
  void
  evaluate( array1d< real64 > const & state,
            array1d< real64 > & values,
            array2d< real64 > & derivatives ) const
  {
    if( py_evaluate_func )
      (*py_evaluate_func)( state, values );
    else
      GEOS_THROW( GEOS_FMT( "{}: Evaluation function not set.", getName()), InputError );

    if( py_derivative_func )
      (*py_derivative_func)( state, derivatives );
    else
      GEOS_THROW( GEOS_FMT( "{}: Derivative function not set.", getName()), InputError );
  }
  /**
   * @brief Retrieves constant point data.
   * @return const unordered_map<longIndex, array1d<real64>>&
   *         A constant reference to the point data mapping.
   */
  GEOS_HOST_DEVICE
  inline
  unordered_map< longIndex, array1d< real64 > > const &
  getPointData() const { return m_pointData; }
  /**
   * @brief Retrieves modifiable point data.
   * @return unordered_map<longIndex, array1d<real64>>&
   *         A reference to the point data mapping, which can be modified by the caller.
   */
  GEOS_HOST_DEVICE
  inline
  unordered_map< longIndex, array1d< real64 > > &
  getPointData() { return m_pointData; }
  /**
   * @brief Retrieves constant hypercube data.
   * @return const unordered_map<longIndex, array1d<real64>>&
   *         A constant reference to the hypercube data mapping.
   */
  GEOS_HOST_DEVICE
  inline
  unordered_map< longIndex, array1d< real64 > > const &
  getHypercubeData() const { return m_hypercubeData; }
  /**
   * @brief Retrieves modifiable hypercube data.
   * @return unordered_map<longIndex, array1d<real64>>&
   *         A reference to the hypercube data mapping, which can be modified by the caller.
   */
  GEOS_HOST_DEVICE
  inline
  unordered_map< longIndex, array1d< real64 > > &
  getHypercubeData() { return m_hypercubeData; }
  /**
   * @brief Get the axes minimums
   * @return a reference to an array of axes minimums
   */
  arrayView1d< real64 const > getAxisMinimums() const { return m_axisMinimums.toViewConst(); }
  /**
   * @brief Get the axes maximums
   * @return a reference to an array of axes maximums
   */
  arrayView1d< real64 const > getAxisMaximums() const { return m_axisMaximums.toViewConst(); }
  /**
   * @brief Get the axes discretization points numbers
   * @return a reference to an array of axes discretization points numbers
   */
  arrayView1d< integer const > getAxisPoints() const { return m_axisPoints.toViewConst(); }
  /**
   * @brief Get the axes step sizes
   * @return a reference to an array of axes step sizes
   */
  arrayView1d< real64 const > getAxisSteps() const { return m_axisSteps.toViewConst(); }
  /**
   * @brief Get the axes step sizes inversions
   * @return a reference to an array of axes step sizes inversions
   */
  arrayView1d< real64 const > getAxisStepInvs() const { return m_axisStepInvs.toViewConst(); }
  /**
   * @brief Get the table axes hypercube index multiplicators
   * @return a reference to an array of table axes hypercube index multiplicators
   */
  arrayView1d< __uint128_t const > getAxisHypercubeMults() const { return m_axisHypercubeMults.toViewConst(); }

#if defined(GEOS_USE_PYGEOSX)
  /**
   * @brief Return PySolver type.
   * @return Return PySolver type.
   */
  virtual PyTypeObject * getPythonType() const override
  {
    return python::getPyPythonFunctionType();
  }
#endif
private:
  /// Array [numDims] of axis minimum values
  array1d< real64 > m_axisMinimums;

  /// Array [numDims] of axis maximum values
  array1d< real64 > m_axisMaximums;

  /// Array [numDims] of axis discretization points numbers
  array1d< integer > m_axisPoints;

  /// Array [numDims] of axis interval lengths (axes are discretized uniformly)
  array1d< real64 > m_axisSteps;

  /// Array [numDims] of inversions of axis interval lengths (axes are discretized uniformly)
  array1d< real64 > m_axisStepInvs;

  /// Array [numDims] of hypercube index mult factors for each axis
  array1d< __uint128_t > m_axisHypercubeMults;

  /**
   * @brief adaptive point storage: the values of operators at requested supporting points
   * Storage is grown dynamically in the process of simulation.
   * Only supporting points that are required for interpolation are computed and added
   */
  unordered_map< longIndex, array1d< real64 > > m_pointData;
  /**
   * @brief adaptive hypercube storage: the values of operators at every vertex of reqested hypercubes
   * Storage is grown dynamically in the process of simulation
   * Only hypercubes that are required for interpolation are computed and added
   *
   * In fact it is an excess storage used to reduce memory accesses during interpolation.
   * Here all values of all vertexes of requested hypercube are stored consecutevely and are accessed via a single index
   * Usage of point_data for interpolation directly would require N_VERTS memory accesses (>1000 accesses for 10-dimensional space)
   *  *
   */
  unordered_map< longIndex, array1d< real64 > > m_hypercubeData;
};

} /* namespace geos */

#endif /* GEOS_FUNCTIONS_PYTHONFUNCTION_HPP_ */
