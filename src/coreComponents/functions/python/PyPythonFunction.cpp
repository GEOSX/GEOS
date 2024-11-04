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

// Python.h must be included first.
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "LvArray/src/python/pythonHelpers.hpp"
#include "functions/PythonFunction.hpp"
#include "PyPythonFunctionType.hpp"
#include "dataRepository/python/PyGroupType.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PyPythonFunction is not initialized.", nullptr )

namespace geos
{
namespace python
{

struct PyPythonFunction
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geos::PythonFunction.";

  geos::PythonFunction<> * group;
};


/**
 * @brief Allocate a new PyPythonFunction object
 */
static PyObject * PyPythonFunction_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOS_UNUSED_VAR( args, kwds );
  PyPythonFunction *self;

  // Allocate memory for the Python object
  self = (PyPythonFunction *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->group = nullptr;  // Initialize the pyFunc to nullptr (not yet assigned)
  }

  return (PyObject *)self;
}
/**
 * @brief String representation of PyPythonFunction object
 */
static PyObject * PyPythonFunction_repr( PyObject * const obj ) noexcept
{
  PyPythonFunction const * const pyPythonFunc = reinterpret_cast< PyPythonFunction * >( obj );
  if( pyPythonFunc == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyPythonFunc );

  std::string const path = pyPythonFunc->group->getPath();  // Assuming PythonFunction has getPath()
  std::string const type = LvArray::system::demangle( typeid( *(pyPythonFunc->group) ).name() );
  std::string const repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );
}
/**
 * @brief Expose the setAxes function to Python
 */
static PyObject * PyPythonFunction_setAxes( PyPythonFunction * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  PythonFunction<> * func = self->group;

  if( func == nullptr )
  {
    PyErr_SetString( PyExc_RuntimeError, "Not a valid PythonFunction instance." );
    return nullptr;
  }

  // Parse arguments
  integer numDims, numOps;
  PyObject *axisMinList, *axisMaxList, *axisPointsList;

  if( !PyArg_ParseTuple( args, "iiOOO", &numDims, &numOps, &axisMinList, &axisMaxList, &axisPointsList ))
  {
    return nullptr;
  }

  // Check list lengths
  if( PyList_Size( axisMinList ) != numDims || PyList_Size( axisMaxList ) != numDims || PyList_Size( axisPointsList ) != numDims )
  {
    PyErr_SetString( PyExc_ValueError, "List lengths do not match numDims." );
    return nullptr;
  }

  // Convert Python lists to internal array representations
  real64_array axisMinimums( numDims );
  real64_array axisMaximums( numDims );
  integer_array axisPoints( numDims );

  for( integer i = 0; i < numDims; ++i )
  {
    axisMinimums[i] = PyFloat_AsDouble( PyList_GetItem( axisMinList, i ));
    axisMaximums[i] = PyFloat_AsDouble( PyList_GetItem( axisMaxList, i ));
    axisPoints[i] = PyLong_AsLong( PyList_GetItem( axisPointsList, i ));
  }

  // Call the C++ function to set axes
  func->setAxes( numDims, numOps, axisMinimums, axisMaximums, axisPoints );

  Py_RETURN_NONE;
}

/**
 * @brief Expose the setEvaluateFunction function to Python
 */
static PyObject * PyPythonFunction_setEvaluateFunction( PyPythonFunction * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  PythonFunction<> * func = self->group;

  if( func == nullptr )
  {
    PyErr_SetString( PyExc_RuntimeError, "Not a valid PythonFunction instance." );
    return nullptr;
  }

  PyObject * pyFuncObj;

  if( !PyArg_ParseTuple( args, "O", &pyFuncObj ))
  {
    return nullptr;
  }

  // Ensure pyFuncObj is callable
  if( !PyCallable_Check( pyFuncObj ))
  {
    PyErr_SetString( PyExc_TypeError, "Argument must be callable." );
    return nullptr;
  }

  // Call the C++ function to set the evaluate function
  func->setEvaluateFunction( pyFuncObj );

  Py_RETURN_NONE;
}

static PyMethodDef PyPythonFunction_methods[] = {
  { "setAxes", (PyCFunction) PyPythonFunction_setAxes, METH_VARARGS, "Set axes for the Python function." },
  { "setEvaluateFunction", (PyCFunction) PyPythonFunction_setEvaluateFunction, METH_VARARGS, "Set Python evaluate function and optionally a derivative function." },
  { nullptr, nullptr, 0, nullptr }  // Sentinel
};

/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PyPythonFunctionType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.PythonFunction",
  .tp_basicsize = sizeof(PyPythonFunction),
  .tp_repr = PyPythonFunction_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PyPythonFunction::docString,
  .tp_methods = PyPythonFunction_methods,
  .tp_base = getPyGroupType(),
  .tp_new = PyPythonFunction_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

/**
 * Return the PyTypeObject for PyPythonFunction
 */
PyTypeObject * getPyPythonFunctionType()
{
  return &PyPythonFunctionType;
}

} // namespace python
} // namespace geos
