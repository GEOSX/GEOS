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
#include "PyWrapper.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->wrapper == nullptr, PyExc_RuntimeError, "The PyWrapper is not initialized.", nullptr )

namespace geos
{
namespace python
{

struct PyWrapper
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geos::dataRepository::WrapperBase.";

  dataRepository::WrapperBase * wrapper;
};

/**
 *
 */
static PyObject * PyWrapper_repr( PyObject * const obj )
{
  PyWrapper const * const self = LvArray::python::convert< PyWrapper >( obj, getPyWrapperType() );
  if( self == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( self );

  string const path = self->wrapper->getPath();
  string const type = LvArray::system::demangle( typeid( *(self->wrapper) ).name() );
  string const repr = path + " ( " + type + " )";
  return PyUnicode_FromString( repr.c_str() );
}

static constexpr char const * PyWrapper_valueDocString =
  "value(self)\n"
  "--\n\n"
  "Return a view of the wrapped value, or ``None`` if it cannot be exported to Python.\n"
  "\n"
  "\n"
  "Returns\n"
  "_______\n"
  "Instance of a pylvarray class\n"
  "    If the wrapped type is one of the LvArray types that have a Python wrapper type.\n"
  "1D numpy.ndarray\n"
  "    If the wrapped type is a numeric constant. The returned array is a shallow copy and "
  "has a single entry.\n"
  "str\n"
  "    If the wrapped type is a string this returns a copy of the string.\n"
  "list of str\n"
  "    If the wrapped type is a `LvArray::Array< string, 1, ... > or a `std::vector< string >`. "
  "This is a copy.\n"
  "None\n"
  "    If the wrapped type is not covered by any of the above.";
static PyObject * PyWrapper_value( PyWrapper * const self, PyObject * const args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  LVARRAY_UNUSED_VARIABLE( args );

  PyObject * const ret = self->wrapper->createPythonObject( );

  // The return value can be a nullptr for two reasons. The first is if the wrapped object is not
  // exportable to Python, in which case we want to return 'None'. The second is if an error occurred
  // in which case we want to propogate it by also returning a nullptr.
  if( ret == nullptr && PyErr_Occurred() == nullptr )
  {
    Py_RETURN_NONE;
  }

  return ret;
}

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef PyWrapperMethods[] = {
  { "value", (PyCFunction) PyWrapper_value, METH_VARARGS, PyWrapper_valueDocString },
  { nullptr, nullptr, 0, nullptr } // Sentinel
};

static PyTypeObject PyWrapperType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.Wrapper",
  .tp_basicsize = sizeof( PyWrapper ),
  .tp_itemsize = 0,
  .tp_repr = PyWrapper_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = PyWrapper::docString,
  .tp_methods = PyWrapperMethods,
  .tp_new = PyType_GenericNew,
};

END_ALLOW_DESIGNATED_INITIALIZERS

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PyObject * createNewPyWrapper( dataRepository::WrapperBase & wrapper )
{
  // Create a new Group and set the dataRepository::Group it points to.
  PyObject * const ret = PyObject_CallFunction( reinterpret_cast< PyObject * >( getPyWrapperType() ), "" );
  PyWrapper * const retWrapper = reinterpret_cast< PyWrapper * >( ret );
  if( retWrapper == nullptr )
  {
    return nullptr;
  }

  retWrapper->wrapper = &wrapper;

  return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PyTypeObject * getPyWrapperType()
{ return &PyWrapperType; }

} // namespace python
} // namespace geos
