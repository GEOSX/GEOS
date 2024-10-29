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

#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "fileIO/timeHistory/HistoryCollection.hpp"

#include "PyHistoryCollectionType.hpp"
#include "dataRepository/python/PyGroupType.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PyHistoryCollection is not initialized.", nullptr )


namespace geos
{
namespace python
{


struct PyHistoryCollection
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to HistoryCollection.";

  geos::HistoryCollection * group;
};


static PyObject * PyHistoryCollection_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOS_UNUSED_VAR( args, kwds );
  PyHistoryCollection *self;

  self = (PyHistoryCollection *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->group = nullptr;
  }

  return (PyObject *)self;
}


static PyObject * PyHistoryCollection_repr( PyObject * const obj ) noexcept
{
  PyHistoryCollection const * const pyHistoryCollection = LvArray::python::convert< PyHistoryCollection >( obj, getPyHistoryCollectionType() );
  if( pyHistoryCollection == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyHistoryCollection );

  string repr;
  string const path = pyHistoryCollection->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pyHistoryCollection->group) ).name() );
  repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );

}



static PyObject * collect( PyHistoryCollection * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;

  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geos::DomainPartition & domain = self->group->getGroupByPath< DomainPartition >( "/Problem/domain" );

  int cycleNumber = int(round( time/dt ));
  try
  {
    self->group->execute( time, dt, cycleNumber, 0, 0, domain );
  }
  catch( std::out_of_range const & e )
  {
    std::cout<<"Target not found. Impossible output."<<std::endl;
  }
  Py_RETURN_NONE;
}

static PyMethodDef PyHistoryCollection_methods[] = {
  { "collect", (PyCFunction) collect, METH_VARARGS, "wrapper to routine HistoryCollection::execute" },
  { nullptr, nullptr, 0, nullptr }      /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PyHistoryCollectionType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.HistoryCollection",
  .tp_basicsize = sizeof( PyHistoryCollection ),
  .tp_itemsize = 0,
  .tp_repr = PyHistoryCollection_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PyHistoryCollection::docString,
  .tp_methods = PyHistoryCollection_methods,
  .tp_base = getPyGroupType(),
  .tp_new = PyHistoryCollection_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

PyTypeObject * getPyHistoryCollectionType()
{ return &PyHistoryCollectionType; }

}
}
