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

#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "PySolver.hpp"
#include "dataRepository/python/PyGroupType.hpp"
#include "PySolverType.hpp"
#include "mesh/DomainPartition.hpp"

#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PySolver is not initialized.", nullptr )

namespace geos
{
namespace python
{

struct PySolver
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geos::PhysicsSolverBase.";

  geos::PhysicsSolverBase *group;
};


static PyObject * PySolver_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOS_UNUSED_VAR( args, kwds );
  PySolver *self;

  self = (PySolver *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->group = nullptr;
  }

  return (PyObject *)self;
}


static PyObject * PySolver_repr( PyObject * const obj ) noexcept
{
  PySolver const * const pySolver = LvArray::python::convert< PySolver >( obj, getPySolverType() );
  if( pySolver == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pySolver );

  string const path = pySolver->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pySolver->group) ).name() );
  string const repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );
}



static PyObject * execute( PySolver * self, PyObject * args )
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

  self->group->execute( time, dt, cycleNumber, 0, 0, domain );

  Py_RETURN_NONE;
}


static PyObject * reinit( PySolver * self, PyObject *args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOS_UNUSED_VAR( args );

  self->group->reinit();

  Py_RETURN_NONE;
}

static PyObject * cleanup( PySolver * self, PyObject *args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOS_UNUSED_VAR( args );

  double time;
  if( !PyArg_ParseTuple( args, "d", &time ) )
  {
    return nullptr;
  }

  geos::DomainPartition & domain = self->group->getGroupByPath< DomainPartition >( "/Problem/domain" );
  self->group->cleanup( time, 0, 0, 0.0, domain );

  Py_RETURN_NONE;
}


static PyMethodDef PySolver_methods[] = {
  { "execute", (PyCFunction) execute, METH_VARARGS, "solver Step" },
  { "reinit", (PyCFunction) reinit, METH_NOARGS, "re-initialize certain variable depending on the solver being used"},
  { "cleanup", (PyCFunction) cleanup, METH_VARARGS, "Call cleanup step"},
  { nullptr, nullptr, 0, nullptr }      /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PySolverType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.Solver",
  .tp_basicsize = sizeof( PySolver ),
  .tp_itemsize = 0,
  .tp_repr = PySolver_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PySolver::docString,
  .tp_methods = PySolver_methods,
  .tp_base = getPyGroupType(),
  .tp_new = PySolver_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS


PyTypeObject * getPySolverType()
{ return &PySolverType; }

}
}
