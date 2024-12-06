#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "fileIO/Outputs/CatalystOutput.hpp"
#include "PyCatalystOutputType.hpp"
#include "dataRepository/python/PyGroupType.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PyCatalystOutput is not initialized.", nullptr )

namespace geos
{
namespace python
{

struct PyCatalystOutput
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to CatalystOutput.";

  geos::CatalystOutput * group;
};


static PyObject * PyCatalystOutput_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOS_UNUSED_VAR( args, kwds );
  PyCatalystOutput *self;

  self = (PyCatalystOutput *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->group = nullptr;
  }

  return (PyObject *)self;
}


static PyObject * PyCatalystOutput_repr( PyObject * const obj ) noexcept
{
  PyCatalystOutput const * const pyCatalystOutput = LvArray::python::convert< PyCatalystOutput >( obj, getPyCatalystOutputType() );
  if( pyCatalystOutput == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyCatalystOutput );

  string repr;
  string const path = pyCatalystOutput->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pyCatalystOutput->group) ).name() );
  repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );

}


static PyObject * output( PyCatalystOutput * self, PyObject * args )
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
    std::cout << "Target not found. Impossible output."<< std::endl;
  }
  Py_RETURN_NONE;
}

static PyObject * reinit( PyCatalystOutput * self, PyObject *args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOS_UNUSED_VAR( args );

  self->group->reinit();


  Py_RETURN_NONE;
}


static PyMethodDef PyCatalystOutput_methods[] = {
  { "output", (PyCFunction) output, METH_VARARGS, "wrapper to routine CatalystOutput::execute"},
  { "reinit", (PyCFunction) reinit, METH_VARARGS, "reinitialization function"},
  { nullptr, nullptr, 0, nullptr }      /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PyCatalystOutputType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.CatalystOutput",
  .tp_basicsize = sizeof( PyCatalystOutput ),
  .tp_itemsize = 0,
  .tp_repr = PyCatalystOutput_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PyCatalystOutput::docString,
  .tp_methods = PyCatalystOutput_methods,
  .tp_base = getPyGroupType(),
  .tp_new = PyCatalystOutput_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

PyTypeObject * getPyCatalystOutputType()
{
  return &PyCatalystOutputType;
}

}
}
