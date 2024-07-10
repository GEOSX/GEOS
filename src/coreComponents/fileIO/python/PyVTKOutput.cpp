#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "fileIO/Outputs/VTKOutput.hpp"
#include "PyVTKOutputType.hpp"
#include "dataRepository/python/PyGroupType.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PyVTKOutput is not initialized.", nullptr )

namespace geos
{
namespace python
{

struct PyVTKOutput
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to VTKOutput.";

  geos::VTKOutput * group;
};


static PyObject * PyVTKOutput_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOS_UNUSED_VAR( args, kwds );
  PyVTKOutput *self;

  self = (PyVTKOutput *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->group = nullptr;
  }

  return (PyObject *)self;
}


static PyObject * PyVTKOutput_repr( PyObject * const obj ) noexcept
{
  PyVTKOutput const * const pyVTKOutput = LvArray::python::convert< PyVTKOutput >( obj, getPyVTKOutputType() );
  if( pyVTKOutput == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyVTKOutput );

  string repr;
  string const path = pyVTKOutput->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pyVTKOutput->group) ).name() );
  repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );

}


static PyObject * output( PyVTKOutput * self, PyObject * args )
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

static PyObject * setOutputDir( PyVTKOutput * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  PyObject * unicodePath;

  if( !PyArg_ParseTuple( args, "U", &unicodePath ) )
  {
    return nullptr;
  }

  LvArray::python::PyObjectRef<> asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if( asciiPath == nullptr )
  {
    return nullptr;
  }

  char const * const path = PyBytes_AsString( asciiPath );
  if( path == nullptr )
  {
    return nullptr;
  }

  try
  {
    self->group->setOutputDirectory( path );
    self->group->postInputInitialization();
  }
  catch( std::out_of_range const & e )
  {
    std::cout << "Target not found." <<std::endl;
  }
  Py_RETURN_NONE;
}

static PyObject * setOutputFileRootName( PyVTKOutput * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  PyObject * fileNameRoot;

  if( !PyArg_ParseTuple( args, "U", &fileNameRoot ) )
  {
    return nullptr;
  }

  LvArray::python::PyObjectRef<> asciiFileNameRoot { PyUnicode_AsASCIIString( fileNameRoot ) };
  if( asciiFileNameRoot == nullptr )
  {
    return nullptr;
  }

  char const * const filename = PyBytes_AsString( asciiFileNameRoot );
  if( filename == nullptr )
  {
    return nullptr;
  }

  try
  {
    self->group->setPlotFileRoot( filename );
    self->group->postInputInitialization();
  }
  catch( std::out_of_range const & e )
  {
    std::cout << "Target not found." << std::endl;
  }
  Py_RETURN_NONE;
}


static PyObject * reinit( PyVTKOutput * self, PyObject *args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOS_UNUSED_VAR( args );

  self->group->reinit();


  Py_RETURN_NONE;
}


static PyMethodDef PyVTKOutput_methods[] = {
  { "output", (PyCFunction) output, METH_VARARGS, "wrapper to routine VTKOutput::execute"},
  { "setOutputDir", (PyCFunction) setOutputDir, METH_VARARGS, "wrapper to change directory output"},
  { "setOutputFileRootName", (PyCFunction) setOutputFileRootName, METH_VARARGS, "wrapper to change the root name of the output files and subfolders"},
  { "reinit", (PyCFunction) reinit, METH_VARARGS, "reinitialization function"},
  { nullptr, nullptr, 0, nullptr }      /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PyVTKOutputType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.VTKOutput",
  .tp_basicsize = sizeof( PyVTKOutput ),
  .tp_itemsize = 0,
  .tp_repr = PyVTKOutput_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PyVTKOutput::docString,
  .tp_methods = PyVTKOutput_methods,
  .tp_base = getPyGroupType(),
  .tp_new = PyVTKOutput_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

PyTypeObject * getPyVTKOutputType()
{
  return &PyVTKOutputType;
}

}
}
