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
#include "pygeosx.hpp"
#include "dataRepository/python/PyGroup.hpp"
#include "dataRepository/python/PyGroupType.hpp"
#include "physicsSolvers/python/PySolverType.hpp"
#include "fileIO/python/PyHistoryCollectionType.hpp"
#include "fileIO/python/PyHistoryOutputType.hpp"
#include "fileIO/python/PyVTKOutputType.hpp"
#include "mainInterface/initialization.hpp"
#include "LvArray/src/python/PyArray.hpp"

// System includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#pragma GCC diagnostic pop
#include <chrono>


namespace geos
{

std::unique_ptr< GeosxState > g_state;
bool g_alreadyInitialized = false;

PyObject * init( PyObject * const pyArgv, bool const performSetup, long const pythonMPIRank=-1 )
{
  LvArray::python::PyObjectRef<> iterator{ PyObject_GetIter( pyArgv ) };
  if( iterator == nullptr )
  {
    return nullptr;
  }

  // Parse the python list into a std::vector< std::string >
  std::vector< string > stringArgs;
  for( LvArray::python::PyObjectRef<> item{ PyIter_Next( iterator ) }; item != nullptr; item = PyIter_Next( iterator ) )
  {
    LvArray::python::PyObjectRef<> ascii { PyUnicode_AsASCIIString( item ) };
    if( ascii == nullptr )
    {
      return nullptr;
    }

    char const * const stringValue = PyBytes_AsString( ascii );
    if( stringValue == nullptr )
    {
      return nullptr;
    }

    stringArgs.push_back( stringValue );
  }

  if( PyErr_Occurred() )
  {
    return nullptr;
  }

  std::vector< char * > argv( stringArgs.size() + 1 );
  for( std::size_t i = 0; i < stringArgs.size(); ++i )
  {
    argv[ i ] = const_cast< char * >( stringArgs[ i ].data() );
  }

  // Perform the initial setup if requested.
  if( performSetup )
  {
    basicSetup( argv.size() - 1, argv.data(), false );
    g_alreadyInitialized = true;

    // Verify that the ranks match, there is no recovering from this if incorrect.
    GEOS_ERROR_IF_NE( pythonMPIRank, MpiWrapper::commRank() );
  }

  try
  {
    // Must first delete the existing state.
    g_state = nullptr;

    g_state = std::make_unique< GeosxState >( parseCommandLineOptions( argv.size() - 1, argv.data() ) );
    g_state->initializeDataRepository();
  }
  catch( InputError const & e )
  {
    g_state = nullptr;
    PyErr_SetString( PyExc_RuntimeError, e.what() );
    return nullptr;
  }
  catch( NotAnError const & e )
  {
    g_state = nullptr;
    Py_RETURN_NONE;
  }

  return python::createNewPyGroup( g_state->getProblemManagerAsGroup() );
}

static constexpr char const * initializeDocString =
  "initialize(rank, args)\n"
  "--\n\n"
  "Initialize GEOSX, this must be the first module call.\n\n"
  "Should only be called once. All calls after the first will\n"
  "raise a `RuntimeError`. To reinitialize GEOSX for a new problem,\n"
  "use the `reinit` function.\n"
  "\n"
  "Parameters\n"
  "__________\n"
  "rank : int\n"
  "    The MPI rank of the process.\n"
  "args : list of strings\n"
  "    The list of command line arguments to pass to GEOSX.\n"
  "\n"
  "Returns\n"
  "_______\n"
  "Group\n"
  "    The ProblemManager.";
PyObject * initialize( PyObject * self, PyObject * args ) noexcept
{
  GEOS_UNUSED_VAR( self );

  PYTHON_ERROR_IF( g_alreadyInitialized, PyExc_RuntimeError, "You have already called initialize, call reinitialize.", nullptr );

  long pythonMPIRank;
  PyObject * list;
  if( !PyArg_ParseTuple( args, "lO", &pythonMPIRank, &list ) )
  {
    return nullptr;
  }

  return init( list, true, pythonMPIRank );
}

static constexpr char const * reinitDocString =
  "reinit(args)\n"
  "--\n\n"
  "Reinitialize GEOSX with a new set of command-line arguments.\n"
  "\n"
  "Parameters\n"
  "__________\n"
  "args : list of strings\n"
  "    The list of command line arguments to pass to GEOSX.\n"
  "\n"
  "Returns\n"
  "_______\n"
  "Group\n"
  "    The ProblemManager.";
PyObject * reinit( PyObject * self, PyObject * args ) noexcept
{
  GEOS_UNUSED_VAR( self );

  PyObject * list;
  if( !PyArg_ParseTuple( args, "O", &list ) )
  {
    return nullptr;
  }

  return init( list, false );
}

static constexpr char const * applyInitialConditionsDocString =
  "apply_initial_conditions()\n"
  "--\n\n"
  "Apply the initial conditions.\n"
  "\n"
  "Returns\n"
  "_______\n"
  "None\n";
PyObject * applyInitialConditions( PyObject * self, PyObject * args ) noexcept
{
  GEOS_UNUSED_VAR( self, args );

  PYTHON_ERROR_IF( g_state == nullptr, PyExc_RuntimeError, "state must be initialized", nullptr );

  try
  {
    g_state->applyInitialConditions();
  }
  catch( std::exception const & e )
  {
    PyErr_SetString( PyExc_RuntimeError, e.what() );
    return nullptr;
  }

  Py_RETURN_NONE;
}

static constexpr char const * runDocString =
  "run()\n"
  "--\n\n"
  "Enter the event loop.\n"
  "\n"
  "Returns\n"
  "_______\n"
  "int\n"
  "    The state of the simulation. If the simulation has ended the value is ``COMPLETED``. If the "
  "simulation still has steps left to run the value is ``READY_TO_RUN``.";
PyObject * run( PyObject * self, PyObject * args ) noexcept
{
  GEOS_UNUSED_VAR( self, args );

  PYTHON_ERROR_IF( g_state == nullptr, PyExc_RuntimeError, "state must be initialized", nullptr );

  try
  {
    g_state->run();
  }
  catch( std::exception const & e )
  {
    PyErr_SetString( PyExc_RuntimeError, e.what() );
    return nullptr;
  }

  return PyLong_FromLong( static_cast< int >( g_state->getState() ) );
}

static constexpr char const * getStateDocString =
  "getState()\n"
  "--\n\n"
  "_______\n"
  "int\n"
  "    Return the state of GEOS. 0 : UNINITIALIZED, 1 : INITIALIZED, 2 : READY_TO_RUN, 3 : COMPLETED";

PyObject * getState( PyObject * self, PyObject * args ) noexcept
{
  GEOS_UNUSED_VAR( self, args );

  if( g_state == nullptr )
  {
    return PyLong_FromLong( 0 );
  }
  else
  {
    return PyLong_FromLong( static_cast< int >( g_state->getState() ) );
  }
}

static constexpr char const * finalizeDocString =
  "_finalize()\n"
  "--\n\n"
  "Finalize GEOSX. After this no calls into pygeosx or to MPI are allowed.\n"
  "\n"
  "Returns\n"
  "_______\n"
  "None\n";
PyObject * finalize( PyObject * self, PyObject * args ) noexcept
{
  GEOS_UNUSED_VAR( self, args );

  if( g_state == nullptr )
  {
    Py_RETURN_NONE;
  }

  g_state = nullptr;
  basicCleanup();
  Py_RETURN_NONE;
}

} // namespace geos


/**
 * Add geos::State enums to the given module. Return the module, or nullptr on failure
 */

static bool addConstants( PyObject * module )
{
  std::array< std::pair< long, char const * >, 4 > const constants = { {
    { static_cast< long >( geos::State::COMPLETED ), "COMPLETED" },
    { static_cast< long >( geos::State::INITIALIZED ), "INITIALIZED" },
    { static_cast< long >( geos::State::UNINITIALIZED ), "UNINITIALIZED" },
    { static_cast< long >( geos::State::READY_TO_RUN ), "READY_TO_RUN" }
  } };

  for( std::pair< long, char const * > const & pair : constants )
  {
    PYTHON_ERROR_IF( PyModule_AddIntConstant( module, pair.second, pair.first ), PyExc_RuntimeError,
                     "couldn't add constant", false );
  }

  return true;
}

/**
 * Add exit handler to a module. Register the module's 'finalize' function,
 * which should take no arguments, with the `atexit` standard library module.
 */
static bool addExitHandler( PyObject * module )
{
  LvArray::python::PyObjectRef<> atexit_module { PyImport_ImportModule( "atexit" ) };

  if( atexit_module == nullptr )
  {
    return false;
  }

  LvArray::python::PyObjectRef<> atexit_register_pyfunc { PyObject_GetAttrString( atexit_module, "register" ) };
  if( atexit_register_pyfunc == nullptr )
  {
    return false;
  }

  LvArray::python::PyObjectRef<> finalize_pyfunc { PyObject_GetAttrString( module, "_finalize" ) };
  if( finalize_pyfunc == nullptr )
  {
    return false;
  }

  if( !PyCallable_Check( atexit_register_pyfunc ) || !PyCallable_Check( finalize_pyfunc ) )
  {
    return false;
  }

  LvArray::python::PyObjectRef<> returnval { PyObject_CallFunctionObjArgs( atexit_register_pyfunc, finalize_pyfunc.get(), nullptr ) };

  return returnval != nullptr;
}



BEGIN_ALLOW_DESIGNATED_INITIALIZERS

/**
 *
 */

static PyMethodDef pygeosxFuncs[] = {
  { "initialize", geos::initialize, METH_VARARGS, geos::initializeDocString },
  { "reinit", geos::reinit, METH_VARARGS, geos::reinitDocString },
  { "apply_initial_conditions", geos::applyInitialConditions, METH_NOARGS, geos::applyInitialConditionsDocString },
  { "run", geos::run, METH_NOARGS, geos::runDocString },
  { "getState", geos::getState, METH_NOARGS, geos::getStateDocString },
  { "_finalize", geos::finalize, METH_NOARGS, geos::finalizeDocString },
  { nullptr, nullptr, 0, nullptr }        /* Sentinel */
};



static constexpr char const * pygeosxDocString =
  "Python driver for GEOSX.";

/**
 * Initialize the module object for Python with the exported functions
 */
static struct PyModuleDef pygeosxModuleFunctions = {
  PyModuleDef_HEAD_INIT,
  .m_name = "pygeosx",
  .m_doc = pygeosxDocString,
  .m_size = -1,
  .m_methods = pygeosxFuncs
};


END_ALLOW_DESIGNATED_INITIALIZERS


/**
 * @brief Initialize the module with functions, constants, and exit handler
 */
PyMODINIT_FUNC
PyInit_pygeosx()
{
  import_array();

  LvArray::python::PyObjectRef<> module{ PyModule_Create( &pygeosxModuleFunctions ) };
  if( module == nullptr )
  {
    return nullptr;
  }

  if( !addExitHandler( module ) )
  {
    PYTHON_ERROR_IF( PyErr_Occurred() == nullptr, PyExc_RuntimeError,
                     "couldn't add exit handler", nullptr );

    return nullptr;
  }

  if( !addConstants( module ) )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geos::python::getPyGroupType(), "Group" ) )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geos::python::getPyWrapperType(), "Wrapper" ) )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geos::python::getPySolverType(), "Solver" ) )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geos::python::getPyHistoryCollectionType(), "HistoryCollection" ) )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geos::python::getPyHistoryOutputType(), "HistoryOutput" ) )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geos::python::getPyVTKOutputType(), "VTKOutput" ) )
  {
    return nullptr;
  }

  // Add the LvArray submodule.
  if( !LvArray::python::addPyLvArrayModule( module ) )
  {
    return nullptr;
  }

  // Since we return module we don't want to decrease the reference count.
  return module.release();
}



//=====================================================================================================================================
