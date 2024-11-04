#
# Specify all headers
#
set( functions_headers
     FunctionBase.hpp
     FunctionManager.hpp
     TableFunction.hpp
     PythonFunction.hpp
   )

#
# Specify all sources
#
set( functions_sources
     FunctionBase.cpp
     FunctionManager.cpp
     TableFunction.cpp
     MultivariableTableFunction.cpp
   )

if( ENABLE_MATHPRESSO )
  list( APPEND functions_headers
        SymbolicFunction.hpp
        CompositeFunction.hpp )
  list( APPEND functions_sources
        SymbolicFunction.cpp
        CompositeFunction.cpp )
endif()


set( dependencyList ${parallelDeps} dataRepository )

if( ENABLE_MATHPRESSO )
   list( APPEND dependencyList mathpresso )
endif()

if( ENABLE_PYGEOSX )
  list( APPEND functions_headers
        python/PyPythonFunctionType.hpp )
  list( APPEND functions_sources
        python/PyPythonFunction.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       functions
                 SOURCES    ${functions_sources}
                 HEADERS    ${functions_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( functions PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS functions LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif()
