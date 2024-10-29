#
# Specify all headers
#
set( functions_headers
     FunctionBase.hpp
     FunctionManager.hpp
     TableFunction.hpp
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


set( dependencyList ${parallelDeps} codingUtilities dataRepository fileIO )

if( ENABLE_MATHPRESSO )
   list( APPEND dependencyList mathpresso )
endif()

blt_add_library( NAME       functions
                 SOURCES    ${functions_sources}
                 HEADERS    ${functions_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )

target_include_directories( functions PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif()
