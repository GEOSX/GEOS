#
# Specify all headers
#
set( codingUtilities_headers
     RTTypes.hpp
     Parsing.hpp
     SFINAE_Macros.hpp
     UnitTestUtilities.hpp
     Utilities.hpp
     traits.hpp
   )

#
# Specify all sources
#
set( codingUtilities_sources
     Parsing.cpp
     RTTypes.cpp
   )

set( dependencyList ${parallelDeps} common fast_float )
geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       codingUtilities
                 SOURCES    ${codingUtilities_sources}
                 HEADERS    ${codingUtilities_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

# Avoid compiling with nvcc which sometimes crashes on fast_float
set_source_files_properties( Parsing.cpp PROPERTIES LANGUAGE CXX )

target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS codingUtilities LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )


if( GEOS_ENABLE_TESTS )
  add_subdirectory(tests)
endif( )
