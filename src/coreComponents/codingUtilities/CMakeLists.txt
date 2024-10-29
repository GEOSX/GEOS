#
# Specify all headers
#
set( codingUtilities_headers
     EnumStrings.hpp
     Parsing.hpp
     SFINAE_Macros.hpp
     StringUtilities.hpp
     UnitTestUtilities.hpp
     Utilities.hpp
     traits.hpp
   )

#
# Specify all sources
#
set( codingUtilities_sources
     Parsing.cpp
     StringUtilities.cpp
     )

set( dependencyList ${parallelDeps} common fast_float )

blt_add_library( NAME       codingUtilities
                 SOURCES    ${codingUtilities_sources}
                 HEADERS    ${codingUtilities_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )

# Avoid compiling with nvcc which sometimes crashes on fast_float
set_source_files_properties( Parsing.cpp PROPERTIES LANGUAGE CXX )

target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )


if( GEOS_ENABLE_TESTS )
  add_subdirectory(tests)
endif( )
