
set( python_headers pygeosx.hpp )

set( python_sources pygeosx.cpp )

set( dependencyList ${extraComponentsLinkList} ${externalComponentsLinkList} mainInterface Python3::Python )

blt_add_library( NAME pygeosx
                 SOURCES ${python_sources}
                 HEADERS ${python_headers}
                 DEPENDS_ON ${dependencyList}
                 SHARED TRUE
                 CLEAR_PREFIX TRUE
               )

add_custom_command( TARGET pygeosx
                    PRE_BUILD
                    COMMAND ${Python3_EXECUTABLE} -m venv ${CMAKE_BINARY_DIR}/lib/PYGEOSX --system-site-packages )

add_custom_command( TARGET pygeosx
                    POST_BUILD
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/lib/pygeosx.so ${CMAKE_BINARY_DIR}/lib/PYGEOSX/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/ )

add_custom_command( TARGET pygeosx
                    POST_BUILD
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/lib/pylvarray.so ${CMAKE_BINARY_DIR}/lib/PYGEOSX/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/ )

target_include_directories( pygeosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
target_include_directories( pygeosx PUBLIC ${CMAKE_BINARY_DIR}/include )

geosx_add_code_checks( PREFIX pygeosx )

install(TARGETS pygeosx LIBRARY DESTINATION lib)

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

# Python formatting
if ( ENABLE_YAPF )
    set( pygeosx_scripts )
    file( GLOB_RECURSE pygeosx_scripts "*py" )

    # Note: blt throws an error if sources doesn't include a c-file, so include dummy.cpp
    blt_add_code_checks( PREFIX  pygeosx_scripts_yapf_style
                         SOURCES ${pygeosx_scripts} ${CMAKE_SOURCE_DIR}/coreComponents/dummy.cpp
                         YAPF_CFG_FILE ${PROJECT_SOURCE_DIR}/yapf.cfg )

endif()


