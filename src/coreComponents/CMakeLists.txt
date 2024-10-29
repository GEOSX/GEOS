set( subdirs
     math
     common
     codingUtilities
     dataRepository
     schema
     functions
     constitutive
     mesh
     denseLinearAlgebra
     linearAlgebra
     fieldSpecification
     finiteElement
     finiteVolume
     discretizationMethods
     fileIO
     physicsSolvers
     events
     mainInterface )

unset( parallelDeps )

if ( ENABLE_OPENMP )
  list( APPEND parallelDeps openmp )
endif()

if ( ENABLE_CUDA )
  list( APPEND parallelDeps cuda )
endif()

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND parallelDeps CUDA::nvToolsExt )
endif()

if( ENABLE_HIP )
  list( APPEND parallelDeps blt::hip )
endif()

if( ENABLE_MPI )
  list( APPEND parallelDeps mpi )
endif()

if ( GEOS_BUILD_OBJ_LIBS )
  set( LVARRAY_BUILD_OBJ_LIBS TRUE CACHE BOOL "" FORCE )
endif()

add_subdirectory( LvArray )

if( ENABLE_PYGEOSX )
  list( APPEND coreLibs pylvarray )
endif()

foreach( lib ${subdirs} )
  if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
    list( APPEND coreLibs ${lib} )
  endif()
endforeach()

# if we're building full static libs and not obj_libs
if( NOT GEOS_BUILD_SHARED_LIBS )
  set( geosx_core_list "" )
  foreach( lib ${coreLibs} )
    list( APPEND geosx_core_list ${GEOS_LINK_PREPEND_FLAG} ${lib} ${GEOS_LINK_POSTPEND_FLAG} )
  endforeach()
  set ( geosx_core_libs "${geosx_core_list}" CACHE INTERNAL "" )
endif( )

foreach( lib ${subdirs} )
  add_subdirectory( ${lib} )
endforeach()

if( GEOS_BUILD_SHARED_LIBS AND GEOS_BUILD_OBJ_LIBS )
  message( "Building shared geosx_core library with object coreComponents, executables link to geosx_core" )
  blt_add_library ( NAME geosx_core
                    SOURCES dummy.cpp
                    DEPENDS_ON  mainInterface physicsSolvers HDF5::HDF5
                    SHARED TRUE )

elseif( GEOS_BUILD_SHARED_LIBS AND NOT GEOS_BUILD_OBJ_LIBS )
  message( "Building shared geosx_core library with static coreComponents, executables link to geosx_core" )
  blt_combine_static_libraries( NAME          geosx_core
                                SOURCE_LIBS   ${coreLibs}
                                LIB_TYPE      SHARED
                                LINK_PREPEND  ${GEOS_LINK_PREPEND_FLAG}
                                LINK_POSTPEND ${GEOS_LINK_POSTPEND_FLAG} )

elseif( NOT GEOS_BUILD_SHARED_LIBS AND GEOS_BUILD_OBJ_LIBS )
  message( "Building object coreComponents, executables link to coreComponents" )
else()
  message( "Building static coreComponents, executables link to coreComponents" )
endif()

if ( TARGET geosx_core )
  # Seems to be required on some CMake versions (e.g. 3.16) to get enforce device linking
  if( ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA" )
    set_target_properties( geosx_core PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS TRUE )
  endif()

  # To install the library with the runtime path used during the building
  set_target_properties( geosx_core PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )

  # To install the shared library
  install( TARGETS geosx_core LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

  target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
endif()

geosx_add_code_checks( PREFIX coreComponents
                       EXCLUDES cmake constitutive/PVTPackage )

if( ENABLE_UNCRUSTIFY )
  add_test( NAME testUncrustifyCheck
            COMMAND sh -c "${CMAKE_MAKE_PROGRAM} coreComponents_uncrustify_check main_uncrustify_check 2>&1 >/dev/null"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
endif()

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

