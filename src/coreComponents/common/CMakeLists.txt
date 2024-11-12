#
# Specify all headers
#
set( common_headers
     ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
     format/table/TableLayout.hpp
     format/table/TableFormatter.hpp
     format/table/TableData.hpp
     format/Format.hpp
     format/StringUtilities.hpp
     logger/Logger.hpp
     BufferAllocator.hpp
     DataLayouts.hpp
     DataTypes.hpp
     FieldSpecificationOps.hpp
     GEOS_RAJA_Interface.hpp
     GeosxMacros.hpp
     MemoryInfos.hpp
     MpiWrapper.hpp
     Path.hpp
     Span.hpp
     Stopwatch.hpp
     Timer.hpp
     Tensor.hpp
     TimingMacros.hpp
     TypeDispatch.hpp
     initializeEnvironment.hpp
     LifoStorage.hpp
     LifoStorageCommon.hpp
     LifoStorageHost.hpp
     FixedSizeDeque.hpp
     FixedSizeDequeWithMutexes.hpp
     MultiMutexesLock.hpp
     PhysicsConstants.hpp
     Units.hpp
   )

if ( ENABLE_CUDA )
  list( APPEND common_headers LifoStorageCuda.hpp )
endif( )

#
# Specify all sources
#
set( common_sources
     format/table/TableLayout.cpp
     format/table/TableFormatter.cpp
     format/table/TableData.cpp
     format/StringUtilities.cpp
     logger/Logger.cpp
     BufferAllocator.cpp
     MemoryInfos.cpp
     MpiWrapper.cpp
     Path.cpp
     initializeEnvironment.cpp
     Units.cpp
   )

set( dependencyList ${parallelDeps} lvarray RAJA chai )

if (TARGET conduit::conduit)
  set( dependencyList ${dependencyList} conduit::conduit )
endif()

if (TARGET conduit)
  set( dependencyList ${dependencyList} conduit )
endif()

if (TARGET fmt::fmt-header-only)
  set( dependencyList ${dependencyList} fmt::fmt-header-only )
endif()

if (TARGET fmt)
  set( dependencyList ${dependencyList} fmt )
endif()

if (TARGET pugixml::pugixml)
  set( dependencyList ${dependencyList} pugixml::pugixml )
endif()

if (TARGET pugixml)
  set( dependencyList ${dependencyList} pugixml )
endif()

if ( ENABLE_MKL )
  list( APPEND dependencyList mkl )
endif()

if( ENABLE_CALIPER )
  list( APPEND dependencyList caliper )
endif()

if( ENABLE_CALIPER )
  list( APPEND dependencyList caliper )

  if (TARGET adiak::adiak)
    list( APPEND dependencyList adiak::adiak )
  endif()

  if (TARGET adiak)
    list( APPEND dependencyList adiak )
  endif()
endif()

blt_add_library( NAME       common
                 SOURCES    ${common_sources}
                 HEADERS    ${common_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( common PUBLIC ${CMAKE_BINARY_DIR}/include )
target_include_directories( common PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS common LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
  add_subdirectory( format/table/unitTests )
  add_subdirectory( format/unitTests )
endif()
