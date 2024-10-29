# Specify all headers
set( dataRepository_headers
     BufferOps.hpp
     BufferOpsDevice.hpp
     BufferOps_inline.hpp
     ConduitRestart.hpp
     DefaultValue.hpp
     ExecutableGroup.hpp
     Group.hpp
     HistoryDataSpec.hpp
     InputFlags.hpp
     KeyIndexT.hpp
     KeyNames.hpp
     MappedVector.hpp
     ObjectCatalog.hpp
     ReferenceWrapper.hpp
     RestartFlags.hpp
     Utilities.hpp
     Wrapper.hpp
     WrapperBase.hpp
     wrapperHelpers.hpp
     xmlWrapper.hpp
     DataContext.hpp
     GroupContext.hpp
     WrapperContext.hpp )

# Specify all sources
set( dataRepository_sources
     BufferOpsDevice.cpp
     ConduitRestart.cpp
     ExecutableGroup.cpp
     Group.cpp
     Utilities.cpp
     WrapperBase.cpp
     xmlWrapper.cpp
     DataContext.cpp
     GroupContext.cpp
     WrapperContext.cpp )

set( dependencyList ${parallelDeps} codingUtilities )

if( ENABLE_PYGEOSX )
  list( APPEND dataRepository_headers
        python/PyGroup.hpp
        python/PyGroupType.hpp
        python/PyWrapper.hpp )
  list( APPEND dataRepository_sources
        python/PyGroup.cpp
        python/PyWrapper.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()


blt_add_library( NAME       dataRepository
                 SOURCES    ${dataRepository_sources}
                 HEADERS    ${dataRepository_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS} )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )


if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )
