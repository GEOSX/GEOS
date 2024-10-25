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
     LogLevelsInfo.hpp
     LogLevelsRegistry.hpp
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
     LogLevelsRegistry.cpp
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

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       dataRepository
                 SOURCES    ${dataRepository_sources}
                 HEADERS    ${dataRepository_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
install( TARGETS dataRepository LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )


if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )
