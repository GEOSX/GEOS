# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: dataRepository

Provides the building blocks for the data structure of GEOS objects.
Also contains a wrapper to process entries from an xml file into data types.
#]]

#
# Specify all headers
#
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

#
# Specify all sources
#
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
