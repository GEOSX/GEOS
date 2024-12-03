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
Package: mainInterface

Contains high-level unit tests that depends typically depends on mainInterface components like the
ProblemManager.
#]]

add_subdirectory( toolchain )
add_subdirectory( testingUtilities )
add_subdirectory( xmlTests )
add_subdirectory( virtualElementTests )
add_subdirectory( linearAlgebraTests )
add_subdirectory( dataRepositoryTests )
add_subdirectory( fieldSpecificationTests )
add_subdirectory( constitutiveTests )
add_subdirectory( meshTests )
add_subdirectory( finiteVolumeTests )
add_subdirectory( fileIOTests )
if( GEOS_ENABLE_FLUIDFLOW )
  add_subdirectory( fluidFlowTests )
endif()
add_subdirectory( wellsTests )
add_subdirectory( wavePropagationTests ) 
