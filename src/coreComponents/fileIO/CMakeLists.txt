# Specify all headers
set( fileIO_headers
     Outputs/BlueprintOutput.hpp
     Outputs/OutputBase.hpp
     Outputs/OutputManager.hpp
     Outputs/OutputUtilities.hpp
     Outputs/PythonOutput.hpp
     Outputs/RestartOutput.hpp
     Table/TableLayout.hpp
     Table/TableFormatter.hpp
     Table/TableData.hpp
     Outputs/TimeHistoryOutput.hpp
     timeHistory/HDFFile.hpp
     timeHistory/HistoryCollectionBase.hpp
     timeHistory/BufferedHistoryIO.hpp
     timeHistory/PackCollection.hpp
     timeHistory/HDFHistoryIO.hpp
     timeHistory/HistoryCollection.hpp )

# Specify all sources
set( fileIO_sources
     Outputs/BlueprintOutput.cpp
     Outputs/OutputBase.cpp
     Outputs/OutputManager.cpp
     Outputs/OutputUtilities.cpp
     Outputs/PythonOutput.cpp
     Outputs/RestartOutput.cpp
     Table/TableLayout.cpp
     Table/TableFormatter.cpp
     Table/TableData.cpp
     Outputs/TimeHistoryOutput.cpp
     timeHistory/HDFFile.cpp
     timeHistory/HistoryCollectionBase.cpp
     timeHistory/PackCollection.cpp
     timeHistory/HDFHistoryIO.cpp )

set( dependencyList ${parallelDeps} mesh constitutive hdf5 )
if( ENABLE_PYGEOSX )
  list( APPEND fileIO_headers
        python/PyHistoryCollectionType.hpp
        python/PyHistoryOutputType.hpp
        python/PyVTKOutputType.hpp )
  list( APPEND fileIO_sources
        python/PyHistoryCollection.cpp
        python/PyHistoryOutput.cpp
        python/PyVTKOutput.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()


if( ENABLE_MPI )
  add_subdirectory( coupling/hdf5_interface )
  list( APPEND dependencyList mpi hdf5_interface )
  list( APPEND fileIO_headers coupling/ChomboCoupler.hpp Outputs/ChomboIO.hpp )
  list( APPEND fileIO_sources coupling/ChomboCoupler.cpp Outputs/ChomboIO.cpp )
endif()

if( ENABLE_SILO )
  list( APPEND dependencyList silo )
  list( APPEND fileIO_headers 
               silo/SiloFile.hpp 
               Outputs/SiloOutput.hpp )
  list( APPEND fileIO_sources 
               silo/SiloFile.cpp
               Outputs/SiloOutput.cpp )
endif( )

if( ENABLE_VTK )
  list( APPEND fileIO_headers
        vtk/VTKPVDWriter.hpp
        vtk/VTKVTMWriter.hpp
        vtk/VTKPolyDataWriterInterface.hpp
        Outputs/VTKOutput.hpp
        )
  list( APPEND fileIO_sources
        vtk/VTKPVDWriter.cpp
        vtk/VTKVTMWriter.cpp
        vtk/VTKPolyDataWriterInterface.cpp
        Outputs/VTKOutput.cpp
        )
    list( APPEND dependencyList VTK::IOLegacy VTK::IOXML )
endif()

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND dependencyList CUDA::nvToolsExt )
endif()

blt_add_library( NAME       fileIO
                 SOURCES    ${fileIO_sources}
                 HEADERS    ${fileIO_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS} )


target_include_directories( fileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( Table/unitTests )
endif()



