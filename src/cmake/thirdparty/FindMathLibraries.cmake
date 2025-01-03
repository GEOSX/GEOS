
if( ENABLE_MKL )
    if ( NOT DEFINED MKL_INCLUDE_DIRS )
        message( FATAL_ERROR "MKL is enabled but MKL_INCLUDE_DIRS is not defined." )
    elseif ( NOT DEFINED MKL_LIBRARIES )
        message( FATAL_ERROR "MKL is enabled but MKL_LIBRARIES is not defined." )
    endif()

    message( STATUS "MKL found" )

    set( BLAS_LIBRARIES "${MKL_LIBRARIES}" CACHE STRING "" FORCE )
    set( LAPACK_LIBRARIES "${MKL_LIBRARIES}" CACHE STRING "" FORCE)

elseif( ENABLE_ESSL )
    if ( NOT DEFINED ESSL_INCLUDE_DIRS )
        message( FATAL_ERROR "ESSL is enabled but ESSL_INCLUDE_DIRS is not defined." )
    elseif ( NOT DEFINED ESSL_LIBRARIES )
        message( FATAL_ERROR "ESSL is enabled but ESSL_LIBRARIES is not defined." )
    endif()

    set( BLAS_LIBRARIES "${ESSL_LIBRARIES}" CACHE STRING "" FORCE )
    set( LAPACK_LIBRARIES "${ESSL_LIBRARIES}" CACHE STRING "" FORCE )

endif()


if( NOT DEFINED BLAS_LIBRARIES )
    find_package( BLAS REQUIRED )
endif()

if( NOT DEFINED LAPACK_LIBRARIES )
    find_package( LAPACK REQUIRED )
endif()

# # Create the BLAS link line (used by hypre)
# set( BLAS_SHARED_LIBRARY_DIRS "" )
# foreach( lib ${BLAS_LIBRARIES} )
#     if ( NOT EXISTS ${lib} )
#         message( ERROR "BLAS library '${lib}' does not exist!" )
#     endif()

#     if ( NOT lib MATCHES "\.a$" )
#         get_filename_component( dirname ${lib} DIRECTORY )
#         list( APPEND BLAS_SHARED_LIBRARY_DIRS ${dirname} )
#     endif()
# endforeach()

# list( REMOVE_DUPLICATES BLAS_SHARED_LIBRARY_DIRS )
# string( REPLACE ";" ":" BLAS_SHARED_LIBRARY_DIRS "${BLAS_SHARED_LIBRARY_DIRS}" )

# string( REPLACE ";" " " BLAS_LIBRARIES_STRING "${BLAS_LIBRARIES}" )
# set( BLAS_LINK_LINE "-Wl,-rpath,${BLAS_SHARED_LIBRARY_DIRS} ${BLAS_LIBRARIES_STRING}" CACHE STRING "" )

message( STATUS "BLAS_LIBRARIES = ${BLAS_LIBRARIES}" )
# message( STATUS "BLAS_LINK_LINE = ${BLAS_LINK_LINE}" )

# # Create the LAPACK link line (used by hypre)
# set( LAPACK_SHARED_LIBRARY_DIRS "" )
# foreach( lib ${LAPACK_LIBRARIES} )
#     if ( NOT EXISTS ${lib} )
#         message( ERROR "LAPACK library '${lib}' does not exist!" )
#     endif()

#     if ( NOT lib MATCHES "\.a$" )
#         get_filename_component( dirname ${lib} DIRECTORY )
#         list( APPEND LAPACK_SHARED_LIBRARY_DIRS ${dirname} )
#     endif()
# endforeach()

# list( REMOVE_DUPLICATES LAPACK_SHARED_LIBRARY_DIRS )
# string (REPLACE ";" ":" LAPACK_SHARED_LIBRARY_DIRS "${LAPACK_SHARED_LIBRARY_DIRS}" )

# string( REPLACE ";" " " LAPACK_LIBRARIES_STRING "${LAPACK_LIBRARIES}" )
# set ( LAPACK_LINK_LINE "-Wl,-rpath,${LAPACK_SHARED_LIBRARY_DIRS} ${LAPACK_LIBRARIES_STRING}" CACHE STRING "" )

message( STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}" )
# message( STATUS "LAPACK_LINK_LINE = ${LAPACK_LINK_LINE}" )
