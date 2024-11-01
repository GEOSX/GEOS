#
# Specify all headers
#
set( math_headers
     interpolation/Interpolation.hpp
     extrapolation/Extrapolation.hpp)

blt_add_library( NAME         math
                 SOURCES               
                 HEADERS      ${math_headers}
                 DEPENDS_ON            
                 SHARED FALSE
               )
               
target_include_directories( math INTERFACE ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS math LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
