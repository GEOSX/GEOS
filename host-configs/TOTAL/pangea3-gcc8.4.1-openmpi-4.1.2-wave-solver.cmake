# hostconfig to build only the wave solver on pangea3
#
#
include( ./pangea3-gcc8.4.1-openmpi-4.1.2.cmake )

set ( GEOS_ENABLE_CONTACT OFF )
set ( GEOS_ENABLE_FLUIDFLOW OFF )
set ( GEOS_ENABLE_INDUCEDSEISMICITY OFF )
set ( GEOS_ENABLE_MULTIPHYSICS OFF )
set ( GEOS_ENABLE_SIMPLEPDE OFF )
set ( GEOS_ENABLE_SOLIDMECHANICS OFF )
set ( GEOS_ENABLE_SURFACEGENERATION OFF )

