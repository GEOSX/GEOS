# -*- coding: utf-8 -*-

# Don't modify this file.
# Make a copy of this file renamed with your LC username.  Don't add your version to the git repo.
import platform
lassen = 'lassen' in platform.node()

#Should the particleStressHistory be in a sbatch wrapper so if the os hangs it isn't using job time
  
if lassen:
  geosPath='/usr/WS1/homel1/GEOS/build-lassen-gcc@8.3.1-relwithdebinfo/bin/geosx'
  testRunDirectory='/p/gpfs1/homel1/geosxRuns/test/'
  defaultRunDirectory='/p/gpfs1/homel1/geosxRuns/test/'
  defaultBank='cbronze'
else:
  geosPath='/usr/workspace/homel1/GEOS/build-quartz-gcc@12-release/bin/geosx'
  testRunDirectory='/p/lustre1/homel1/geosxRuns/test/'
  defaultRunDirectory='/p/lustre1/homel1/geosxRuns/test/'
  defaultBank='imcomp'