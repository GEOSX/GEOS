[![DOI](https://zenodo.org/badge/131810578.svg)](https://zenodo.org/badge/latestdoi/131810578)
[![codecov](https://codecov.io/github/GEOS-DEV/GEOS/graph/badge.svg?token=0VTEHPQG58)](https://codecov.io/github/GEOS-DEV/GEOS)
![CI](https://github.com/GEOS-DEV/GEOS/actions/workflows/ci_tests.yml/badge.svg)
[![docs](https://readthedocs.com/projects/geosx-geosx/badge/?version=latest)](https://geosx-geosx.readthedocs-hosted.com/en/latest/)

Welcome to the GEOS project!
-----------------------------
GEOS is a simulation framework for modeling coupled flow, transport, and geomechanics
in the subsurface.  The code provides advanced solvers for a number of target applications,
including
  - carbon sequestration,
  - geothermal energy,
  - and similar systems.  

A key focus of the project is achieving scalable performance on current and next-generation
high performance computing systems.  We do this through a portable programming model and research into scalable algorithms.

You may want to browse our
[publications](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Publications.html)
page for more details on the HPC, numerics,
and applied engineering components of this effort.

Documentation
---------------------

Our documentation is hosted [here](https://geosx-geosx.readthedocs-hosted.com/en/latest/?).


Who develops GEOS?
-------------------
GEOS is an open source project and is developed by a community of researchers at
several institutions.  The bulk of the code has been written by contributors from
four main organizations:
  - Lawrence Livermore National Laboratory,
  - Stanford University,
  - TotalEnergies,
  - Chevron

See our
[authors](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Contributors.html)
and
[acknowledgements](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Acknowledgements.html)
page for more details.  

How does GEOS relate to the earlier GEOS code?
------------------------------
GEOS is the offshoot of an earlier code developed at LLNL also called GEOS.  The new
code differs from our previous efforts in two important ways:
  - This new code GEOS uses a fundamentally different programming model to achieve
    high performance on the complicated chip architectures common on today's
    HPC systems.  This code is ready for exascale-class systems as they are delivered.
  - The new code has been released as an open-source effort to encourage collaboration
    within the research and industrial community.  See the release notes below
    for details of the [LGPL 2.1 License](./LICENSE) that has been adopted.


Release
-------

For release details and restrictions, please read the [LICENSE](./LICENSE) file.

For copyrights, please read the [COPYRIGHT](./COPYRIGHT ) file.

For contributors, please read the [CONTRIBUTORS](./CONTRIBUTORS ) file.

For acknowledgements, please read the [ACKNOWLEDGEMENTS](./ACKNOWLEDGEMENTS ) file.

For notice, please read the [NOTICE](./NOTICE ) file.

`LLNL-CODE-812638`  `OCEC-18-021`
