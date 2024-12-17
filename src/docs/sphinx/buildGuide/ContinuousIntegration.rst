.. _Continuous_Integration_process:

Continuous Integration process
==============================

To save building time, the third party libraries (that do not change so often) and GEOS are build separately.

Everytime a pull is requested in the TPL repository, docker images are generated and deployed on `dockerhub <https://hub.docker.com/r/geosx>`_.
The repository names (`ubuntu18.04-gcc8 <https://hub.docker.com/r/geosx/ubuntu18.04-gcc8>`_,
`centos7.7.1908-clang9.0.0 <https://hub.docker.com/r/geosx/centos7.5.1804-clang6.0.1>`_, `centos7.6.1810-gcc8.3.1-cuda10.1.243 <https://hub.docker.com/r/geosx/centos7.6.1810-gcc8.3.1-cuda10.1.243>`_ etc.)
obviously reflect the OS and the compiler flavour used.
For each image, the unique tag ``${PULL_REQUEST_NUMBER}-${BUILD_NUMBER}`` (defined as ``${{ github.event.number }}-${{ github.run_number }}`` in github actions) is used so we can connect the related code source in a rather convenient way.
Each docker contains the ``org.opencontainers.image.created`` and ``org.opencontainers.image.revision`` labels to provide additional information.

It is necessary to set ``build.args.GEOS_TPL_TAG`` (`e.g.` something like ``235-52``) in the `.devcontainer/devcontainer.json <https://github.com/GEOS-DEV/GEOS/blob/develop/.devcontainer/devcontainer.json>`_ file, to build against one selected version of the TPL.

It must be mentioned that one and only one version of the compiled TPL tarball is stored per pull request (older ones are removed automatically).
Therefore, a client building against a work in progress PR may experience a 404 error sooner or later.

Building docker images
----------------------

Our continuous integration process builds the TPL and GEOS against two operating systems (ubuntu and centos) and two compilers (clang and gcc).
The docker files use `multi-stage builds <https://docs.docker.com/develop/develop-images/multistage-build/>`_ in order to minimise the sizes of the images.

* First stage installs and defines all the elements that are commons to both TPL and GEOS (for example, MPI and c++ compiler, BLAS, LAPACK, path to the installation directory...).
* As a second stage, we install everything needed to build (`not run`) the TPLs.
  We keep nothing from this second step for GEOS, except the compiled TPL themselves.
  For example, a fortran compiler is needed by the TPL but not by GEOS: it shall be installed during this step, so GEOS won't access a fortran compiler (it does not have to).
* Last stage copies the compiled TPL from second stage and installs the elements only required by GEOS (there are few).

.. _Docker_images_contract:

Docker images contract
----------------------

GEOS will find a compiled version of the third party libraries.

As part of the contract provided by the TPL, the docker images also defines several environment variables.
The

.. code-block:: sh

    GEOS_TPL_DIR

variable contains the absolute path of the installation root directory of the third party libraries.
GEOS must use it when building.

Other variables are classical absolute path compiler variables.

.. code-block:: sh

    CC
    CXX
    MPICC
    MPICXX

And the absolute path the mpirun (or equivalent) command.

.. code-block:: sh

    MPIEXEC

The following ``openmpi`` environment variables allow it to work properly in the docker container.
But there should be no reason to access or use them explicitly.

.. code-block:: sh

    OMPI_CC=$CC
    OMPI_CXX=$CXX
