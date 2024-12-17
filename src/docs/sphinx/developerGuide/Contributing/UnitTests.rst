**************************************
Unit Testing
**************************************

Unit testing is integral to the GEOS development process. While not all components naturally lend themselves to unit testing (for example a physics solver) every effort should be made to write comprehensive quality unit tests.

Each sub-directory in ``coreComponents`` should have a ``unitTests`` directory containing the test sources. Each test consists of a ``cpp`` file whose name begins with ``test`` followed by a name to describe the test. Please read over the `LvArray unit test documentation <https://lvarray.readthedocs.io/en/latest/testing.html>`_ as it gives an intro to the Google Test framework and a set of best practices.

GEOS Specific Recommendations
------------------------------

An informative example is ``testSinglePhaseBaseKernels`` which tests the single phase flow mobility and accumulation kernels on a variety of inputs.

.. literalinclude:: ../../../../coreComponents/unitTests/fluidFlowTests/testSinglePhaseMobilityKernel.cpp
   :language: c++
   :start-after: // Sphinx start after test mobility
   :end-before: // Sphinx end before test mobility

*[Source: coreComponents/physicsSolvers/fluidFlow/unitTests/testSinglePhaseMobilityKernel.cpp]*

What makes this such a good test is that it depends on very little other than kernels themselves. There is no need to involve the data repository or parse an XML file. Sometimes however this is not possible, or at least not without a significant duplication of code. In this case it is better to embed the XML file into the test source as a string instead of creating a separate XML file and passing it to the test as a command line argument or hard coding the path. One example of this is ``testLaplaceFEM`` which tests the laplacian solver. The embedded XML is shown below.

.. literalinclude:: ../../../../coreComponents/unitTests/fluidFlowTests/testCompMultiphaseFlow.cpp
   :language: c++
   :start-after: // Sphinx start after input XML
   :end-before: // Sphinx end before input XML

*[Source: coreComponents/physicsSolvers/fluidFlow/unitTests/testCompMultiphaseFlow.cpp]*

MPI
---
Often times it makes sense to write a unit test that is meant to be run with multiple MPI ranks. This can be accomplished by simply adding the ``NUM_MPI_TASKS`` parameter to ``geos_add_test`` in the CMake file. For example

::

  geos_add_test( NAME testWithMPI
                 COMMAND testWithMPI
                 NUM_MPI_TASKS ${NUMBER_OF_MPI_TASKS} )

With this addition ``make test`` or calling ``ctest`` directly will run ``testWithMPI`` via something analogous to ``mpirun -n NUMBER_OF_MPI_TASKS testWithMPI``.
