.. _Overview:

Packages Overview
##################################################

.. contents:: Table of contents
    :depth: 2


This page presents an overview of the main code components of GEOS for a developer. You will find a dependency scheme representing all the relationship between packages and brief documentation of the core packages with links to their main components documentation.


.. figure:: dependency.png
    :align: center
    :figclass: align-center

    GEOS packages dependencies scheme.

**TODO: change figure path + automatic generation on make docs**


codingUtilities
^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../../../coreComponents/codingUtilities/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


common
^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../../../coreComponents/common/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


constitutive
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/constitutive/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]

More information:
    - :ref:`constitutiveModelsDoc`
    - **ref/link to PVT package** see `PVTPAckage <../coreComponents/constitutive/PVTPackage/docs/main.html>`_


constitutiveDrivers
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/constitutiveDrivers/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]

More information:
    - :ref:`PVTDriver`   **example of PVT driver and how to add unit test**


dataRepository
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/dataRepository/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]

More information:
    - :ref:`dataRepository`


denseLinearAlgebra
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/denseLinearAlgebra/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


discretizationMethods
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/discretizationMethods/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


events
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/events/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


fieldSpecification
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/fieldSpecification/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


fileIO
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/fileIO/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


finiteElement
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/finiteElement/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


More information:
    - :ref:`kernelInterface`


finiteVolume
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/finiteVolume/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


functions
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/functions/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


linearAlgebra
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/linearAlgebra/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]

More information:
    - :ref:`DoFManager`


mainInterface
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/mainInterface/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


math
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/math/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]



mesh
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/mesh/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]


More information:
    - :ref:`Meshes`  **Page from user guide, include it here ?**
    - :ref:`MeshHierarchy`


physicsSolvers
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/physicsSolvers/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]

More information:
    - :ref:`AddingNewSolver`

schema
^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../../coreComponents/schema/CMakeLists.txt
    :start-after: #[[
    :end-before: #]]

More information:
    - :ref:`InputSchemaGenerationPar`