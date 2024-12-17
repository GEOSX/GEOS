:orphan:

============
Outputs
============

This section describes how outputs are handled by GEOS

The outputs are defined in a ``<Outputs>`` XML block.

There are three available formats to output the results of a simulation: SILO_, VTK_, and Time History output into simple dataset HDF5_ files which are consumable by post-processing scripts..

************************
Defining an output
************************

SILO Output
===========

The SILO output is defined through the ``<Silo>`` XML node (subnode of ``<Outputs> XML block``) as shown here:

.. code-block:: xml 

  <Outputs>
    <Silo name="siloOutput"/>
  </Outputs>
  
  
The parameter options are listed in the following table:

.. include:: /docs/sphinx/datastructure/Silo.rst


VTK Output
===========

The VTK output is defined through the ``<VTK>`` XML node (subnode of ``<Outputs> XML block``) as shown here:

.. code-block:: xml 

  <Outputs>
    <VTK name="vtkOutput"/>
  </Outputs>

The parameter options are listed in the following table:

.. include:: /docs/sphinx/datastructure/VTK.rst

TimeHistory Output
==================

The TimeHistory output is defined through the ``<TimeHistory>`` XML node (subnode of ``<Outputs> XML block``) as shown here:

.. code-block:: xml

  <Outputs>
    <TimeHistory name="timeHistoryOutput" sources="{/Tasks/collectionTask}" filename="timeHistory" />
  </Outputs>

The parameter options are listed in the following table:

.. include:: /docs/sphinx/datastructure/TimeHistory.rst

In order to properly collect and output time history information the following steps must be accomplished:

#. Specify one or more collection tasks using the :ref:`TasksManager`.
#. Specify a `TimeHistory Output` using the collection task(s) as source(s).
#. Specify an event in the :ref:`EventManager` to trigger the collection task(s).
#. Specify an event in the :ref:`EventManager` to trigger the output.

Note: Currently if the collection and output events are triggered at the same simulation time, the one specified first will also trigger first. Thus in order to output time history for the current time in this case, always specify the time history collection events prior to the time history output events.

************************
Triggering the outputs
************************

The outputs can be triggered using the :ref:`EventManager`.
It is recommended to use a ``<PeriodicEvent>`` to output results with a defined frequency:

.. code-block:: xml

  <PeriodicEvent name="outputs"
                 timeFrequency="5000.0"
                 targetExactTimestep="1"
                 target="/Outputs/siloOutput" />

The keyword ``target`` has to match with the name of the ``<Silo>``, ``<VTK>``, or ``<TimeHistory>`` node.

****************************
Visualisation of the outputs
****************************

We suggest the use of VisIT_, Paraview_, and MatPlotLib_ to visualize the outputs.

Visualizing Silo outputs with VisIT
===================================

If the ``<Silo>`` XML node was defined, GEOS writes the results in a folder called ``siloFiles``.

In VisIT :

1. File > Open file...
2. On the right panel, browse to the ``siloFiles`` folder.
3. On the left panel, select the file(s) you want to visualize. Usually, one file is written according the
   frequency defined in the ``timeFrequency`` keyword of the Event that has triggered the output.
4. To load fields, use the "Add" button and browse to the fields you want to plot.
5. To plot fields, use the "Draw" button.

Please consult the VisIT_ documentation for further explanations on its usage.


Visualizing VTK outputs with VisIT
===================================

If the ``<VTK>`` XML node was defined, GEOS writes the results in a folder named after the ``plotFileRoot`` attribute (default = ``vtkOutputs``).
For problems with multiple active regions (e.g. :ref:`TutorialHydraulicFractureWithAdvancedXML`), additional work may be required to ensure that vtk files can be read by VisIt.
Options include:

1. Using a VisIt macro / python script to convert default multi-region vtk files (see ``GEOS/src/coreComponents/python/visitMacros/visitVTKConversion.py``).
2. Using the ``outputRegionType`` attribute to output separate sets of files per region.  For example:

.. code-block:: xml

  <Problem>
    <Events>
      <!-- Use nested events to trigger both vtk outputs -->
      <PeriodicEvent
        name="outputs"
        timeFrequency="2.0"
        targetExactTimestep="0">
        <PeriodicEvent
          name="outputs_cell"
          target="/Outputs/vtkOutput_cell"/>
        <PeriodicEvent
          name="outputs_surface"
          target="/Outputs/vtkOutput_surface"/>
      </PeriodicEvent>
    </Events>
    
    <Outputs>
      <VTK
        name="vtkOutput_cell"
        outputRegionType="cell"
        plotFileRoot="vtk_cell"/>

      <VTK
        name="vtkOutput_surface"
        outputRegionType="surface"
        plotFileRoot="vtk_surface"/>
    </Outputs>
  </Problem>


In VisIt, the prepared files can be visualized as follows:

1. Click the Open icon (or select File/Open from the top menu).
2. Enter the vtk output folder in the left-hand panel (these will be separate if you use option 2 above).
3. Select the desired ``*.vtm`` group in the right-hand panel (by double-clicking or selecting OK).
4. (Repeat for any additional datasets/regions.)
5. In the main window, select the desired dataset from the Active Sources dropdown menu.
6. Click on the Add icon, and select a parameter to visualize from the mesh, psuedocolor, or vector plot options
7. At some point you may be prompted to create a Database Correlation.  If you select Yes, then VisIt will create a new time slider, which can be used to change the time state while keeping the various aspects of the simulation synchronized.
8. Finally, click the Draw button and use the time slider to visualize the output.

Please consult the VisIT_ documentation for further explanations on its usage.


Visualizing VTK outputs with Paraview
=====================================

If the ``<VTK>`` XML node was defined, GEOS writes a folder and a ``.pvd`` file named after the string defined
in ``name`` keyword.

The ``.pvd`` file contains references to the ``.pvtu`` files. One ``.pvtu`` file is output according the frequency defined in the ``timeFrequency`` keyword of the Event that has triggered the output.

One ``.pvtu`` contains references to ``.vtu`` files. There is as much ``.vtu`` file as there were MPI processes
used for the computation.

All these files can be opened with paraview. To have the whole results for every output time steps, you can
open the ``.pvd`` file.

Visualizing TimeHistory outputs with MatPlotLib
===============================================

If the ``<TimeHistory>`` XML node was defined, GEOS writes a file named after the string defined
in the ``filename`` keyword and formatted as specified by the string defined in the ``format``
keyword (only HDF5_ is currently supported).

The TimeHistory file contains the collected time history information from each specified time history collector.
This information includes datasets for the time itself, any metadata sets describing index association with specified
collection sets, and the time history information itself.

It is recommended to use MatPlotLib_ and format-specific accessors (like H5PY for HDF5_) to access and easily plot the
time history datat.

.. _SILO: https://wci.llnl.gov/simulation/computer-codes/silo
.. _VTK: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
.. _HDF5: https://portal.hdfgroup.org/display/HDF5/HDF5
.. _VisIT: https://wci.llnl.gov/simulation/computer-codes/visit/downloads
.. _Paraview: https://www.paraview.org/
.. _MatPlotLib: https://matplotlib.org/
