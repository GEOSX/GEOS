.. _TasksManager:

Tasks Manager
=============

The GEOS tasks manager allows a user to specify tasks to be executed. These tasks are compatible targets for the :ref:`EventManager`.

The tasks manager is configured via the ``Tasks`` block in an input .xml file, i.e.:

.. code-block:: xml

   <Tasks>
     <PackCollection name="historyCollection" objectPath="nodeManager" fieldName="Velocity" />
   </Tasks>


Tasks Manager Configuration
---------------------------

Task
***************************
The children of the Tasks block define different Tasks to be triggered by events specified in the :ref:`EventManager` during the execution of the simulation. At present the only supported task is the ``PackCollection`` used to collect time history data for output by a TimeHistory output.

.. include:: /docs/sphinx/datastructure/Tasks.rst

PackCollection
***************************
The ``PackCollection`` Task is used to collect time history information from fields. Either the entire field or specified named sets of indices in the field can be collected.

.. include:: /docs/sphinx/datastructure/PackCollection.rst

Note: The time history information collected via this task is buffered internally until it is output by a linked TimeHistory Output.


***************************
Triggering the Tasks
***************************
Tasks can be triggered using the :ref:`EventManager`.
Recurring tasks sould use a ``<PeriodicEvent>`` and one-time tasks should use a `<SoloEvent>`:

.. code-block:: xml

  <PeriodicEvent name="historyCollectEvent"
                 timeFrequency="1.0"
                 targetExactTimeset="1"
                 target="/Tasks/historyCollection" />

The keyword ``target`` has to match the ``name`` of a Task specified as a child of the ``<Tasks>`` block.

