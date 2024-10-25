.. _EventManager:

Event Management
===============================================================================

The goal of the GEOS event manager is to be flexible with regards to event type, application order, and method of triggering.  The event manager is configured via the ``Event`` block in an input .xml file, i.e.:

.. code-block:: xml

  <Events maxTime="1.0e-2">
    <PeriodicEvent name="event_a"
                   target="/path/to/event"
                   forceDt="1" />
    <HaltEvent name="event_b"
               target="/path/to/halt_target"
               maxRunTime="1e6" />
  </Events>


Event Execution Rules
---------------------------------------------

The EventManager will repeatedly iterate through a list of candidate events specified via the Events block **in the order they are defined in the xml**.  When certain user-defined criteria are met, they will trigger and perform a task.  The simulation ``cycle`` denotes the number of times the primary event loop has completed, ``time`` denotes the simulation time at the beginning of the loop, and ``dt`` denotes the global timestep during the loop.

During each cycle, the EventManager will do the following:

1. Loop through each event and obtain its timestep request by considering:

   a. The maximum dt specified via the target's GetTimestepRequest method
   b. The time remaining until user-defined points (e.g. application start/stop times)
   c. Any timestep overrides (e.g. user-defined maximum dt)
   d. The timestep request for any of its children

2. Set the cycle dt to the smallest value requested by any event

3. Loop through each event and:

   a. Calculate the event ``forecast``, which is defined as the expected number of cycles until the event is expected to execute.
   b. ``if (forecast == 1)`` the event will signal its target to prepare to execute.  This is useful for preparing time-consuming I/O operations.
   c. ``if (forecast <= 0)`` the event will call the Execute method on its target object

4. Check to see if the EventManager exit criteria have been met


After exiting the main event loop, the EventManager will call the ``Cleanup`` method for each of its children (to produce final plots, etc.).  Note: if the code is resuming from a restart file, the EventManager will pick up exactly where it left off in the execution loop.


Event Manager Configuration
----------------------------

Event
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The children of the Event block define the events that may execute during a simulation.  These may be of type ``HaltEvent``, ``PeriodicEvent``, or ``SoloEvent``.  The exit criteria for the global event loop are defined by the attributes ``maxTime`` and ``maxCycle`` (which by default are set to their max values).  If the optional logLevel flag is set, the EventManager will report additional information with regards to timestep requests and event forecasts for its children.

.. include:: /docs/sphinx/datastructure/Events.rst
    :start-line: 3

PeriodicEvent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is the most common type of event used in GEOS.  As its name suggests, it will execute periodically during a simulation.  It can be triggered based upon a user-defined ``cycleFrequency`` or ``timeFrequency``.

If cycleFrequency is specified, the event will attempt to execute every X cycles.  Note: the default behavior for a PeriodicEvent is to execute every cycle.  The event forecast for this case is given by: ``forecast = cycleFrequency - (cycle - lastCycle)`` .

If timeFrequency is specified, the event will attempt to execute every X seconds (this will override any cycle-dependent behavior).  By default, the event will attempt to modify its timestep requests to respect the timeFrequency (this can be turned off by specifying targetExactTimestep="0").  The event forecast for this case is given by: ``if (dt > 0), forecast = (timeFrequency - (time - lastTime)) / dt, otherwise forecast=max``

By default, a PeriodicEvent will execute throughout the entire simulation.  This can be restricted by specifying the beginTime and/or endTime attributes.  Note: if either of these values are set, then the event will modify its timestep requests so that a cycle will occur at these times (this can be turned off by specifying targetExactStartStop="0").

The timestep request event is typically determined via its target.  However, this value can be overridden by setting the ``forceDt`` or ``maxEventDt`` attributes.

.. include:: /docs/sphinx/datastructure/PeriodicEvent.rst
    :start-line: 3

SoloEvent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This type of event will execute once once the event loop reaches a certain cycle (targetCycle) or time (targetTime).  Similar to the PeriodicEvent type, this event will modify its timestep requests so that a cycle occurs at the exact time requested (this can be turned off by specifying targetExactTimestep="0").  The forecast calculations follow an similar approach to the PeriodicEvent type.

.. include:: /docs/sphinx/datastructure/SoloEvent.rst
    :start-line: 3

HaltEvent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This event type is designed to track the wall clock.  When the time exceeds the value specified via maxRunTime, the event will trigger and set a flag that instructs the main EventManager loop to cleanly exit at the end of the current cycle.  The event for cast for this event type is given by: ``forecast = (maxRuntime - (currentTime - startTime)) / realDt``

.. include:: /docs/sphinx/datastructure/HaltEvent.rst
    :start-line: 3



Other Event Features
---------------------------------------------

Event Progress Indicator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Because the event manager allows the user to specify the order of events, it could introduce ambiguity into the timestamps of output files.  To resolve this, we pass two arguments to the target's Execute method:

1. eventCounter (integer) - the application index for the event (or sub-event)
2. eventProgress (real64) - the percent completion of the event loop, paying attention to events whose targets are associated with physics (from the start of the event, indicated via target->GetTimestepBehavior())

For example, consider the following Events block:

.. code-block:: xml

  <Events maxTime="1.0e-2">
    <PeriodicEvent name="outputs"
                   timeFrequency="1e-6"
                   targetExactTimestep="0"
                   target="/Outputs/siloOutput">
    <PeriodicEvent name="solverApplications_a"
                   forceDt="1.0e-5"
                   target="/Solvers/lagsolve" />
    <PeriodicEvent name="solverApplications_b"
                   target="/Solvers/otherSolver" />
    <PeriodicEvent name="restarts"
                   timeFrequency="5.0e-4"
                   targetExactTimestep="0"
                   target="/Outputs/restartOutput"/>
  </Events>

In this case, the events solverApplications_a and solverApplications_b point target physics events.  The eventCounter, eventProgress pairs will be: outputs (0, 0.0), solverApplications_a (1, 0.0), solverApplications_b (2, 0.5), and restarts (3, 1.0).  These values are supplied to the target events via their Execute methods for use.  For example, for the name of a silo output file will have the format: "%s_%06d%02d" % (name, cycle, eventCounter), and the time listed in the file will be ``time = time + dt*eventProgress``



Nested Events
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The event manager allows its child events to be nested.  If this feature is used, then the manager follows the basic execution rules, with the following exception:  When its criteria are met, an event will first execute its (optional) target.  It will then estimate the forecast for its own sub-events, and execute them following the same rules as in the main loop.  For example:

.. code-block:: xml

  <Events maxTime="1.0e-2">
    <PeriodicEvent name="event_a"
                   target="/path/to/target_a" />

    <PeriodicEvent name="event_b"
                   timeFrequency="100">

      <PeriodicEvent name="subevent_b_1"
                     target="/path/to/target_b_1"/>

      <PeriodicEvent name="subevent_b_2"
                     target="/path/to/target_b_2"/>
    <PeriodicEvent/>
  </Events>

In this example, event_a will trigger during every cycle and call the Execute method on the object located at /path/to/target_a.  Because it is time-driven, event_b will execute every 100 s.  When this occurs, it will execute it will execute its own target (if it were defined), and then execute subevent_b_1 and subevent_b_2 in order. Note: these are both cycle-driven events which, by default would occur every cycle.  However, they will not execute until each of their parents, grandparents, etc. execution criteria are met as well.

