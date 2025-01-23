.. _LogLevelDocumentation:

Log levels documentation
========================

Add a log level
---------------

To add a log level, you must respect the following structure and add it to the appropriate ``LogLevelsInfos.hpp`` :

.. code-block:: c++

    struct MyMessage
    {
        static constexpr int getMinLogLevel() { return 2; }
        static constexpr std::string_view getDescription() { return msg; }
    };

If there is no ``LogLevelsInfos.hpp`` in the corresponding folder, you can create a ``LogLevelsInfos.hpp``

Example of usage
----------------

To log a message with a log level, 4 macros are defined in LogLevelsInfo.hpp located in dataRepository:

* ``GEOS_LOG_LEVEL_INFO( logInfoStruct, msg )``: Output messages based on current ``Group``'s log level.
* ``GEOS_LOG_LEVEL_INFO_RANK_0( logInfoStruct, msg )``: Output messages (only on rank 0) based on current ``Group``'s log level.
* ``GEOS_LOG_LEVEL_INFO_BY_RANK( logInfoStruct, msg )``: Output messages (with one line per rank) based on current ``Group``'s log level.
* ``GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfoStruct, msg )``: Output messages (only on rank 0) based on current ``Group``'s log level without the line return.

An exemple of adding and using a log level in a ``group``:

.. code-block:: c++

    MyGroup::MyGroup( string const & name );
    {
        // To use a log level, make sure it is declared in the constructor of the group
        addLogLevel< logInfo::MyMessage >();
    }

    void MyGroup::outputMessage()
    {
        // effectively output the message taking into account the log level of the group instance
        GEOS_LOG_LEVEL_INFO( logInfo::MyMessage, "output some message" );
    }
