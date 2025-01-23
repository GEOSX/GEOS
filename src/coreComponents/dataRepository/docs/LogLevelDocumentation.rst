.. _LogLevelDocumentation:

Log levels documentation
========================

Add a log level
---------------

To add a log level, you must respect the following structure and add it to the appropriate ``LogLevelsInfos.hpp`` :

.. code-block:: c++

    struct LogName
    {
    static constexpr int getMinLogLevel() { return logLevel; }
    static constexpr std::string_view getDescription() { return msg ; }
    };

If there is no ``LogLevelsInfos.hpp`` in the corresponding folder, you can create a ``LogLevelsInfos.hpp``

Example of usage
----------------

To use a log level, make sure it is declared in the constructor :

.. code-block:: c++

    addLogLevel< logInfo::StructName >();

To log a message with a log level, 4 macros are defined in LogLevelsInfo.hpp located in dataRepository :

* GEOS_LOG_LEVEL_INFO( logInfoStruct, msg )
* GEOS_LOG_LEVEL_INFO_RANK_0( logInfoStruct, msg ) 
* GEOS_LOG_LEVEL_INFO_BY_RANK( logInfoStruct, msg ) 
* GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfoStruct, msg ) 

.. code-block:: c++

    class Base : public dataRepository::Group
    {
        public:
        Base( string const & name );
        {
            addLogLevel< logInfo::Message >();
        }

        void outputMessage()
        {
            GEOS_LOG_LEVEL_INFO( logInfo::Message, "output some description" );
        }

    };