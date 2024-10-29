#!/bin/bash


# Parse input arguments
INPUT_ARGS=""
INPUT_COUNTER=0
COMPILED_XML_NAME=""
COMPILED_XML_NAME_OVERRIDE=""
PARAMETER_ARGS=""
NEW_ARGS=""
USE_PYGEOSX=1
PYGEOS_WARNINGS=0
SCRIPT_DIR=$(dirname "$0")
PYGEOSX=$SCRIPT_DIR/../lib/PYGEOSX/bin/python

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -i|--input)
        INPUT_ARGS="$INPUT_ARGS -i $2"
        INPUT_COUNTER=$(( INPUT_COUNTER + 1 ))
        COMPILED_XML_NAME=$2.preprocessed
        shift
        ;;
        -c|--compiled_name)
        COMPILED_XML_NAME_OVERRIDE=$2
        shift
        ;;
        -p|--parameter)
        PARAMETER_ARGS="$PARAMETER_ARGS -p $2 $3"
        shift
        shift
        ;;
        -u|--use-pygeosx)
        USE_PYGEOSX=$2
        shift
        ;;
        -w|--pygeosx-warnings)
        PYGEOS_WARNINGS=$2
        shift
        ;;
        -h|--help)
        echo ""
        echo "Preprocessor options:"
        $SCRIPT_DIR/preprocess_xml --help
        echo ""
        echo "GEOSX options:"
        $SCRIPT_DIR/geosx --help
        exit
        ;;
        *)
        NEW_ARGS="$NEW_ARGS $key $2"
        shift
        ;;
    esac
    shift
done


# Choose the compiled xml name
if [ "$INPUT_COUNTER" -gt "1" ]
then
   COMPILED_XML_NAME="composite.xml.preprocessed"
fi

if [ ! -z "$COMPILED_XML_NAME_OVERRIDE" ]
then
   COMPILED_XML_NAME=$COMPILED_XML_NAME_OVERRIDE
fi


# Check for pygeosx
if [ "$USE_PYGEOSX" -eq "1" ]
then
   if [ -f $PYGEOSX ]
   then
      if [ "$PYGEOS_WARNINGS" -eq "1" ]
      then
         echo "Using pygeosx to preprocess the xml file"
      fi
   else
      if [ "$PYGEOS_WARNINGS" -eq "1" ]
      then
         echo "Pygeosx installation not found... reverting to non-pygeosx version"
      fi
      USE_PYGEOSX=0
   fi
fi


# Preprocess the xml files
if [ "$USE_PYGEOSX" -eq "1" ]
then
   $PYGEOSX $SCRIPT_DIR/pygeosx_preprocess.py $INPUT_ARGS $PARAMETER_ARGS $NEW_ARGS -c $COMPILED_XML_NAME -s $SCRIPT_DIR/../../src/coreComponents/schema/schema.xsd
else
   # As a backup, manually call the preprocessor and then continue with GEOSX
   $SCRIPT_DIR/preprocess_xml $INPUT_ARGS $PARAMETER_ARGS -c $COMPILED_XML_NAME -s $SCRIPT_DIR/../../src/coreComponents/schema/schema.xsd
   $SCRIPT_DIR/geosx -i $COMPILED_XML_NAME $NEW_ARGS
fi

