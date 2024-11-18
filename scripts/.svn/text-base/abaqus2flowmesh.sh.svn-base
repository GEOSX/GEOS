#!/bin/bash

#ARGS (1) Abaqus (2) fault definitions XML
if [ $# -ne 2 ]; then
   echo "usage: $0 <ABAQUS input filename> <fault definition file: space-delimited tuples of [x,y,z,length,depth,dip(deg),strike(deg)]>"
   exit 1
fi 
#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$FF_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export FF_EXEC=/usr/gapps/silo/4.8/chaos_5_x86_64_ib/bin/silodiff
    if [ -e "$FF_EXEC" ]; then
	echo "Setting default FlowMeshToNodeSet path to $FF_EXEC"
    else
	echo "Cannot find default FlowMeshToNodeSet path: $FF_EXEC"
	exit 1
    fi
fi
#DO INPUT FILES EXIST
if [ ! -f $1 ]; then
    echo $PWD
    echo "$1 (1) does not exist - cannot execute!"
    exit 1
fi
if [ ! -f $2 ];then
    echo $PWD
    echo "$2 (2) does not exist - cannot execute!"
    exit 1
fi
#START EXECUTING
$FF_EXEC 




EXIT_STATUS=`$SILODIFF_EXEC $1 $2 2>&1 | tee atstmplog | grep ERROR | wc | awk '{if($1 != 0){print "1"}else{print "0"}}'`
if [ "$EXIT_STATUS" -ne "0" ];then
    echo "$PWD - Failed due to error during silodiff; probably a problem with one of the files"
    cat atstmplog
    rm atstmplog
    exit 1
fi

EXIT_STATUS=`$SILODIFF_EXEC $1 $2 -A 1e-16 | tee atstmplog | wc | awk '{print $1}'`
if [ "$EXIT_STATUS" -ne "0" ];then
    echo "$PWD - Failed due to differences found"
    cat atstmplog
    rm atstmplog
    exit 1
else
    rm atstmplog
    exit 0
fi
