#!/bin/bash

#ARGS (1) file 1 (2) file 2 (3) optional: argument to silodiff
if [ $# -lt 2 ]; then
   echo "usage: ARGS (1) first file (2) second file (3) optional: argument to silodiff"
   exit 1
fi 

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$SILODIFF_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export SILODIFF_EXEC=/usr/gapps/GEOS/external_libs_icc_dbg/bin/silodiff
    if [ -e "$SILODIFF_EXEC" ]; then
	echo "$0 setting default SILODIFF path to $SILODIFF_EXEC"
    else
	SILODIFF_EXEC=/usr/gapps/silo/4.8/chaos_5_x86_64_ib/bin/silodiff
	if [ -e "$SILODIFF_EXEC" ]; then
	    echo "$0 setting default SILODIFF path to $SILODIFF_EXEC"
	else
	    echo "$0 cannot find default SILODIFF path: $SILODIFF_EXEC"
	    exit 1
	fi
    fi
fi

if [ ! -f $1 ]; then
    echo $PWD
    echo "$1 (1) does not exist - cannot silodiff!"
    exit 1
fi
if [ ! -f $2 ];then
    echo $PWD
    echo "$2 (2) does not exist - cannot silodiff!"
    exit 1
fi

EXIT_STATUS=`$SILODIFF_EXEC -E _silolibinfo $1 $2 2>&1 | tee atstmplog | grep ERROR | wc | awk '{if($1 != 0){print "1"}else{print "0"}}'`
if [ "$EXIT_STATUS" -ne "0" ];then
    echo "$PWD - Failed due to error during silodiff; probably a problem with one of the files"
    cat atstmplog
    rm atstmplog
    exit 1
fi

#EXIT_STATUS=`$SILODIFF_EXEC -E _silolibinfo $1 $2 -R 1e-15 | tee atstmplog | wc | awk '{print $1}'`
EXIT_STATUS=`$SILODIFF_EXEC $3 -E _silolibinfo $1 $2 | tee atstmplog | wc | awk '{print $1}'`
if [ "$EXIT_STATUS" -ne "0" ];then
    echo "$PWD - Failed due to differences found"
    cat atstmplog
    rm atstmplog
    exit 1
else
    rm atstmplog
    exit 0
fi
