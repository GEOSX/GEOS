#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -A imcomp
#SBATCH -p pdebug

# Set fileName=xxx (no spaces), where the input file is pfw_input_xxx.py 

fileNames=(
	<INPUT_FILE_NAME_1>
	<INPUT_FILE_NAME_2>
	<...>
)
# ==========================================================================================================================================
# This should be the location of the input file and anything else you need to copy over:
fileLocation='<PATH_TO_INPUT_FILES>'

# This is where you want to run the simulation, which should be on a large parallel file system (lustre or workspace).  
# This directory should exist.  sub-directory with fileName will be created
runLocation='<PATH_TO_RUN_LOCATION>' # Quartz and Ruby

# This gets the username so the right localdefs file can be used.
userName="$(whoami)"

# You shouldn't need to modify any of this.  It checks if the fileName you provided isn't a null string (so you don't delete the parent directory!)
# then it copies files needed to run, changes directory, and runs the particle file writer.  That script will create a batch file and submit the job
# if those options are enabled:
#
# For certain input options (like using CT data) you'll also need to copy over any other needed files (CT data, input tables, etc.)
for fileName in "${fileNames[@]}"
do
	if [ -n "$fileName" ]
	then
		if [ $# -eq 0 ] || [ $1 -eq 1 ];
		then
			num_tasks="1"
			echo "Running job on 1 process: ${fileName}"
		else
			num_tasks="$1"
			echo "Running job on ${num_tasks} processes: ${fileName}"
		fi 

		aborted=false
		if [ -d $runLocation/$fileName/ ] && [ -z "$SLURM_JOBID" ];
		then
			echo "Directory ${runLocation}/${fileName} exists."
			while true; do
				read -p "Do you wish to overwite? " yn
				case $yn in
					[Yy]* ) echo "Overwriting..."; break;;
					[Nn]* ) echo "Aborted overwrite..."; aborted=true; break;;
					* ) echo "Please answer yes (Y/y) or no (N/n).";;
				esac
			done
		fi

		if [ $aborted = true ]
		then
			continue
		fi

		rm -rf $runLocation/$fileName/											# delete old results for the same fileName!!!
		mkdir -p $runLocation/$fileName/										# create the run/output directory

		cp $fileLocation/pfw_input_$fileName.py $runLocation/$fileName          # copy the input file
		cp particleFileWriter.py $runLocation/$fileName           				# copy the preprocessor
		cp pfw_check.py $runLocation/$fileName                    				# copy the autoRestart script
		cp pfw_geometryObjects.py $runLocation/$fileName          				# copy the geometry object functions
		cp userDefs_$userName.py $runLocation/$fileName           				# copy the local path information	

		cd $runLocation/$fileName                                 				# move to the run location
		if [ $# -eq 0 ] || [ $1 -eq 1 ];
		then
			python3 particleFileWriter.py pfw_input_$fileName
		else
			srun -n ${num_tasks} python3 particleFileWriter.py pfw_input_$fileName          # launch the VML
		fi
		echo # Print empty line for legibility
	fi
done
