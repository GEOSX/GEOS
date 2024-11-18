#!/bin/bash

########################################
#ARGS (1) executable path (2) lscratch sub-directory location for working - ignored if LSCRATCH_WD set (3) min processes per job (4) max processes per job (5) minimum number of seconds to run (6) maximum number of seconds to run (7..N) required files to be copied
########################################
if [ $# -lt 7 ]; then
   echo "usage: ARGS (1) executable path (2) lscratch sub-directory location for working - ignored if LSCRATCH_WD set (3) min processes per job (4) max processes per job (5) minimum number of seconds to run (6) maximum number of seconds to run (7..N) required files to be copied" >&2
   exit 1
fi 

if [ -z "$SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export SCRIPT_PATH=~/scripts
    echo "setting default for SCRIPT_PATH: $SCRIPT_PATH" >&2
fi

########################################
#(1) SETUP YOUR WORKING DIRECTORY
########################################
#get the location at which to run (working directory)
#create if necesary
if [ -z "$LSCRATCH_WD" ]; then
    #if the working directory is not explicit, use the best throughput file system on the CZ
    export LSCRATCH_WD=`perl $SCRIPT_PATH/lc_filesystemsummary.pl | tail -n 1 | awk '{print $2}'`"/$USER/$2"
    echo "setting default for LSCRATCH_WD: $LSCRATCH_WD" >&2
fi
#copy files to the working directory
if [ -d $LCRATCH_WD ]; then
    rm -rf $LSCRATCH_WD
fi
mkdir -p $LSCRATCH_WD
for ((i=7; i<=$#; i++))
do
    cp ${!i} $LSCRATCH_WD 
done

########################################
#(2) DECIDE THE CLUSTER TO USE
########################################
#get the list of clusters for which I have banks ... cache this
if [ -e ~/_mshare ];then
    echo "using cached ~/_mshare" >&2
else
    mshare -p ALL -u $USER > ~/_mshare
fi
#get the configuration of the clusters for which I have banks ... cache this, too
if [ -e ~/_lcclusters ];then
    echo "using cached ~/_lcclusters" >&2
else
    for i in `cat ~/_mshare | awk '/^Partition/{print $2}'`
    do
	perl -e 'my $cluster = shift;if ($cluster eq "oslic") {exit 0;} open(NEWS,"news job.lim.$cluster|");while(<NEWS>){if(/pbatch/){my $line = $_;chomp($line);$line=~s/^\s+//;my @arr = split(/\s+/,$line); if(($arr[$#arr] eq "hours")||($arr[$#arr] eq "hrs")){my $tt = "weekend,weekday"; my $npj = 1;my $ppn = 12;if($cluster ne "aztec"){$ppn = 16; $npj = $arr[$#arr - 2];if(/weekday/){$tt = "weekday";}else{$tt = "weekend";}}print "$cluster\t$tt\t$npj\t$ppn\t$arr[$#arr - 1]\n";}}}close NEWS' -f $i >> ~/_lcclusters
    done
fi
#what day is it?
export WEEKDAY=`perl -e 'open(IN,"date|");my $line = <IN>; my @arr = split(/\s+/,$line); $line = $arr[0]; if($line eq "Sat"){print "0";} elsif($line eq "Sun"){print "0";}else{print "1"};close IN;'`
#get the best cluster to use
echo "running perl $SCRIPT_PATH/lc_availability.pl $3 $4 $5 $6 $WEEKDAY" >&2
export AVAILABILITY=`perl $SCRIPT_PATH/lc_availability.pl $3 $4 $5 $6 $WEEKDAY`
echo "result: $AVAILABILITY" >&2
export PARTITION=`echo $AVAILABILITY | awk '{print $1}'`
export NNODES=`echo $AVAILABILITY | awk '{print $2}'`
export NPROCS=`echo $AVAILABILITY | awk '{print $3}'`
export WTIME=`echo $AVAILABILITY | awk '{print $4}'`
#get the bank to use
export BANK=`perl -e 'my $bank = ""; my $nbank = 0; my $partition = shift; my $user = shift; open(IN,"<~/_mshare"); my $on = 0; while(<IN>){if(/^Partition/){if(/$partition/){$on==1}else{$on=0;}}elsif($on==1){if(/^$user/){my $line = $_; chomp($line); $line =~s/^\s+//; chomp($line); my @arr = split(/\s+/,$line); if($arr[3] > $nbank){$nbank = $arr[3]; $bank = $arr[1];}}}}close IN;print $bank;' -f $PARTITION $USER`
export INPUTFILENAME=`perl -e 'use File::Basename;my $path = shift;my $filename = fileparse($path); print $filename' -f $7`
echo "using input filename as $INPUTFILENAME" >&2

echo "#!/bin/csh" > $LSCRATCH_WD/msub
echo "#MSUB -N $2" >> $LSCRATCH_WD/msub
echo "#MSUB -l walltime=$WTIME" >> $LSCRATCH_WD/msub
echo "#MSUB -l nodes=$NNODES" >> $LSCRATCH_WD/msub
echo "#MSUB -l partition=$PARTITION" >> $LSCRATCH_WD/msub
echo "#MSUB -q pbatch" >> $LSCRATCH_WD/msub
echo "#MSUB -A $BANK" >> $LSCRATCH_WD/msub
echo "#MSUB -r n" >> $LSCRATCH_WD/msub
echo "#MSUB -m be" >> $LSCRATCH_WD/msub
echo "date > jstart" >> $LSCRATCH_WD/msub
echo "cd $LSCRATCH_WD" >> $LSCRATCH_WD/msub
echo "srun -n $NPROCS $1 -i $INPUTFILENAME > prunout" >> $LSCRATCH_WD/msub
echo "date > jend" >> $LSCRATCH_WD/msub

exit
