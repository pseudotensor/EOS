#!/bin/bash
# example script to run many EOS jobs on Orange cluster at KIPAC

########## 1 ##############
# Enter some HEAD directory
# First run install0.sh .
# sh ./install0.sh .
# this creates an eosfull directory.

########## 2 ##############
# Enter directory: eosfull/datadir/
#
# Then get these files either from repository or copy from a known location:
#
# scp -rp jon@ki-rh42.slac.stanford.edu:"/data/jon/svneostest/helmeoscode/scripts/chunkbunch.sh /data/jon/svneostest/helmeoscode/scripts/runchunkone.sh /data/jon/svneostest/helmeoscode/scripts/runchunkn.sh" .
# 

########## 3 ##############
# might want to remove old eoschunk directories:
# rm -rf eoschunkc*

# to kill lots of jobs:
# bkill -q normal 0

########## 4 ##############
# choose system/batch type
useorange=0
uselonestar=1

########## 5 ##############
# run this script by doing:
# sh chunkbunch.sh <DATADIR> <TOTALCHUNKS>
#e.g. 
# sh chunkbunch.sh /lustre/ki/orange/jmckinne/eosfull/datadir/ 100
# sh chunkbunch.sh . 160
#
# for 25,000 chunks, at 1 chunk per minute and 400 chunks will take 1 hour.  Doing 100 chunks will take about 4 hours if cluster free.



###########################
# Now begin to chunk



# Check that arguments given
# http://tldp.org/LDP/abs/html/testconstructs.html
if test -z "$1"
then
    echo "No directory given."
    echo "sh chunkbunch.sh <DATADIR> <TOTALCHUNKS>"
    exit 1
else
    echo "Directory is $1"
fi

if test -z "$2"
then
    echo "No TOTALCHUNKS given."
    echo "sh chunkbunch.sh <DATADIR> <TOTALCHUNKS>"
    exit 1
else
    echo "TOTALCHUNKS is $2"
fi


cd $1
export DATADIR=`pwd`
export HELMDIR=$DATADIR/../svneos/helmeoscode/
export SAFESTARTDIR=~/chunkdump/
export TOTALCHUNKS=$2


if [ $useorange -eq 1 ]
then
    truenumprocs=1
    MAXJOBS=200
fi
if [ $uselonestar -eq 1 ]
then
    truenumprocs=4
    # Note that lonestar has maximum of 40 queued jobs
    MAXJOBS=40
    # So TOTALCHUNKS/truenumprocs<40 has to be true
fi

numjobs=$(($TOTALCHUNKS/$truenumprocs))
if [ $numjobs -gt $MAXJOBS ]
then
    echo "Exceeded $MAXJOBS jobs with requesting of $numjobs jobs"
    exit 1
fi

# BSUB commands
# constant commands over chunks
jobprefix="eoschunk"
outputfile="stdout.out"
errorfile="stderr.err"


### LOOP OVER CHUNKS and start jobs
for CHUNK in `seq 1 $truenumprocs $TOTALCHUNKS`
do

    ##############
    # Setup safe directory to start from
    #
    # enter valid LSF directory
    mkdir -p $SAFESTARTDIR
    cd $SAFESTARTDIR
    chmod ug+x $DATADIR/runchunkone.sh
    chmod ug+x $DATADIR/runchunkn.sh

    # then jobnumber indicates starting chunk when doing multiple subchunks
    jobnumber=$CHUNK
    jobname=${jobprefix}c${jobnumber}tc${TOTALCHUNKS}

    #############
    # Start run
    #
    if [ $useorange -eq 1 ]
    then
        ###############
        # Setup job number/name/directory
	JOBDIR=${DATADIR}/${jobname}
	#
	bsub -n 1 -x -R span[ptile=1] -q kipac-ibq -J $jobname -o $outputfile -e $errorfile -a openmpi $DATADIR/runchunkone.sh $JOBDIR
    fi
    if [ $uselonestar -eq 1 ]
    then
	# http://www.tacc.utexas.edu/services/userguides/lonestar/
	# run with script using: bsub < bsub.batch
        # run on command line: (e.g.)  bsub -I -n 4 -W 0:05 -q development ibrun ./a.out 
        # to check jobs: showq -u or bjobs
	# 4 cores at a time
        # queues are: serial,normal,high,hero,development
	# Program to run is: "ibrun ./a.out" for MPI/parallel run
	# ptile=1 always so only 1 job started, but with exclusive (-x) access to node.
	bsub -B -N -u jmckinne@stanford.edu -P TG-AST080025N -x -W 24:00 -n $truenumprocs -x -o $outputfile -e $errorfile -R span[ptile=1] -q normal -J $jobname $DATADIR/runchunkn.sh  $CHUNK $truenumprocs $jobprefix $jobnumber $TOTALCHUNKS $DATADIR $jobname $HELMDIR
    fi

    ############
    # go back to data directory
    cd $DATADIR
    # report that done with submitting job
    echo "Done with CHUNK=$CHUNK out of $TOTALCHUNKS"
done

# report done with submitting all jobs
echo "Done with $TOTALCHUNKS chunks"


