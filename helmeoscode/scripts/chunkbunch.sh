#!/bin/bash
# example script to run many EOS jobs on Orange cluster at KIPAC

########## 1 ##############
# First run install0.sh
# sh ./install0.sh

########## 2 ##############
# get this file either from repository or copy from a known location:
#
# scp jon@ki-rh42:/data/jon/svneostest/helmeoscode/scripts/chunkbunch.sh jon@ki-rh42:/data/jon/svneostest/helmeoscode/scripts/runchunkone.sh .
# 

########## 3 ##############
# might want to remove old eoschunk directories:
# rm -rf eoschunkc*

########## 4 ##############
# run this script by doing:
# sh chunkbunch.sh <DATADIR> <TOTALCHUNKS>
#e.g. 
# sh chunkbunch.sh /lustre/ki/orange/jmckinne/eosfull/datadir/ 100
#
# for 25,000 chunks, at 1 chunk per minute and 400 chunks will take 1 hour.  Doing 100 chunks will take about 4 hours if cluster free.



# Now begin to chunk

export DATADIR=$1
export HELMDIR=$DATADIR/../svneos/helmeoscode/
export TOTALCHUNKS=$2

# BSUB commands
# constant commands over chunks
truenumprocs=1
jobprefix="eoschunk"
outputfile="stdout.out"
errorfile="stderr.err"


### LOOP OVER CHUNKS and start jobs
for CHUNK in `seq 1 $TOTALCHUNKS`
do
    ###############
    # Setup job number/name/directory
    jobnumber=$CHUNK
    jobname=${jobprefix}c${jobnumber}tc${TOTALCHUNKS}
    JOBDIR=${DATADIR}/${jobname}

    # Create new directory
    rm -rf $JOBDIR
    mkdir -p $JOBDIR
    
    ################
    #   Copy/link relevant files (must be in "source" directory where previously used copyjonhelm.sh in install0.sh).  So this creates links to links except binary is copied
    sh $HELMDIR/copyjonhelm.sh $JOBDIR

    # change to job directory
    cd $JOBDIR

    # Setup binary input parameters
    echo "$CHUNK $TOTALCHUNKS" > eoschunk.dat

    ##############
    # run eos code
    #
    # enter valid LSF directory
    mkdir -p ~/chunkdump/
    cd ~/chunkdump/
    bsub -n 1 -x -R span[ptile=$truenumprocs] -q kipac-ibq -J $jobname -o $outputfile -e $errorfile -a openmpi $DATADIR/runchunkone.sh $JOBDIR
    # go back to data directory
    cd $DATADIR
    # report that done with submitting job
    echo "Done with CHUNK=$CHUNK out of $TOTALCHUNKS with jobname=$jobname and jobnumber=$jobnumber"
done

# report done with submitting all jobs
echo "Done with $TOTALCHUNKS chunks"


