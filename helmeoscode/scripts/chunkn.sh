#!/bin/bash

# chunkn.sh starts a single job
# chunkn.sh is called by chunkbunch.sh to start many jobs.
# one can call chunkn.sh separately by doing:
# sh chunkn.sh "<list>" $TOTALCHUNKS $DATADIR $HELMDIR $SAFESTARTDIR $jobprefix $jobnumber $useorange $uselonestar $truenumprocs
# E.g.:
# sh chunkn.sh "1 49 33 17" 160 /work/01014/tg802609/eosfull/datadirfull/ /work/01014/tg802609/eosfull/svneos/helmeoscode/ ~/chunkdump/ "eoschunk" 666 0 1 4
# The "666" number is just for the batch job ID and can be anything.
# This is useful becomes for some reason code crashes and one can restart those chunks that failed.  One can check what failed by ensuring non-zero and regular size for eos.dat:
# ls -alrt eoschunkc*/eos.dat

CHUNKLIST=$1
TOTALCHUNKS=$2
DATADIR=$3
HELMDIR=$4
SAFESTARTDIR=$5
jobprefix=$6
jobnumber=$7
useorange=$8
uselonestar=$9
truenumprocs=${10}



    ##############
    # Setup safe directory to start from
    #
    # enter valid LSF directory
mkdir -p $SAFESTARTDIR
cd $SAFESTARTDIR
chmod ug+x $DATADIR/runchunkone.sh
chmod ug+x $DATADIR/runchunkn.sh

    # then jobnumber indicates starting chunk when doing multiple subchunks
jobname=${jobprefix}c${jobnumber}tc${TOTALCHUNKS}

# set stdout and stderr files
outputfilebase="stdout.out"
outputfile=$DATADIR/${outputfilebase}.${jobname}
errorfilebase="stderr.err"
errorfile=$DATADIR/${errorfilebase}.${jobname}



# Loop over subchunks to create directory and copy files.
# Want this done *BEFORE* batch job submitted so batch has all files necessary setup.  Maybe would change binary or other files after batch submitted and don't want to wait for queued jobs to start just to do that.

for SUBCHUNK in ${CHUNKLIST}
do
    subjobnumber=$SUBCHUNK
    subjobname=${jobprefix}c${subjobnumber}tc${TOTALCHUNKS}
    SUBJOBDIR=${DATADIR}/${subjobname}

    # Create new directory
    cd $DATADIR
    rm -rf $SUBJOBDIR
    mkdir -p $SUBJOBDIR

    ################
    #   Copy/link relevant files (must be in "source" directory where previously used copyjonhelm.sh in install0.sh).  So this creates links to links except binary is copied
    cd $DATADIR
    sh $HELMDIR/copyjonhelm.sh $SUBJOBDIR

    # copy over runtime scripts (NOTE: bsub job submitted for many SUBJOBDIR's)
    #cp $DATADIR/runchunkone.sh $SUBJOBDIR
    #cp $DATADIR/runchunkn.sh $SUBJOBDIR

done


# back to datadir before start job with everything in $SUBJOBDIR necessary to start job without anymore copying of general files
# Can't have bsub called in SUBJOBDIR so stderr/stdout are there instead of in $DATADIR since bsub called on behalf of multiple directories
cd $DATADIR
#cd $SUBJOBDIR


    #############
    # Start run
    #
if [ $useorange -eq 1 ]
then
        ###############
        # Setup job number/name/directory
    JOBDIR=${DATADIR}/${jobname}
	#
    bsub -n 1 -x -R span[ptile=1] -q kipac-ibq -J $jobname -o $outputfile -e $errorfile -a openmpi $JOBDIR/runchunkone.sh $JOBDIR
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
        # "$CHUNKLIST" has to be in quotes to preserve fact that single argument, otherwise gets expanded
    bsub -B -N -u jmckinne@stanford.edu -P TG-AST080025N -x -W 47:59 -n $truenumprocs -x -o $outputfile -e $errorfile -R span[ptile=1] -q normal -J $jobname $DATADIR/runchunkn.sh "$CHUNKLIST" $TOTALCHUNKS $DATADIR $jobprefix
fi

    ############
    # go back to data directory
cd $DATADIR
