#!/bin/bash

# This is used to start job in one directory but run in another.  Otherwise, Orange complains like:
# ERROR: Local output directory /lustre/ki/orange/jmckinne/eosfull/datadir/eoschunkc74tc100/ invalid for batch jobs.
#
OLDDIR=`pwd`

# read-in arguments
CHUNKLIST=$1
TOTALCHUNKS=$2
DATADIR=$3
jobprefix=$4



# check that all old processes have died
NUMTRIES=10
SUCCESSKILL=0
for jj in `seq 1 $NUMTRIES`
do
    NUMPIDS=`pidof helmeos.exe | wc -w`
    if [ $NUMPIDS -eq 0 ]
	then
	SUCCESSKILL=1
	break
    else
	sleep 3
        # aggressively kill any instances of this program
	killall -s 9 helmeos.exe
	sleep 3
    fi
done

if [ $SUCCESSKILL -eq 0 ]
then
    echo "Could not destroy old processes."
    exit 1
fi


# Loop over subchunks
for SUBCHUNK in ${CHUNKLIST}
do
    subjobnumber=$SUBCHUNK
    subjobname=${jobprefix}c${subjobnumber}tc${TOTALCHUNKS}
    SUBJOBDIR=${DATADIR}/${subjobname}

    # change to subjob directory
    cd $SUBJOBDIR

    # Setup binary input parameters
    echo "$SUBCHUNK $TOTALCHUNKS" > eoschunk.dat

    # create unique file name for kill and ps lookup
    # NO, now use same filename so can easily identify when programs all stop
    #cp ./helmeos.exe ./helmeos.exe.$subjobname
    #killall helmeos.exe.$subjobname

    cd $SUBJOBDIR
    # in case didn't remove old directory, remove old files
    rm -rf eos.head eosdetails.dat eosother.dat eos.dat eoscoulomb.dat eosazbar.dat
    # run in background, so other subchunks can get started.
    nohup ./helmeos.exe &
    cd $OLDDIR
done


#############################################
# check that all old processes have completed
NUMTRIES=1000000
SUCCESSKILL=0
for jj in `seq 1 $NUMTRIES`
do
    NUMPIDS=`pidof helmeos.exe | wc -w`
    if [ $NUMPIDS -eq 0 ]
	then
	SUCCESSKILL=1
	break
    else
	# sleep for 1 minute before checking again
	sleep 60
    fi
done

if [ $SUCCESSKILL -eq 0 ]
then
    echo "Processes never stopped."
    exit 1
fi



# DONE

