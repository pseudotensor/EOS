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
# scp -rp jon@ki-rh42.slac.stanford.edu:"/data/jon/svneostest/helmeoscode/scripts/chunkbunch.sh /data/jon/svneostest/helmeoscode/scripts/chunkn.sh /data/jon/svneostest/helmeoscode/scripts/runchunkone.sh /data/jon/svneostest/helmeoscode/scripts/runchunkn.sh" .
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
# sh chunkbunch.sh . 160 "1 13 14 15 16 17 33 49 53 54 55 56 65 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160"
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


if test -z "$3"
then
    #COUNTFAKECHUNKLIST=`seq 1 $truenumprocs $TOTALCHUNKS`
    echo "No MAINCHUNKLIST given, so generating chunks from 1 to TOTALCHUNKS=$TOTALCHUNKS in SUBCHUNKS of size truenumprocs=$trunumprocs comprising of $numjobs jobs"
    CHUNKPOOL=`seq 1 $TOTALCHUNKS`
else
    CHUNKPOOL=$3
    echo "Using TOTALCHUNKS=$TOTALCHUNKS but with arbitrary pool of CHUNKS rather than direct sequence"
fi






# BSUB commands
# constant commands over chunks
jobprefix="eoschunk"

# set execute permissions
chmod ug+x chunkn.sh

## get number of chunks and number of groups of chunks to do
echo "$CHUNKPOOL" > poollist.txt
numchunks=`wc -w poollist.txt | awk '{print $1}'`
numgroups=$(( ($numchunks / $truenumprocs)+1))


# check that number of jobs does not exceed maximum for this system
if [ $numgroups -gt $MAXJOBS ]
then
    echo "Exceeded $MAXJOBS jobs with requesting of $numgroups jobs"
    exit 1
fi


## setup array of pool of chunks to do
i=1
for fil in $CHUNKPOOL
do
    chunkarray[$i]=$fil
    i=$(($i+1))
done


### LOOP OVER CHUNKS and start jobs
ii=1
for CHUNK in `seq 1 $numgroups`
do

    # for name:
    jobnumber=$CHUNK

    # determine last chunk number
    #LASTCHUNK=$(($CHUNK+$truenumprocs-1))
    #CHUNKLIST=`seq $CHUNK $LASTCHUNK`

    # determine CHUNKLIST
    start=$ii
    finish=$(($start+$truenumprocs))
    CHUNKLIST=`awk '{for(j='$start';j<'$finish';j++)printf "%s ",$j;print ""}' poollist.txt`
    # increment $ii
    ii=$finish

    # run
    sh chunkn.sh "$CHUNKLIST" $TOTALCHUNKS $DATADIR $HELMDIR $SAFESTARTDIR $jobprefix $jobnumber $useorange $uselonestar $truenumprocs

    # report that done with submitting job
    echo "Done with CHUNK=$CHUNK out of $numgroups groups"
done

# report done with submitting all jobs
echo "Done with $numgroups groups in up to $TOTALCHUNKS chunks"



# "1 13 14 15 16 17 33 49 53 54 55 56 65 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160"
