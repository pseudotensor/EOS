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
# NOTE: No longer have to do this if using copyjonhelm.sh
#
# scp -rp jon@ki-rh42.slac.stanford.edu:"/data/jon/svneostest/helmeoscode/scripts/chunkbunch.sh /data/jon/svneostest/helmeoscode/scripts/chunkn.sh /data/jon/svneostest/helmeoscode/scripts/runchunkone.sh /data/jon/svneostest/helmeoscode/scripts/runchunkn.sh /data/jon/svneostest/helmeoscode/scripts/collatechunks.sh" .
# 

########## 3 ##############
# might want to remove old eoschunk directories:
# rm -rf eoschunkc*

# to kill lots of jobs:
# bkill -q normal 0


########## 4 ##############
# run this script by doing:
# sh chunkbunch.sh <DATADIR> <TOTALCHUNKS> <SYSTEMTYPE>
#e.g. 
# sh chunkbunch.sh 2 /lustre/ki/orange/jmckinne/eosfull/datadir/ 100
# sh chunkbunch.sh 2 . 160
#
# sh chunkbunch.sh 2 . 160 "1 13 14 15 16 17 33 49 53 54 55 56 65 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160"
#
# for 25,000 chunks, at 1 chunk per minute and 400 chunks will take 1 hour.  Doing 100 chunks will take about 4 hours if cluster free.
#
#
# To operate in MPI mode, use system type 3 or 4 or 5
# In that case, number of chunks in CHUNKLIST (not totalchunks unless no CHUNKLIST given) is used as truenumprocs
#
# sh chunkbunch.sh 5 . 160
#


###########################
# Now begin to chunk



# Check that arguments given
# http://tldp.org/LDP/abs/html/testconstructs.html
if test -z "$1"
then
    echo "No system given."
    echo "sh chunkbunch.sh <SYSTEMTYPE> <DATADIR> <TOTALCHUNKS> <CHUNKLIST>"
    echo "<CHUNKLIST> is optional argument.  Only have to give 3 arguments."
    exit 1
else
    if [ $1 -eq 0 ]
    then
	useorange=0
	uselonestar=0
	USEMPI=0
	echo "System is local system"
    fi
    if [ $1 -eq 1 ]
    then
	useorange=1
	uselonestar=0
	USEMPI=0
	echo "System is Orange"
    fi
    if [ $1 -eq 2 ]
    then
	useorange=0
	uselonestar=1
	USEMPI=0
	echo "System is Lonestar"
    fi

    if [ $1 -eq 3 ]
    then
	useorange=0
	uselonestar=0
	USEMPI=1
	echo "System is local system using MPI"
    fi
    if [ $1 -eq 4 ]
    then
	useorange=1
	uselonestar=0
	USEMPI=1
	echo "System is Orange using MPI"
    fi
    if [ $1 -eq 5 ]
    then
	useorange=0
	uselonestar=1
	USEMPI=1
	echo "System is Lonestar using MPI"
    fi
fi


##############
if test -z "$2"
then
    echo "No directory given."
    echo "sh chunkbunch.sh <SYSTEMTYPE> <DATADIR> <TOTALCHUNKS> <CHUNKLIST>"
    echo "<CHUNKLIST> is optional argument.  Only have to give 3 arguments."
    exit 1
else
    INDIR=$2
    echo "Directory is $INDIR"
fi

##############
if test -z "$3"
then
    echo "No TOTALCHUNKS given."
    echo "sh chunkbunch.sh <SYSTEMTYPE> <DATADIR> <TOTALCHUNKS> <CHUNKLIST>"
    echo "<CHUNKLIST> is optional argument.  Only have to give 3 arguments."
    exit 1
else
    export TOTALCHUNKS=$3
    echo "TOTALCHUNKS is $TOTALCHUNKS"
fi



################
if test -z "$4"
then
    #COUNTFAKECHUNKLIST=`seq -s " " 1 $truenumprocs $TOTALCHUNKS`
    echo "No CHUNKLIST given, so generating chunks from 1 to TOTALCHUNKS=$TOTALCHUNKS"
    echo "No CHUNKLIST given, which is ok.  If want, can give CHUNKLIST:"
    echo "sh chunkbunch.sh <SYSTEMTYPE> <DATADIR> <TOTALCHUNKS> <CHUNKLIST>"
    echo "<CHUNKLIST> is optional argument.  Only have to give 3 arguments."
    CHUNKPOOL=`seq -s " " 1 $TOTALCHUNKS`
else
    CHUNKPOOL=$4
    echo "Using TOTALCHUNKS=$TOTALCHUNKS but with arbitrary pool of CHUNKS rather than direct sequence"
fi






#########
# Change to data directory

cd $INDIR
export DATADIR=`pwd`
#export HELMDIR=$DATADIR/../svneos/helmeoscode/
# assume ran copyjonhelm.sh into datadir already
export HELMDIR=`pwd`
export SAFESTARTDIR=~/chunkdump/



##########
# get number of chunks and number of groups of chunks to do
echo "$CHUNKPOOL" > poollist.txt
numchunks=`wc -w poollist.txt | awk '{print $1}'`



##############
# Setup truenumprocs and MAXJOBS
#
if [ $USEMPI -eq 0 ]
then
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
    if [ $useorange -eq 0 ] &&
	[ $uselonestar -eq 0 ]
    then
        # Assume on 4-core node system like ki-rh42.slac.stanford.edu
	truenumprocs=4
	MAXJOBS=10000000
    fi

else
    # MPI mode:
    truenumprocs=$numchunks
    MAXJOBS=1
fi




# BSUB command stuff:
# constant commands over chunks
jobprefix="eoschunk"


# set execute permissions
chmod ug+x chunkn.sh


## get number of groups of chunks to do
# below is my version of ceil() for bash
numgroups=$(( $numchunks / $truenumprocs ))
testvar=$(( $numgroups*$truenumprocs ))
if [ $testvar -lt $numchunks ]
then
    numgroups=$(( $numgroups+1))
fi


# check that number of jobs does not exceed maximum for this system
if [ $numgroups -gt $MAXJOBS ]
then
    echo "Exceeded $MAXJOBS jobs with requesting of $numgroups jobs"
    exit 1
fi


## setup array of pool of chunks to do (NOT USED)
#i=1
#for fil in $CHUNKPOOL
#do
#    chunkarray[$i]=$fil
#    i=$(($i+1))
#done



### LOOP OVER CHUNKS and start jobs
ii=1
for CHUNK in `seq -s " " 1 $numgroups`
do

    # for name:
    jobnumber=$CHUNK

    # determine last chunk number
    #LASTCHUNK=$(($CHUNK+$truenumprocs-1))
    #CHUNKLIST=`seq -s " " $CHUNK $LASTCHUNK`

    # determine CHUNKLIST
    start=$ii
    finish=$(($start+$truenumprocs))
    CHUNKLIST=`awk '{for(j='$start';j<'$finish';j++)printf "%s ",$j;print ""}' poollist.txt`
    # increment $ii
    ii=$finish

    # run
    sh chunkn.sh "$CHUNKLIST" $TOTALCHUNKS $DATADIR $HELMDIR $SAFESTARTDIR $jobprefix $jobnumber $USEMPI $useorange $uselonestar $truenumprocs

    # report that done with submitting job
    echo "Done with CHUNK=$CHUNK out of $numgroups groups"
done

# report done with submitting all jobs
echo "Done with $numgroups groups in up to $TOTALCHUNKS chunks"



# "1 13 14 15 16 17 33 49 53 54 55 56 65 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160"
