#!/bin/bash

# Used to collate chunk results given by chunkbunch.sh

# 
########### 1 ##############
# assume start in datadir

########### 2 ##############
#
# NO LONGER required with new copyjonhelm.sh:
# scp -rp jon@ki-rh42.slac.stanford.edu:"/data/jon/svneostest/helmeoscode/scripts/collatechunks.sh" .

########### 3 ##############
# run like:
#
# sh collatechunks.sh . <TOTALCHUNKS>
# e.g. sh collatechunks.sh . 160



# Check that arguments given
# http://tldp.org/LDP/abs/html/testconstructs.html
if test -z "$1"
then
    echo "No directory given."
    echo "sh collatechunks.sh <DATADIR> <TOTALCHUNKS>"
    exit 1
else
    echo "Directory is $1"
fi

if test -z "$2"
then
    echo "No TOTALCHUNKS given."
    echo "sh collatechunks.sh <DATADIR> <TOTALCHUNKS>"
    exit 1
else
    echo "TOTALCHUNKS is $2"
fi


cd $1
export DATADIR=`pwd`
export TOTALCHUNKS=$2
jobprefix="eoschunk"

rm -rf eos.final.head ; touch eos.final.head
rm -rf eos.final.dat ; touch eos.final.dat
rm -rf eosother.final.dat ; touch eosother.final.dat
rm -rf eoscoulomb.final.dat ; touch eoscoulomb.final.dat
rm -rf eosazbar.final.dat ; touch eosazbar.final.dat

### LOOP OVER CHUNKS and collate
for CHUNK in `seq -s " " 1 $truenumprocs $TOTALCHUNKS`
do

    # then jobnumber indicates starting chunk when doing multiple subchunks
    jobnumber=$CHUNK
    jobname=${jobprefix}c${jobnumber}tc${TOTALCHUNKS}
    JOBDIR=${DATADIR}/${jobname}

    cat $JOBDIR/eos.head > eos.final.head
    cat $JOBDIR/eos.dat >> eos.final.dat
    cat $JOBDIR/eosother.dat >> eosother.final.dat
    cat $JOBDIR/eoscoulomb.dat >> eoscoulomb.final.dat
    cat $JOBDIR/eosazbar.dat >> eosazbar.final.dat

    echo "Doing $CHUNK out of $TOTALCHUNKS"
done

#DONE
