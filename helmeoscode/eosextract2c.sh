#!/bin/bash
# eosextract2c.sh <docopytoremote> <remotedir> <eos-extension for new files>
# e.g. eosextract2c.sh 1 remotedirtest1 test1
# e.g. eosextract2c.sh 0 /data/matlab/remotedirtest5 ""

DOCOPY=$1
REMOTEDIR=$2
POSTFIX=$3
LOCALDIR=./

####################################################
#
# copy over the data to Matlab user
#
####################################################
if [ $DOCOPY -eq 1 ]
then
    echo Begin transfering EOS files
    ssh -x matlab@ki-rh42.slac.stanford.edu mkdir -p $REMOTEDIR
    scp $LOCALDIR/eos.dat $LOCALDIR/eos.head matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR
fi

echo Begin transfering eosparms.head
# always transfer eosparms.head
scp $LOCALDIR/eosparms.head matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR

####################################################
#
# Run Matlab
#
####################################################
echo Calling Matlab
# run matlab script to extract eos in HARM form
ssh -x matlab@ki-rh42.slac.stanford.edu sh /home/matlab/eosconvertrun.sh $REMOTEDIR
#
echo End Matlab call
#
####################################################
#
# Copy back Matlab-created data
#
####################################################
#
# monotonized original table:
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosmonodegen.dat $LOCALDIR/eosmonodegen${POSTFIX}.dat
#
# primary table:
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosnew.dat $LOCALDIR/eosnew${POSTFIX}.dat
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosdegennew.dat $LOCALDIR/eosdegennew${POSTFIX}.dat
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosnew.head $LOCALDIR/eosnew${POSTFIX}.head
#
# extras table:
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosextranew.dat $LOCALDIR/eosextranew${POSTFIX}.dat
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosextradegennew.dat $LOCALDIR/eosextradegennew${POSTFIX}.dat
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosextranew.head $LOCALDIR/eosextranew${POSTFIX}.head

#
echo End transfering files
