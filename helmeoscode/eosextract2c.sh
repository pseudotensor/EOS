#!/bin/bash
# eosextract2c.sh <remotedir> <eos-extension for new files>
# e.g. eosextract2c.sh remotedirtest1 test1
# copy over the data to Matlab user
echo Begin transfering files
REMOTEDIR=$1
LOCALDIR=./
ssh -x matlab@ki-rh42.slac.stanford.edu mkdir -p $REMOTEDIR
scp $LOCALDIR/eos.dat $LOCALDIR/eos.head $LOCALDIR/eosparms.head matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR
echo Calling Matlab
# run matlab script to extract eos in HARM form
ssh -x matlab@ki-rh42.slac.stanford.edu sh /home/matlab/eosconvertrun.sh $REMOTEDIR
#
echo End Matlab call
#
#
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosmonodegen.dat $LOCALDIR/eosmonodegen$2.dat
#
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosnew.dat $LOCALDIR/eosnew$2.dat
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosdegennew.dat $LOCALDIR/eosdegennew$2.dat
scp matlab@ki-rh42.slac.stanford.edu:$REMOTEDIR/eosnew.head $LOCALDIR/eosnew$2.head
#
echo End transfering files
