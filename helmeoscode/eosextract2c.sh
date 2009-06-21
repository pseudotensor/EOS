# copy over the data to Matlab user
echo Begin transfering files
#DIR=/u1/matlab/eosdata/
DIR=/data/matlab/eosdata2/
LOCALDIR=./
ssh -x matlab@ki-rh42.slac.stanford.edu mkdir -p $DIR
scp $LOCALDIR/eos.dat $LOCALDIR/eos.head $LOCALDIR/eosparms.head matlab@ki-rh42.slac.stanford.edu:$DIR
echo Calling Matlab
# run matlab script to extract eos in HARM form
ssh -x matlab@ki-rh42.slac.stanford.edu sh /home/matlab/eosconvertrun.sh $DIR
#
echo End Matlab call
#
#
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosmonodegen.dat $LOCALDIR/eosmonodegen$1.dat
#
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosnew.dat $LOCALDIR/eosnew$1.dat
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosdegennew.dat $LOCALDIR/eosdegennew$1.dat
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosnew.head $LOCALDIR/eosnew$1.head
#
echo End transfering files
