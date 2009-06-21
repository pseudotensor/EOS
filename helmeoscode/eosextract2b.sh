# copy over the data to Matlab user
echo Begin Matlab call
DIR=/u1/matlab/eosdata/
scp eos.dat eos.head eosparms.head matlab@ki-rh42.slac.stanford.edu:$DIR
# run matlab script to extract eos in HARM form
ssh -x matlab@ki-rh42.slac.stanford.edu sh /home/matlab/eosconvertrun.sh $DIR
#
echo End Matlab call
#
#
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosmonodegen.dat eosmonodegen$1.dat
#
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosnew.dat eosnew$1.dat
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosdegennew.dat eosdegennew$1.dat
scp matlab@ki-rh42.slac.stanford.edu:$DIR/eosnew.head eosnew$1.head
#
