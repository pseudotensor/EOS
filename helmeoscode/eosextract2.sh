# copy over the data to Matlab user
echo Begin Matlab call
scp eos.dat eos.head eosparms.head matlab@ki-rh42.slac.stanford.edu:eosdata
# run matlab script to extract eos in HARM form
ssh -x matlab@ki-rh42.slac.stanford.edu sh /home/matlab/eosconvertrun.sh /home/matlab/eosdata/
#
echo End Matlab call
#
#
scp matlab@ki-rh42.slac.stanford.edu:eosdata/eosmonodegen.dat eosmonodegen$1.dat
#
scp matlab@ki-rh42.slac.stanford.edu:eosdata/eosnew.dat eosnew$1.dat
scp matlab@ki-rh42.slac.stanford.edu:eosdata/eosdegennew.dat eosdegennew$1.dat
scp matlab@ki-rh42.slac.stanford.edu:eosdata/eosnew.head eosnew$1.head
#
