# copy over the data to Matlab user
echo Begin Matlab call
scp eos.dat eos.head eosparms.head matlab@relativity:eosdata
# run matlab script to extract eos in HARM form
ssh -x matlab@relativity sh /home/matlab/eosconvertrun.sh /home/matlab/eosdata/
#
echo End Matlab call
#
#
scp matlab@relativity:eosdata/eosmonodegen.dat eosmonodegen$1.dat
#
scp matlab@relativity:eosdata/eosnew.dat eosnew$1.dat
scp matlab@relativity:eosdata/eosdegennew.dat eosdegennew$1.dat
scp matlab@relativity:eosdata/eosnew.head eosnew$1.head
#
