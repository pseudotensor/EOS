#!bin/bash

# Script interpolates Shen EOS to intertesting form usable by jon_sheneos.f in helm code

MATLABUSER=matlab
#MATLABHOST=ki-rh42.slac.stanford.edu
MATLABHOST=relativity.cfa.harvard.edu
REMOTEDIR=/home/matlab/sheneosdata/

# copy over the data to Matlab user
echo Begin Matlab call for Shen EOS interp
scp sheneos.tab sheneos.yp0 sheneos.t00 ${MATLABUSER}@${MATLABHOST}:${REMOTEDIR}
# run matlab script to interpolate Shen eos
ssh -x ${MATLABUSER}@${MATLABHOST} sh /home/matlab/matlab/sheneosinterprunmatlab.sh ${REMOTEDIR}
#
echo End Matlab call to interpolation Shen EOS
#
#
scp ${MATLABUSER}@${MATLABHOST}:${REMOTEDIR}/sheneos.dat sheneos$1.dat
scp ${MATLABUSER}@${MATLABHOST}:${REMOTEDIR}/sheneos.head sheneos$1.head
#
