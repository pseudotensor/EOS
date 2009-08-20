#!bin/bash

# Script interpolates Shen EOS to intertesting form usable by jon_sheneos.f in helm code

# Run like: sh sheneosinterp.sh ""

echo "Using name prefix $1"

echo "Run this in helmeoscode dir where script is located.  Script changes to correct directory if standard installation."

echo "Ensure to have copied /sheneosinterprunmatlab.sh to Matlab user's ~/matlab/ directory"
echo "Also ensure matlabscripts/shen_interp.m is copied over to Matlab user's ~/matlab/ directory."
echo "Also ensure that latest functions that have a .mex.f cousin are in Matlab user's ~/matlab/ directory."

MATLABUSER=matlab
#MATLABHOST=ki-rh42.slac.stanford.edu
MATLABHOST=relativity.cfa.harvard.edu
REMOTEDIR=/home/matlab/sheneosdata/
# below is relative to helmeoscode/ directory
LOCALSHENDIR=../../sheneos.tables.orig_and_processed/

# change to Shen directory where original should be and processed tables will appear
cd $LOCALSHENDIR


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
