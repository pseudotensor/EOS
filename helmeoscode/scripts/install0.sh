#!/bin/bash
# example script to run many EOS jobs on Orange cluster at KIPAC
#
# NOW follow (basic things done below after http line):
# http://harm.unfuddle.com/projects/5/notebooks/8/pages/37/latest

########## 1 ##########
#
# scp -rp jon@ki-rh42.slac.stanford.edu:"/data/jon/svneostest/helmeoscode/scripts/install0.sh" .
#

########### 2 ############
# Decide whether can use svn or have to use wget to get a potentially older copy of the source.  See below "svn checkout"

# whether to use svn or wget from old copy on Jon's website
usesvn=0

########### 3 ############
# Note that you can avoid installing GRB stellar model stuff if only want to generate HARM table.  Not bad idea if on (say) cluster where only generating table because that file is 1GB in size and takes 1 hour to download on Lonestar.
# This is compared to Shen/LS table stuff that only takes 6 minutes on Lonestar.

# whether to get and setup use of stellar model generation.  Not required for generating HARM table.
dostellarmodel=0

########### 4 ############
# run like: sh install0.sh  /lustre/ki/orange/jmckinne/
# can also just run like: sh install0.sh .
#


# 1) First, Change to directory with lots of space

# http://tldp.org/LDP/abs/html/testconstructs.html
if test -z "$1"
then
    echo "No directory given."
    exit 1
else
    echo "Directory is $1"
fi



cd $1
HEADDIR=`pwd`

mkdir eosfull
cd eosfull
BASEDIR=`pwd`
mkdir datadir
cd datadir
DATADIR=`pwd`

# 2) Get EOS source tree:
cd $BASEDIR
mkdir svneos ; cd svneos
SVNEOSDIR=`pwd`
if [ $usesvn -eq 1 ]
then
    svn checkout https://harm.unfuddle.com/svn/harm_joneos/ .
else
    # or use below if no SVN and believe the below source is updated enough
    wget http://www.slac.stanford.edu/~jmckinne/svneostest.tgz
    tar xvzf svneostest.tgz
    rm -rf svneostest.tgz
fi

cd joneoscode/jonmod
export JONEOSDIR=`pwd`
make clean ; make
#Done with test of compile jonmod

cd $BASEDIR
rm -rf sheneos.tables.orig_and_processed
mkdir sheneos.tables.orig_and_processed
cd sheneos.tables.orig_and_processed
SHENTABLEDIR=`pwd`
wget http://www.slac.stanford.edu/~jmckinne/sheneos.tables.orig_and_processed.tgz
tar xvzf sheneos.tables.orig_and_processed.tgz
rm -rf sheneos.tables.orig_and_processed.tgz

cd $SVNEOSDIR/helmeoscode
HELMDIR=`pwd`
rm -rf helm_table.dat sheneos.dat sheneosnonan.dat sheneos.tab sheneos.head sheneos.t00 sheneos.yp0 ls220_guesses.dat
ln -s $SHENTABLEDIR/helm_table.dat
ln -s $SHENTABLEDIR/sheneos.dat
ln -s $SHENTABLEDIR/sheneosnonan.dat
ln -s $SHENTABLEDIR/sheneos.tab
ln -s $SHENTABLEDIR/sheneos.head
ln -s $SHENTABLEDIR/sheneos.t00
ln -s $SHENTABLEDIR/sheneos.yp0
ln -s $SHENTABLEDIR/ls220_guesses.dat

rm -rf kazeos.parms.dek kazeos.loopvars.dek kazeos.loopvars2.dek kazeos.loopvars1.dek kazeos.loopparms.dek kazeos.dek kazeos3.dek kazeos2.dek kazeos1.dek const.dek
ln -s $JONEOSDIR/kazeos.parms.dek
ln -s $JONEOSDIR/kazeos.loopvars.dek
ln -s $JONEOSDIR/kazeos.loopvars2.dek
ln -s $JONEOSDIR/kazeos.loopvars1.dek
ln -s $JONEOSDIR/kazeos.loopparms.dek
ln -s $JONEOSDIR/kazeos.dek
ln -s $JONEOSDIR/kazeos3.dek
ln -s $JONEOSDIR/kazeos2.dek
ln -s $JONEOSDIR/kazeos1.dek
ln -s $JONEOSDIR/const.dek

make clean ; make

if [ $dostellarmodel -eq 1 ]
then
    # GRBMODEL STUFF
    mkdir $BASEDIR/svngrbmodel
    sh copy2grbmodel.sh $BASEDIR/svngrbmodel
    cd $BASEDIR/svngrbmodel
    SVNGRBMODELDIR=`pwd`
    wget http://www.slac.stanford.edu/~jmckinne/grb_stellarmodels.tgz
    tar xvzf grb_stellarmodels.tgz
    rm -rf grb_stellarmodels.tgz
    cp u/ki/jmckinne/research/ww95_stars/www.supersci.org/data/pre_sn_models/sol_metal/s251s7b\@14233.gz .
    gunzip s251s7b\@14233.gz
fi

# NOW follow (basic things done below after http line):
# http://harm.unfuddle.com/projects/5/notebooks/8/pages/40/latest

# 1)
cd $SVNEOSDIR/helmeoscode
# Change eosparms.f, kazeos.parms.dek, and kazeos.loopparms.dek as required
# Ensure kazeos.parms.dek has 1,1,0 for settings
# Ensure loop ranges are correct
make clean ; make
sh copyjonhelm.sh $DATADIR

# assume ready to create normal data in datadir:
cd $DATADIR


