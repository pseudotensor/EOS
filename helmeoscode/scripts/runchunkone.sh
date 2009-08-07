#!/bin/bash

# This is used to start job in one directory but run in another.  Otherwise, Orange complains like:
# ERROR: Local output directory /lustre/ki/orange/jmckinne/eosfull/datadir/eoschunkc74tc100/ invalid for batch jobs.
#
OLDDIR=`pwd`
LOCALDATADIR=$1

cd $LOCALDATADIR
./helmeos.exe
cd $OLDDIR
# DONE
