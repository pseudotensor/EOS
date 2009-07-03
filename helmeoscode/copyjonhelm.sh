#~/bin/bash
# runjonhelm.sh <DATADIR>
PROG=helmeos.exe
DATADIR=$1
#
DIRSOURCE=`pwd`
DIRDEST=$1
mkdir -p $DIRDEST
cd $DIRDEST
#
ln -s $DIRSOURCE/eosextract*.sh $DIRSOURCE/sheneosinterp.sh $DIRSOURCE/gethelmmatlabs.sh $DIRSOURCE/copy2grbmodel.sh .
ln -s $DIRSOURCE/*.atb $DIRSOURCE/helm_table.dat $DIRSOURCE/ls220_guesses.dat $DIRSOURCE/sheneos.dat $DIRSOURCE/sheneos.head .
cp -a $DIRSOURCE/$PROG .


