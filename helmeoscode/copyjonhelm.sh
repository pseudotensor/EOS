#~/bin/bash
# runjonhelm.sh <DATADIR>
DATADIR=$1
#
DIRSOURCE=`pwd`
DIRDEST=$1
mkdir -p $DIRDEST
cd $DIRDEST
#
ln -s $DIRSOURCE/eosextract*.sh $DIRSOURCE/sheneosinterp.sh $DIRSOURCE/gethelmmatlabs.sh $DIRSOURCE/copy2grbmodel.sh $DIRSOURCE/copyjonhelm.sh .
# don't just link chunk scripts since change for each job, so copy instead
alias cp='cp'
mkdir -p scripts/
cp -a $DIRSOURCE/scripts/*.sh scripts/
# make local link
ln -s scripts/*.sh .
ln -s $DIRSOURCE/*.atb $DIRSOURCE/helm_table.dat $DIRSOURCE/ls220_guesses.dat $DIRSOURCE/sheneos.dat $DIRSOURCE/sheneos.head .
cp -a $DIRSOURCE/helmeos.exe .
cp -a $DIRSOURCE/helmstareos.exe .
cp -a $DIRSOURCE/helmeosc .

alias cp='cp -i'

