#
#DIRSOURCE=~/research/helm/
DIRSOURCE=`pwd`
DIRDEST=$1
mkdir -p $DIRDEST
cd $DIRDEST
rm -rf *.atb helm_table.dat ls220_guesses.dat sheneos.dat sheneos.head
ln -s $DIRSOURCE/*.atb $DIRSOURCE/helm_table.dat $DIRSOURCE/ls220_guesses.dat $DIRSOURCE/sheneos.dat $DIRSOURCE/sheneos.head .
cp -a $DIRSOURCE/*.exe .
cd $DIRSOURCE
#
