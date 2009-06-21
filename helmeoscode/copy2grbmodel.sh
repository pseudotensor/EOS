#
DIRSOURCE=~/research/helm/
DIRDEST=~/research/grbmodel/
cd $DIRDEST
ln -s $DIRSOURCE/*.atb $DIRSOURCE/helm_table.dat $DIRSOURCE/ls220_guesses.dat $DIRSOURCE/sheneos.dat $DIRSOURCE/sheneos.head .
cp -a $DIRSOURCE/*.exe .
cd $DIRSOURCE
#
