DIR=50rhox50tkx50yex1ynux1h
#DIR=100rhox50tkx50yex20ynux1h.part6
#DIR=1test
#DIR=200x200x1x50x1.ynumethod2
#DIR=100x50x1x50x20
#DIR=200x200x1x50x40
#DIR=200x200x1x50x1
#
PROG=helmeos.$DIR.exe
#
make clean ; make
cp helmeos.exe $PROG
mkdir -p $DIR
#cp -a eosextract*.sh sheneosinterp.sh gethelmmatlabs.sh copy2grbmodel.sh $DIR/
#cp -a *.atb $PROG helm_table.dat ls220_guesses.dat sheneos.dat sheneos.head $DIR/
cd $DIR/
ln -s ../eosextract*.sh ../sheneosinterp.sh ../gethelmmatlabs.sh ../copy2grbmodel.sh .
ln -s ../*.atb ../helm_table.dat ../ls220_guesses.dat ../sheneos.dat ../sheneos.head .
cp -a ../$PROG .
./$PROG &> output.txt
cd ..


