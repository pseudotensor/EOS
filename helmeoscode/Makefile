##############
#
#jon_helmstandard.f : Wrapper for full EOS tables
#jon_helm_stellarinput.f : Wrapper for stellar EOS tables
#
#eosparms.f : Setup which EOSs to use
#eosmemory.f : Memory allocated for processing EOS
#
#compute_nuclear_azbar.mex.f : Compute nuclear version of <A> and <Z> from component information (i.e. xprot, xneut, etc.)
#general_compute_functions.f : Various calculations used also in Matlab
#computeyerhob.f : Compute Y_e(rhob)
#jon_helm_azbarset.f : Compute <A>, <Z>, etc.
#enforce_x_consistency.mex.f : Matlab Mex function to check x=\sum x_i
#oldyeshecheck.f : Check Y_e
#enforce_ye_consistency.mex.f : Matlab Mex to check Y_e
#
#
#jon_eosfxt.f : Generates TIMMES EOS result
#gauss_fermi.f : TIMMES aux file
#jon_ls_2p7.f: Generates LS EOS result
#jon_sheneos.f : Generates Shen EOS result
#jon_helm.f : Generates HELM EOS result
#jon_helm_programstart.f : Top-level program
#
#jon_lsbox.f :  Combines Shen/LS with non-nuclear EOS
#
#jon_sneut5.f : Lattimer (Itoh) neutrino physics
#jon_kazQm.f : Simplified Kaz neutrino physics
#
#jon_helm_computekazlike.f : Reformat LS/Shen/HELM/TIMMES result into Kaz format
#jon_helm_outputstyle.f : Output formats for final finals and calls full neutrino physics function
#
#
#
#
#
#
#
#
#
#
##############
#
##############
#
# General instructions for KAZ/HELM/SHEN/MATLABINTERP code use
#
########################
# Debug mex function
########################
#
# 1) http://www.mathworks.com/support/solutions/data/1-17Z0R.html?solution=1-17Z0R
# 2) [compile mex's]
# e.g. mex -v -g  -output compute_nuclear_azbar /home/matlab/matlab/compute_nuclear_azbar.mex.f /home/matlab/matlab/general_compute_functions.f
# 3) matlab -Dgdb
# 4) run -nojvm
# 5) dbmex on
# 6) [Run mex command]
# e.g. [abarnum abar zbar yelocal] = compute_nuclear_azbar(0.9,0.1,0.0,0.1,1,1)
# 7) DO THINGS:
#    stop in compute_nuclear_azbar
#    continue
#    bt
# or:
#    break 95
#    continue
#    print PLHS
# 8) [see where in mex.f problem is and fix]
# 
# 9) check mex bit: mexext
#
########################
# Setup Shen EOS Table
########################
# 1) ki-rh39: cd /u/ki/jmckinne/research/helm/
# Using relativity because 32-bit and problem with 64-bit
# below uses research/eos/shen_interp.m that should be copied using gethelmmatlabs.sh on the system using Matlab
# 2) sh sheneosinterp.sh
# ?????????? below ?????????? -- seems there are no nan's or inf's anymore
# 3) sed 's/NaN/0.0/g' sheneos.dat > sheneosnonnan.dat ; cp sheneosnonnan.dat sheneos.dat
#
#
#
########################
# Setup HARM EOS Table
########################
#    Goto correct directory (helm directory is central site for use of KAZ/HELM/SHEN together):
# 1) cd /u/ki/jmckinne/research/helm/
# 1.5) [Modify eosparms.f, kazeos.parms.dek, kazeos.loopparms.dek]
#    Make combo code and create directory for processing (modify DIR=)
# 2) sh runjonhelm.sh
#    Create eosparms.head
# 3) eosparms.head filled with (e.g.):
# 1E-8 1E-8 1E-8 1E-14 1E-14 1E-14
# can copy from another directory (e.g. cp ../50rhox50tkx50yex1ynux1h.orig/eosparms.head .)
#    Ensure matlab server running
# 4) /data/opt/matlab/etc/lmstart
#    Copy data-directory to another computer (say from ki-rh39 to ki-rh42) before processing
# 5) scp -rp <dir> jon@ki-rh42:<otherdir>
#    Convert to HARM form
# 6) sh eosextract???.sh : Choose the version of script to run. Normal is (with mod of DIR=) (takes < 1 hour for 50x50x50x1x1):
#    sh eosextract2c.sh

########################
# Setup SM GRB Model
########################
# Separately for grbmodel directory:
# 0.5) [Modify eosparms.f, kazeos.parms.dek (1,0,1)]
# 1) make clean ; make
#    Setup grbmodel directory for looking at stellar data and generating inital stellar model data for HARM
# 2) sh copy2grbmodel.sh <SMGRBDESTDIR>    (if want to have helmstareos.exe and EOS tables to be in grbmodel directory)
#
#
# Purpose of other scripts:
# gethelmmatlabs.sh: If not on ki-rh39, then grabs matlab-related code so ki-rh39 is kept as central site for original code
# eostxt2bin.sh : old test
# sheneosinterp.sh: Interpolates Shen table to new form using Matlab (setup for relativity, but could be changed to ki-rh42)
# eosextract.sh: original Matlab processing of output from runjonhelm.sh script into form used by HARM (assumes Matlab code already there)
# eosextract2.sh : As eosextract.sh but for ki-rh42
# eosextract2b.sh : As eosextract2.sh but specifies DIR= where processing should occur on other side (DIR=/u1/matlab/eosdata/ or other names)
# eosextract2c.sh : As eosextract2b.sh but different specific DIR= used (DIR=/data/matlab/eosdata for ki-rh42)
##############


######################
# to compile and install and run:
#
# 0) Choose options like USEMPI, USEGFORTRAN, etc.
#
# 0.5) And of course setup options within code like described above (kazeos.parms.dek and kazeos.loopparms.dek)
#
# 1) Compile:
# cd ../joneoscode/jonmod/ ; make clean ; make ; cd ../../helmeoscode/ ; make clean ; make
#
# 2) Install:
# sh copyjonhelm.sh <installdir>
#
# 3a) Run in 1-CPU mode:
# cd <installdir>
#
# ./helmeos.exe
# OR ./helmeosc <totalchunks> <chunklist>
# OR ./helmstareos.exe
#
# Only helmeosc has MPI capability
#
# 3b) Run in MPI mode (only with helmeosc binary)
# cd <installdir>
# mpirun -np <numprocs> ./helmeosc <chunklist> <totalchunks> <DATADIR> <jobprefix>
# Where number of chunks in chunklist must be equal to numprocs.
# e.g.
# mpirun -np 4 ./helmeosc "1 2 3 4" 4 . "eoschunk"
#
#
# ensure no odd binaries are still running after mpirun failure: killall helmeosc
#
###########################


#############
# Choose whether to use MPI
USEMPI=1

# USEGFORTRAN==1 means any system with gfortran like ki-rh42 or any standard linux system
# USEGFORTRAN==0 means a system without gfortran like lonestar
USEGFORTRAN=1


all: helmeosc helmeos.exe helmstareos.exe

OBJS=jon_helm.o jon_helm_azbarset.o jon_sneut5.o jon_kazQm.o jon_helm_computekazlike.o jon_helm_outputstyle.o jon_helmstandard.o jon_eosfxt.o jon_lsbox.o general_compute_functions.o jon_ls_2p7.o jon_sheneos.o

#KAZDIR=~/research/kazeos/jonmod
KAZDIR=../joneoscode/jonmod/

KAZOBJS=$(KAZDIR)/muekT.o $(KAZDIR)/kazeos_compute.o $(KAZDIR)/rhoe.o $(KAZDIR)/kine.o $(KAZDIR)/rnp.o $(KAZDIR)/tau_calc.o $(KAZDIR)/fermidirac.o $(KAZDIR)/xnuc_calc_new.o $(KAZDIR)/pressure.o $(KAZDIR)/subs.o $(KAZDIR)/ne.o $(KAZDIR)/qpe2.o $(KAZDIR)/qAe2.o $(KAZDIR)/qne2.o $(KAZDIR)/qn2.o $(KAZDIR)/pe.o $(KAZDIR)/gammape.o $(KAZDIR)/gammaAe.o $(KAZDIR)/gammane.o $(KAZDIR)/gamman.o $(KAZDIR)/gammannu.o $(KAZDIR)/gammapenu.o $(KAZDIR)/gammapnu.o $(KAZDIR)/yefit.o $(KAZDIR)/kazloopfunctions.o

KAZDEPS=$(KAZDIR)/kazeos.loopparms.dek $(KAZDIR)/kazeos.loopvars.dek $(KAZDIR)/const.dek $(KAZDIR)/kazeos.dek

DEPS=Makefile $(OBJS) $(KAZOBJS) $(KAZDEPS) kaz_state.dek kaz_fakestate.dek



ifeq ($(USEGFORTRAN),1)
# DEBUG WITH BOUNDS CHECKS
#FFLAGS=-cm -132 -O0 -C -warn nounused -warn nouncalled -warn
FFLAGS=-cm -132 -w90 -g
# autoparallelization doesn't do well
#FFLAGS=-cm -parallel -132 -w90
#FFLAGS=-cm -132 -O0 -C
# NORMAL WITH OPTS
#FFLAGS=-cm -w90 -132
#FC=ifort


#FFLAGS=-ffixed-line-length-none -ff2c
FFLAGS=-ffixed-line-length-none -g
FC=gfortran 
#LINKER=gfortran -ff2c -lstdc++

ifeq ($(USEMPI),1)
# mpicc should work perfectly with linking as long as all paths are set as if source'ed relevant lib paths, include paths, etc.  Should be done by system!
LINKER=mpicc -cc=gfortran -lstdc++
endif
ifeq ($(USEMPI),0)
LINKER=gfortran -lstdc++
endif

# g++ works, but requires extern "C" stuff.  Also not standard on supercomputers.
ifeq ($(USEMPI),1)
CC=mpicc -cc=gcc
CFLAGS=-O2 -Wall -g -DUSINGMPI=$(USEMPI) -DUSINGGFORTRAN=$(USEGFORTRAN)
endif
ifeq ($(USEMPI),0)
#CC=g++
# gcc doesn't require extern "C" stuff, just normal extern prototype with _ at end of function name works.
# gcc also standard even on supercomputers and should have MPI support
CC=gcc
# icc doesn't work
#CC=icc
CFLAGS=-O2 -Wall -g -DUSINGMPI=$(USEMPI) -DUSINGGFORTRAN=$(USEGFORTRAN)
#LINKER=g++ $(CFLAGS) -lstdc++
endif


# end if USEGFORTRAN==1
endif


############ BEGIN USEGFORTRAN==0
ifeq ($(USEGFORTRAN),0)


# DEBUG WITH BOUNDS CHECKS
#FFLAGS=-cm -132 -O0 -C -warn nounused -warn nouncalled -warn
#FFLAGS=-cm -132 -w90 -g
# autoparallelization doesn't do well
#FFLAGS=-cm -parallel -132 -w90
#FFLAGS=-cm -132 -O0 -C
# NORMAL WITH OPTS
#FFLAGS=-cm -w90 -132
#FC=ifort


#FFLAGS=-ffixed-line-length-none -ff2c
FFLAGS=-ffixed-line-length-none -g
FC=f77
#LINKER=ifort -ff2c -lstdc++

ifeq ($(USEMPI),1)
# mpicc should work perfectly with linking as long as all paths are set as if source'ed relevant lib paths, include paths, etc.  Should be done by system!
LINKER=mpicc -cc=g77 $(FFLAGS) -lstdc++ -lg2c
endif
ifeq ($(USEMPI),0)
LINKER=g77 $(FFLAGS) -lstdc++ -lg2c
endif

# g++ works, but requires extern "C" stuff.  Also not standard on supercomputers.
ifeq ($(USEMPI),1)
CC=mpicc -cc=gcc
CFLAGS=-O2 -Wall -g -DUSINGMPI=$(USEMPI) -DUSINGGFORTRAN=$(USEGFORTRAN)
endif
ifeq ($(USEMPI),0)
#CC=g++
# gcc doesn't require extern "C" stuff, just normal extern prototype with _ at end of function name works.
# gcc also standard even on supercomputers and should have MPI support
CC=gcc
# icc doesn't work
#CC=icc
CFLAGS=-O2 -Wall -g -DUSINGMPI=$(USEMPI) -DUSINGGFORTRAN=$(USEGFORTRAN)
#LINKER=g++ $(CFLAGS) -lstdc++
endif


endif



#helmeosc: jon_helm_programstartc.o $(DEPS)
#	$(CC) $(CFLAGS) -o helmeosc jon_helm_programstartc.o $(OBJS) $(KAZOBJS)

helmeosc: jon_helm_programstartc.o $(DEPS)
	$(LINKER) -o helmeosc jon_helm_programstartc.o $(OBJS) $(KAZOBJS)

helmeos.exe: jon_helm_programstart.o $(DEPS)
	$(FC) $(FFLAGS) -o helmeos.exe jon_helm_programstart.o $(OBJS) $(KAZOBJS)

helmstareos.exe: jon_helm_stellarinput.o $(DEPS)
	$(FC) $(FFLAGS) -o helmstareos.exe jon_helm_stellarinput.o $(OBJS) $(KAZOBJS)



clean:
	rm -f $(OBJS) $(KAZOBJS) core* *~ helmeosc helmeos.exe helmstareos.exe

jon_helm_programstartc.o: jon_helm_programstartc.c Makefile
	$(CC)  -c $(CFLAGS)  jon_helm_programstartc.c

jon_helm_programstart.o: jon_helm_programstart.f jon_helm.f vector_eos.dek eosparms.f kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos4.dek   kazeos.loopparms.dek  kaz_state.dek  kaz_fakestate.dek const.dek jon_helmstandard.o jon_eosfxt.o Makefile
	$(FC)  -c $(FFLAGS)  jon_helm_programstart.f

jon_helm.o: jon_helm.f vector_eos.dek eosparms.f kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos4.dek   kazeos.loopparms.dek  kaz_state.dek  kaz_fakestate.dek const.dek jon_helmstandard.o jon_eosfxt.o Makefile
	$(FC)  -c $(FFLAGS)  jon_helm.f

jon_helm_azbarset.o: jon_helm_azbarset.f vector_eos.dek eosparms.f kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos4.dek   kazeos.loopparms.dek  kaz_state.dek const.dek Makefile
	$(FC) $(FFLAGS) -c   jon_helm_azbarset.f


jon_helm_stellarinput.o: jon_helm_stellarinput.f vector_eos.dek eosparms.f kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos4.dek   kazeos.loopparms.dek  kaz_state.dek const.dek jon_helmstandard.o jon_eosfxt.o Makefile
	$(FC) $(FFLAGS) -c   jon_helm_stellarinput.f

jon_helm_computekazlike.o: jon_helm_computekazlike.f vector_eos.dek eosparms.f kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos4.dek   kazeos.loopparms.dek kaz_state.dek const.dek Makefile
	$(FC) $(FFLAGS) -c   jon_helm_computekazlike.f

jon_helm_outputstyle.o: jon_helm_outputstyle.f vector_eos.dek eosparms.f kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos4.dek   kazeos.loopparms.dek  kaz_state.dek const.dek Makefile
	$(FC) $(FFLAGS) -c   jon_helm_outputstyle.f

jon_helmstandard.o: jon_helmstandard.f vector_eos.dek eosparms.f kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos4.dek   kazeos.loopparms.dek  kaz_state.dek const.dek Makefile
	$(FC) $(FFLAGS) -c   jon_helmstandard.f

jon_eosfxt.o: jon_eosfxt.f vector_eos.dek eosparms.f const.dek vector_sneos.dek eos_m4c.inc el_eos.inc force.inc maxwel.inc Makefile
	$(FC) $(FFLAGS) -c   jon_eosfxt.f

jon_ls_2p7.o: jon_ls_2p7.f vector_eos.dek eosparms.f const.dek vector_sneos.dek eos_m4c.inc el_eos.inc force.inc maxwel.inc Makefile
	$(FC) $(FFLAGS) -c   jon_ls_2p7.f

jon_lsbox.o: jon_lsbox.f vector_eos.dek eosparms.f const.dek vector_sneos.dek eos_m4c.inc el_eos.inc force.inc maxwel.inc Makefile
	$(FC) $(FFLAGS) -c   jon_lsbox.f


general_compute_functions.o: general_compute_functions.f vector_eos.dek eosparms.f const.dek vector_sneos.dek eos_m4c.inc el_eos.inc force.inc maxwel.inc Makefile
	$(FC) $(FFLAGS) -c    general_compute_functions.f

jon_sheneos.o: jon_sheneos.f  vector_eos.dek eosparms.f vector_sheneos.dek const.dek Makefile
	$(FC) $(FFLAGS) -c   jon_sheneos.f

#xnuc_calc_new.o: xnuc_calc_new.f eosparms.f const.dek Makefile
#	$(FC) $(FFLAGS) -c   xnuc_calc_new.f

jon_sneut5.o: jon_sneut5.f  eosparms.f const.dek Makefile
	$(FC) $(FFLAGS) -c   jon_sneut5.f

jon_kazQm.o: jon_kazQm.f  eosparms.f const.dek Makefile kaz_fakestate.dek  kaz_state.dek
	$(FC) $(FFLAGS) -c   jon_kazQm.f


#kazloopfunctions.o: kazloopfunctions.f
#	$(FC)  $(FFLAGS) -c  kazloopfunctions.f


#include maketail.kazgeneral.inc


# end of Makefile
