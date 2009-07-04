#all: eos.exe kazeos.o
all: eos.exe

KAZOBJS=muekT.o kazeos_compute.o rhoe.o kine.o rnp.o tau_calc.o xnuc_calc_new.o pressure.o subs.o ne.o qpe2.o qAe2.o qne2.o qn2.o pe.o gammape.o fermidirac.o gammaAe.o gammane.o gamman.o gammannu.o gammapenu.o gammapnu.o yefit.o kazloopfunctions.o

KAZSRCS=muekT.f kazeos_compute.f rhoe.f kine.f rnp.f tau_calc.f xnuc_calc_new.f pressure.f subs.f ne.f qpe2.f qAe2.f qne2.f qn2.f pe.f gammape.f fermidirac.f gammaAe.f gammane.f gamman.f gammannu.f gammapenu.f gammapnu.f yefit.f kazloopfunctions.f

DEPS=$(KAZOBJS) Makefile kazeos.loopparms.dek kazeos.loopvars.dek const.dek kazeos.dek kazeos1.dek kazeos2.dek kazeos3.dek kazeos.parms.dek

# DEBUG WITH BOUNDS CHECKS
#FFLAGS=-cm -132 -O0 -C -warn nounused -warn nouncalled -warn
# autoparallelization doesn't do well
#FFLAGS=-O2 -parallel -cm -132 -w90
FFLAGS=-O2 -cm -132 -w90
# DEBUG WITH BOUNDS CHECKS
#FFLAGS=-cm -132 -O0 -C
# NORMAL WITH OPTS
#FFLAGS=-cm -w90 -132

FC=ifort

eos.exe: eos.o $(DEPS)
	$(FC) $(FFLAGS) -o eos.exe eos.o $(KAZOBJS)

#kazeos.o: $(KAZOBJS)
#	$(FC) $(FFLAGS) -c $(KAZSRCS)
#	$(FC) $(FFLAGS) -c -o kazeos.o $(KAZOBJS)



clean:
	rm -f $(KAZOBJS) core* *~ eos.exe

eos.o: eos.f kazeos.loopparms.dek kazeos.loopvars.dek const.dek kazeos.dek kazeos.parms.dek
	$(FC) $(FFLAGS) -c   eos.f

include maketail.kazgeneral.inc


# end of makefile.



