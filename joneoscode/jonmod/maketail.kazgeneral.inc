muekT.o: muekT.f const.dek
	ifort $(FFLAGS) -c  muekT.f

kazeos_compute.o: kazeos_compute.f const.dek
	ifort $(FFLAGS) -c   kazeos_compute.f

rhoe.o: rhoe.f
	ifort $(FFLAGS) -c   rhoe.f

kine.o: kine.f
	ifort $(FFLAGS) -c   kine.f

rnp.o: rnp.f const.dek
	ifort $(FFLAGS) -c   rnp.f

tau_calc.o: tau_calc.f
	ifort $(FFLAGS) -c   tau_calc.f

xnuc_calc_new.o: xnuc_calc_new.f
	ifort $(FFLAGS) -c   xnuc_calc_new.f


yefit.o: yefit.f const.dek
	ifort $(FFLAGS) -c   yefit.f

pressure.o: pressure.f const.dek
	ifort $(FFLAGS) -c   pressure.f

subs.o: subs.f
	ifort $(FFLAGS) -c   subs.f

ne.o: ne.f const.dek
	ifort $(FFLAGS) -c   ne.f

qpe2.o: qpe2.f
	ifort $(FFLAGS) -c  qpe2.f

qne2.o: qne2.f
	ifort $(FFLAGS) -c   qne2.f

pe.o: pe.f
	ifort $(FFLAGS) -c   pe.f

gammape.o: gammape.f
	ifort $(FFLAGS) -c  gammape.f

gammane.o: gammane.f
	ifort $(FFLAGS) -c  gammane.f

gamman.o: gamman.f
	ifort  $(FFLAGS) -c  gamman.f

gammannu.o: gammannu.f
	ifort  $(FFLAGS) -c  gammannu.f

gammapenu.o: gammapenu.f
	ifort  $(FFLAGS) -c  gammapenu.f

gammapnu.o: gammapnu.f
	ifort  $(FFLAGS) -c  gammapnu.f


kazloopfunctions.o: kazloopfunctions.f
	ifort  $(FFLAGS) -c  kazloopfunctions.f
