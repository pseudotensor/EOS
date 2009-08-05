
c JCM:   Outline structure of code:
c
c
c  eos1point() calls subroutines and functions:
c     1) function xnuccalc(rhob,tk)
c     2) subroutine yetot_compute(rho10,T11,xnuc,yetot) (WILL BE)
c     3) subroutine muekT(rho10,T11,xnuc,etae)
c     4) function rnp(T11,etae)
c     5) function gammape(T11,etae)
c     6) subroutine computeyetot(rhob,tk,npratiofree,xnuc,yetot)
c     7) function pressure_e(rho10,T11,xnuc,etae,npratiofree)
c     8) function pN(rho10,T11,xnuc,etae,npratiofree)
c     9) function pressure_nu(rho10,T11,xnuc,etae,npratiofree)
c     10) function rhoe(T11,etae)
c     11) function ne(T11,etaev)
c
c     Before ending, eos1point() assumes at some point stored things in the following global quantities:
c
c     tautel, tauael, tautmu, tauamu, tauttau, tauatau
c
c     And eos1point() itself sets locally the following globals:
c
c     xnuc, yetot, etae, npratiofree
c     p_photon, p_eleposi, p_N, p_N, p_nu, p_tot
c     rho_photon, rho_eleposi, u_N, rho_N, rho_nu, rho_tot, u_tot
c     s_photon, s_eleposi, s_N, s_nu, s_tot
c     Qm
c     Qmel,Qmmu,Qmtau


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Other functions/subroutines that call other functions/subroutines:
c
c     muekT(): Generates etae: Calls tau_calc(etae), ne(etae), rnp(etae) in loop
c
c     rnp(): Uses results of tau_calc() and calls gammape(etae), gammane(etae), gamman(etae)
c
c     tau_calc(): Calls xnuccalc(),
c                 qdotNe(etae), qdotpair(etae), tauabs(Q,H), tauscatt(H),
c                 qdotNe_ele(etae), qdotNe_antiele(etae), 
c
c     In subs.f:
c
c     qdotNe(): Calls ne(etae), rnp(etae), qpe2(etae), qne2(etae)
c
c     qdotNe_ele(): Calls ne(etae), rnp(etae), qpe2(etae)
c
c     qdotNe_antiele(): Calls ne(etae), rnp(etae), qne2(etae)
c
c
c
c     Notice that muekT() relies on rnp() at many levels, and rnp() is complicated when trying to estimate thermalization or incomplete thermalization method.
c
c
c     Once muekT() is done, the following things are globally set as final quantities as determine by final \eta_e (etae)
c
c     tautel,tautmu,tauttau,tauael,tauamu,tauatau
c
c     New method with Y_e fixed makes rnp() much simpler
c
c
c
c



c     kaz_eos() assumes rhob,tk,hcm,tdynorye,computespecies are globally set
C======================================================================
      subroutine kaz_eos()
      implicit none
C======================================================================
c     Bring-in Kaz single global quantities
      include 'const.dek'
      include 'kazeos.parms.dek'
      include 'kazeos1.dek'
      include 'kazeos2.dek'
C======================================================================



cccccccccccccccccccccccc
c
c     Compute \eta_e (Once muekT is done, we have etae)
c     Assumes rhob,tk,hcm,tdynorye,tdynorynu globally set
c
ccccccccccccccccccccccc

c     kaz EOS assumes 0 degeneracy parameter for protons, neutrons, and neutrinos
      etan=0
      etap=0
      etanu=0



      call muekT(etae,etap,etan,etanu)


c      write(*,*) 'GOT HERE2',etae,etap,etan,etanu

cccccccccccccccccccccccc
c
c     Remaining things are functions of rho,T and etae, so any source for etae can be used
c     Assumes rhob,tk,hcm,tdynorye,tdynorynu globally set
c
ccccccccccccccccccccccc

      call kaz_physics_etae(etae,etap,etan,etanu)




      return
      end
















cccccccccccccccccccccccc
c
c
c     Computes xnuc, yetot, yefree, yebound, npratiofree, npratiobound, kazabar, kazzbar
c
c     For whichrnpmethod==0, assumes yefree and yetot already set
c     Once etae is known, this function is called to obtain solution state for given etae
c
c     For whichrnpmethod==1, then computes everything from yetot=tdynorye without etae
c
c
ccccccccccccccccccccccc

      subroutine kaz_species_etae(etae,etap,etan,etanu)
c      subroutine kaz_physics_etae(rho10,T11)
      implicit none

C======================================================================
c     Bring-in Kaz single global quantities
      include 'const.dek'
      include 'kazeos.parms.dek'
      include 'kazeos1.dek'
      include 'kazeos3.dek'
C======================================================================
      integer i,j,k,l           !dummy indeces
      real*8 tmp(100)
C======================================================================
      real*8 rho10              !matter density in 10^10 g/cc
      real*8 T11                !temperature in 10^11 K
      real*8 etae,etap,etan,etanu
C======================================================================
      real*8 T10                !temperature in 10^10 K
      real*8 tmev               !temperature in MeV
C======================================================================
      real*8 xnuc_small         !variable used in xnuc
      real*8 rnp                !external function
c======================================================================
      real*8 xnuccalc,Ypcalc        !function fraction of free nucleons
C======================================================================
      real*8 mutotxnuc0,mutot
C======================================================================
      real*8 nbtotal



c     Only should be here if computing species
      if(computespecies.eq.0) then
         return
      end if

cccccccccccccccccccccccc
c
c Set often-used variables
c
ccccccccccccccccccccccc
      T11 = tk/1.d11            !in 10^{11}K
      rho10=rhob/1.d10          !in 10^{10} g/cc
      T10 = T11*10.d0           !in 10^{10}K
      tmev=T10/(1.1605d0)       !in MeV


cccccccccccccccccccccccc
c
c Set xnuc
c
ccccccccccccccccccccccc

      xnuc=xnuccalc(rhob,tk)
c     xnuc_small=22.16d0*rho10**(-0.75d0)*T10**1.125d0*exp(-8.209d0/T10)
c     xnuc = min(1.d0, xnuc_small) !check1


cccccccccccccccccccccccc
c
c     Set npratiofree, yefree
c     Next set yetot, kazabar, kazzbar
c
ccccccccccccccccccccccc


      if(whichrnpmethod.eq.0) then

         call computemutotfit_xnuc0(rhob,tk,mutotxnuc0)
         abarbound = mutotxnuc0/amu*mb
         kazaheav=abarbound
         kazzheav=yeheav*kazaheav
         kazxheav=xnuc

         nbtotal = rhob/mb      ! temp var
         yeheav = (kazzheav/kazaheav)
         npratioheav = (1.0-yeheav)/yeheav
         npheav = nbtotal*kazxheav*yeheav
         nnheav = npheav*npratioheav


c     Rest presumes yetot set as guess or part of iteration or completely after iterations

c     yetot as input
         call computeyefreebound(rhob,tk,yetot,xnuc,yefree,yebound)
         call computemutotfit(rhob,xnuc,yefree,tk,mutot)
         kazabar = mutot/amu*mb
         kazzbar = yetot*kazabar

         npratiototal = (1.0-yetot)/yetot
         npratiobound = (1.0-yebound)/yebound
         yeheav=yebound

c     These can only be determined after rnp() iteratively or else assume yetot set
         nptotal = nbtotal*yetot
         npratiototal = (1.0-yetot)/yetot

         nntotal = nptotal*npratiototal
         npfree  = nbtotal*xnuc*yefree
         nnfree  = npfree*npratiofree
         npbound = nptotal-npfree
         nnbound = nntotal-nnfree
         if(npbound<0.0) then
            npbound=0.0
         end if
         if(nnbound<0.0) then
            nnbound=0.0
         end if

c     As partial iteration (if iterating) now recompute npratiofree using new parameters
c     If not iterating, then below changes nothing

c     Get optical depths used to determine rnp
         call tau_calc(whichynumethod,rho10,T11,etae,etap,etan,etanu)
c     Obtain n_n/n_p (free nucleons only)
         npratiofree=rnp(tdynorye,ntauael,ntaustotel,T11,etae,etap,etan,etanu)
         yefree = 1.0/(1.0+npratiofree)

c     finally close loop by computing yetot as new guess if iterating
         call computeyetot(rhob,tk,npratiofree,xnuc,yetot)




      else if(whichrnpmethod.eq.1) then



c     Compute yefree and yebound from tdynorye
c     GODMARK:
c     This yefree,yebound are Kaz-like definitions for his version
c     of nuclei EOS.  When using Kaz w/ nuclear EOS, this is NOT used

         yetot = tdynorye
         npratiototal = (1.0-yetot)/yetot

         call computeyefreebound(rhob,tk,yetot,xnuc,yefree,yebound)

         npratiofree=(1.0-yefree)/yefree
         npratiobound = (1.0-yebound)/yebound
         yeheav=yebound

         call computemutotfit_xnuc0(rhob,tk,mutotxnuc0)
         abarbound = mutotxnuc0/amu*mb
         kazaheav=abarbound
         kazzheav=yeheav*kazaheav
         kazxheav=xnuc

c     Get mutot
         call computemutotfit(rhob,xnuc,yefree,tk,mutot)


         kazabar = mutot/amu*mb
         kazzbar = yetot*kazabar


cccccccccccccccccc
c     Set species number density (should be consistent with similar computation in jon_lsbox.f)



         nbtotal = rhob/mb      ! temp var
         
         nptotal = nbtotal*yetot
         nntotal = nptotal*npratiototal
         
         npfree  = nbtotal*xnuc*yefree
         nnfree  = npfree*npratiofree
         
         npbound = nptotal-npfree
         nnbound = nntotal-nnfree
         if(npbound<0.0) then
            npbound=0.0
         end if
         if(nnbound<0.0) then
            nnbound=0.0
         end if
         
         yeheav = (kazzheav/kazaheav)
         npratioheav = (1.0-yeheav)/yeheav
         npheav = nbtotal*kazxheav*yeheav
         nnheav = npheav*npratioheav



      end if



      

      return
      end







cccccccccccccccccccccccc
c
c Once etae is known, this function is called to obtain solution state for given etae
c
c     Assumes xnuc set globally (e.g. either above for Kaz code or by nuclear EOS in HELM code)
c
c
c
ccccccccccccccccccccccc

      subroutine kaz_physics_neutrinos_etae(etae,etap,etan,etanu)
      implicit none

C======================================================================
c     Bring-in Kaz single global quantities
      include 'const.dek'
      include 'kazeos.parms.dek'
      include 'kazeos1.dek'
      include 'kazeos3.dek'
C======================================================================
      integer i,j,k,l           !dummy indeces
      real*8 tmp(100)
C======================================================================
      real*8 rho10              !matter density in 10^10 g/cc
      real*8 T11                !temperature in 10^11 K
      real*8 etae,etap,etan,etanu
C======================================================================
      real*8 T10                !temperature in 10^10 K
      real*8 tmev               !temperature in MeV
C======================================================================
      real*8 xnuc_small         !variable used in xnuc
C======================================================================
      real*8 xnuccalc,Ypcalc        !function fraction of free nucleons
      real*8 ne_cgs,nptest,ferror
C======================================================================





cccccccccccccccccccccccc
c
c Set often-used variables
c
ccccccccccccccccccccccc
      T11 = tk/1.d11            !in 10^{11}K
      rho10=rhob/1.d10          !in 10^{10} g/cc
      T10 = T11*10.d0           !in 10^{10}K
      tmev=T10/(1.1605d0)       !in MeV


c     DEBUG:
c      write(*,*) 'bad1',xnuc,rho10,kazabar,kazzbar,hcm





ccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute final opticals depths
c
ccccccccccccccccccccccccccccccccccccccccccccccc
      call tau_calc(whichynumethod,rho10,T11,xnuc,etae,etap,etan,etanu)

      

c     DEBUG:
c      write(*,*) 'etaposttaucalc=',etae,etap,etan,etanu





ccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute dY_e/dt (assumes tau_calc() above already called)
c     Also computes rates (gamma???global)
c
c     Also assumes yetot set as global quantity to call this function
c
c     For whichrnpmethod==0, need to have called rnp() first to get
c     thermalye
ccccccccccccccccccccccccccccccccccccccccccccccc




c      call dyedtcalc(tdynorye,yetot
c     1     ,gamman2pglobal,gammap2nglobal,dyedt,thermalye)






      return
      end










cccccccccccccccccccccccc
c
c Once etae is known, this function is called to obtain solution state for given etae
c
c     Requires computespecies be globally set
c
ccccccccccccccccccccccc

      subroutine kaz_physics_etae(etae,etap,etan,etanu)
c      subroutine kaz_physics_etae(rho10,T11)
      implicit none

C======================================================================
c     Bring-in Kaz single global quantities
      include 'const.dek'
      include 'kazeos.parms.dek'
      include 'kazeos1.dek'
      include 'kazeos3.dek'

c      integer computespecies ! whether to compute species

C======================================================================
      integer i,j,k,l           !dummy indeces
      real*8 tmp(100)
C======================================================================
      real*8 rho10              !matter density in 10^10 g/cc
      real*8 T11                !temperature in 10^11 K
      real*8 etae,etap,etan,etanu
C======================================================================
      real*8 T10                !temperature in 10^10 K
      real*8 tmev               !temperature in MeV
C======================================================================
      real*8 xnuc_small         !variable used in xnuc
      real*8 rnp                !external function
c======================================================================
      real*8 Gammaacc           !accretion timescale in sec
      real*8 ratioGG(3)         !Gammareac/Gammaacc
c      real*8 gammapeint         !rate p + e^- --> n + nue  in /sec
      real*8 gammape         !rate p + e^- --> n + nue  in /sec
C======================================================================
      real*8 Prad               !function photon pressure in erg/cc
      real*8 pressure_e         !function electron-positron pressure in erg/cc
      real*8 Pgas               !function nuclei pressure in erg/cc
      real*8 pN                 !function for general nucleon pressure in erg/cc
      real*8 pressure_nu        !function neutrino  pressure in erg/cc
C======================================================================
      real*8 ne                 !function 
      real*8 rhoe               !function MeV^4
      real*8 kine               !kinetic energy of e+d- in MeV^4
C======================================================================
      real*8 etaev              !virtual etae
C======================================================================
      real*8 xnuccalc,Ypcalc        !function fraction of free nucleons
      real*8 ne_cgs,nptest,ferror
C======================================================================
c      real*8 gammapeint         !rate p + e^- --> n + nue  in /sec
c      real*8 gammaneint         !rate n + e^+ --> p + nuebar in /sec
c      real*8 gammanint          !rate n --> p + e ^- + nuebar in /sec
C======================================================================
      real*8 dyedtcalc ! function to compute dY_e/dt [total Y_e]

      real*8 nq,mion,nion



cccccccccccccccccccccccc
c
c Set often-used variables
c
ccccccccccccccccccccccc
      T11 = tk/1.d11            !in 10^{11}K
      rho10=rhob/1.d10          !in 10^{10} g/cc
      T10 = T11*10.d0           !in 10^{10}K
      tmev=T10/(1.1605d0)       !in MeV




c     if computespecies.eq.0 then assume entered this function with
c     xnuc, yetot, yefree, yebound, npratiofree, kazabar, kazzbar
c     already computed into global quantities
      call kaz_species_etae(etae,etap,etan,etanu)


ccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute neutrino-related quantities
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccc
      call kaz_physics_neutrinos_etae(etae,etap,etan,etanu)



cccccccccccccccccccccccc
c
c Compute Pressures (non-neutrino)
c
ccccccccccccccccccccccc


c----------------------------------------------------------------------
c     
c     p_photon=Prad(rho10,T11)
      p_photon=PI**2/30.d0*2.d0*tmev**4 * mev4toecc /3.d0 !in erg/cc
      p_eleposi=pressure_e(rho10,T11,xnuc,etae,npratiofree) !erg/cc
c     
c     p_N=Pgas(rho10,T11)*(1.d0+3.d0*xnuc)/4.d0   !in erg/cc
c     Modified by JCM using Kaz's new formula
c     Below 2 in MeV
c     Ebin = 28.3
c     mN = 938.919
c     Ebinrat = Ebin/(4.d0*mN)
c     dissfactor=(1.d0+(3.d0-Ebinrat)*xnuc)/(1.d0-Ebinrat)
c     p_N=Pgas(rho10,T11)*dissfactor/4.d0   !in erg/cc
c     p_N=Pgas(rho10,T11)*(1.d0+3.d0*xnuc)/4.d0   !in erg/cc
      p_N = pN(rho10,T11,xnuc,etae,npratiofree)

c     TOTAL:
      p_tot=p_photon + p_eleposi + p_N + p_nu



cccccccccccccccccccccccc
c
c     Compute total internal energies + rest-mass (not including baryons)
c     (non-neutrino)
c
ccccccccccccccccccccccc
 
c----------------------------------------------------------------------
      rho_photon = PI**2/30.d0*2.d0*tmev**4 * mev4toecc !in erg/cc
      rho_eleposi=rhoe(T11,etae) * mev4toecc !in erg/cc
      u_N = 1.5d0*p_N           !in erg/cc
      rho_N = rhob*g2erg + u_N  !in erg/cc
      

c     TOTAL:
      rho_tot = rho_photon + rho_eleposi + rho_N + rho_nu
c     Neglect rest-mass of positrons (assumes they are purely thermally generated), neutrinos, (photons)
      u_tot = rho_photon + rho_eleposi + u_N + rho_nu



cccccccccccccccccccccccc
c
c     Compute total entropy densities
c     (non-neutrino)
c
ccccccccccccccccccccccc

c----------------------------------------------------------------------
      s_photon = 4.d0*PI**2/45.d0 * tmev**3*mev3ergKcc/kerg !in erg/K/kb/cc
c     JCM: was used to avoid issue with print out, but that's fixed now
c     if(s_photon<1E-30) s_photon=1E-30
      
      etaev=etae
      tmp(61)=rhoe(T11,etaev)/tmev
      tmp(62)=pressure_e(rho10,T11,xnuc,etaev,npratiofree)/mev4toecc/tmev
      tmp(63)=etae*ne(T11,etaev)
      s_eleposi = dim(tmp(61)+tmp(62)-tmp(63),0.d0)*mev3ergKcc/kerg !1/cc

c     Below does not account for other quantum term
      s_N= 2.50*p_N/mev4toecc/tmev*mev3ergKcc/kerg !in erg/K/kb/cc
c     JCM: Added below
c     See Shen EOS user's guide:
      nq = ((mb*clight*clight)*(kerg*tk)/(2.0d0*PI*(hbarcgs*clight)**2))**(3.0d0/2.0d0)
c     JCM: Below rho_N/mb should really be rho_N/m_N
c     Assuem kazabar already set
      mion = mb*kazabar
      nion = rho_N/mion
      s_N= S_N + rho_N*dlog(2.0d0*nq/nion)
c     JCM: Below scales correctly per baryon, but entropy density has different constant offset compared to Shen/LS, so using above
c      gammaindex=5.0d0/3.0d0
c      nindex=1.0d0/(gammaindex-1.0d0)
c      s_N= S_N + rho_N*dlog(p_N**nindex/rho_N**(nindex+1.0d0))

c     TOTAL:
      s_tot = s_photon + s_eleposi + s_nu + s_N !in erg/K/cc
c     s_ergKcc = s_tot*mev3ergKcc

c     DEBUG:
c      write(*,*) etae,tmev,T11,xnuc,npratiofree,rho10





      return
      end













