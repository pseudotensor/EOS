      subroutine solverhoT (rho10,T11,HoverR,Qm,Qp,Qadv,dQ,emd)

c     Calculates the properties of an NDAF in a GRB or SNe

c     Given the mass em of the central mass in solar mass units, the
c     mass accretion rate emdot in solar masses per year, and the radius
c     r in units of the Schwarzschild radius, all sent through the
c     common block /params/, the subroutine solves for the density rho10
c     in units of 10^10 g/cm^3, the temperature T11 in units of 10^11 K,
c     and other quantities such as the disk height, heating and cooling
c     rates, etc.  The solution is by bracketing and bisection, with
c     starting values rho10guess and T11guess supplied by the user via
c     common bloc /guess/

      implicit double precision (a-h,o-z)
      common/params/alpha1,em,emdot,r
      common/guess/rho10guess,T11guess,ierrsolveT

c     Start with supplied rho10guess and solve for T11

      rho10=rho10guess
      call solveT (T11,rho10,Qm,Qp,Qadv,dQ,emd,HoverR)
c      write (*,*) ' T11, rho10, Qm, Qp, Qadv, dQ, emd, H/R = ',
c     &     T11,rho10,Qm,Qp,Qadv,dQ,emd,HoverR

c     Compare the mass accn rate emd with the required emdot and change
c     rho10 appropriately to bracket the solution for rho10

      do 20 i=1,50
      emdold=emd
      rho10old=rho10
      if ((emd-emdot).lt.0.d0) then  !emdot: needed \dot{M}(t) at r=r_0
         rho10=rho10old*2.d0
      else
         rho10=rho10old/2.d0
      endif
      call solveT (T11,rho10,Qm,Qp,Qadv,dQ,emd,HoverR)
c      write (*,*) ' T11, rho10, Qm, Qp, Qadv, dQ, emd, H/R = ',
c     &     T11,rho10,Qm,Qp,Qadv,dQ,emd,HoverR
      if ((emd-emdot)*(emdold-emdot).lt.0.d0) go to 30
 20   continue

      write (*,*) ' problem in solverhoT ! '

 30   continue

c     The solution has been bracketed.  Now solve for rho10 by bisection

      if (rho10.lt.rho10old) then
         rho10lo=rho10
         emdlo=emd
         rho10hi=rho10old
         emdhi=emdold
      else
         rho10lo=rho10old
         emdlo=emdold
         rho10hi=rho10
         emdhi=emd
      endif

      eps1=0.10
      do while (eps1.gt.1.d-11)
c      do 40 i=1,50
c      rho10=0.5d0*(rho10lo+rho10hi)
      rho10=sqrt(rho10lo*rho10hi)
      call solveT (T11,rho10,Qm,Qp,Qadv,dQ,emd,HoverR)
c      write (*,*) ' T11, rho10, Qm, Qp, Qadv, dQ, emd, H/R = ',
c     &     T11,rho10,Qm,Qp,Qadv,dQ,emd,HoverR
      if ((emd-emdot).lt.0.d0) then
         rho10lo=rho10
         emdlo=emd
      else
         rho10hi=rho10
         emdhi=emd
      endif
c      if ((rho10hi-rho10lo).lt.1.d-6*rho10hi) go to 50
c 40   continue

c      write (*,*) ' problem in solverhoT ! '
      
      eps1=abs(rho10hi-rho10lo)/rho10lo

      end do

 50   continue

c     The solution has been obtained.  Now set rho10guess to the
c     solution so that the next call to this subroutine will have a
c     better starting guess.

      rho10guess=rho10

      return
      end



      subroutine solveT (T11,rho10,Qm,Qp,Qadv,dQ,emd,HoverR)

c     Given mass em, radius r and density rho10, solves for the
c     temperature T11 and the mass accretion rate emd.  Also returns
c     other quantities of interest.

      implicit double precision (a-h,o-z)
      common/params/alpha1,em,emdot,r
      common/guess/rho10guess,T11guess,ierrsolveT

      ierrsolveT=0

c     Start by setting T11 equal to the user-supplied guess T11guess

      T11=T11guess
      call emdQcalc (rho10,T11,Qm,Qp,Qadv,dQ,emd,HoverR)
c      write (*,*) ' T11, rho10, emd = ',T11,rho10,emd

c     Compare heating and cooling rates and change T11 so as to bracket
c     the solution when heating=cooling

      do 20 i=1,50
      dQold=dQ
      T11old=T11
      if (dQold.lt.0.d0) then
         T11=T11/2.d0
      else
         T11=T11*2.d0
      endif
      call emdQcalc (rho10,T11,Qm,Qp,Qadv,dQ,emd,HoverR)
c      write (*,*) ' T11, rho10, emd = ',T11,rho10,emd
      if (dQold*dQ.lt.0.d0) go to 30
 20   continue

      write (*,*) ' problem in solveT ! rho10, T11 = ',rho10,T11
      ierrsolveT=1

 30   continue

c     Solution has been bracketed.  Now solve by bisection.

      if (T11.lt.T11old) then
         T11lo=T11
         dQlo=dQ
         T11hi=T11old
         dQhi=dQold
      else
         T11lo=T11old
         dQlo=dQold
         T11hi=T11
         dQhi=dQ
      endif

      eps2=0.1d0
      do while (eps2.gt.1.d-11)
c      do 40 i=1,50
      T11=0.5d0*(T11lo+T11hi)
      call emdQcalc (rho10,T11,Qm,Qp,Qadv,dQ,emd,HoverR)
c      write (*,*) ' rho10, Qm, Qp, emd, H/R = ',rho10,Qm,Qp,emd,HoverR
      if (dQ.gt.0.d0) then
         T11lo=T11
         dQlo=dQ
      else
         T11hi=T11
         dQhi=dQ
      endif
c      if ((T11hi-T11lo).lt.1.d-6*T11hi) go to 50
c 40   continue

c      write (*,*) ' problem in solveT ! rho10, T11 = ', rho10,T11
      
      eps2=(T11hi-T11lo)/min(abs(T11lo),abs(T11hi))
      end do
      
 50   continue

c     Solution has been obtained.  Set T11guess=T11 so that the
c     subroutine will have a better initial guess when it is called
c     next.

c      write (*,*) ' solveT: rho10, T11, dQ, emd, HoverR = ',
c     &     rho10,T11,dQ,emd,HoverR
      T11guess=T11

      return
      end



      subroutine emdQcalc (rho10,T11,Qm,Qp,Qadv,dQ,emd,HoverR)

c     Given mass em, density rho10 and temperature T11, calculates the
c     heating and cooling rates and the local mass accretion rate

      implicit double precision (a-h,o-z)
      common/params/alpha1,em,emdot,r
      common/guess/rho10guess,T11guess,ierrsolveT
      common/tau/tautel,tautmu,tauttau,tauael,tauamu,tauatau
      common/tau2/tauael_ele,tauael_antiele
      common /check1/etae
      common /check2/xnuc,xnucprev,dr
      common /check3/deg
      common /check4/press
      common /check5/Qphotodis
      common /check6/HoverR0
      common /check7/Qcool

c======================================================================
      real*8 Pelposi            !function, e+e- pressure  in erg/cc

c======================================================================
      real*8 press(4)
      real*8 pmax

c=====================================================================
      real*8 kappar
      real*8 tautot
      real*8 Qrad

C======================================================================
      real*8 xnuccalc               !function fraction of free nucleons
c======================================================================
      real*8 Qcool(4)

c======================================================================
      real*8 tmev               !temperature in MeV
      real*8 tmp(100)

c======================================================================
      real*8 tauael_ele
      real*8 tauael_antiele
      real*8 qminusel_ele
      real*8 qminusel_antiele

      real*8 qdotNe_ele
      real*8 qdotNe_antiele

c======================================================================
C     Xnuc
      
      T10 = T11*10.d0           !in 10^{10}K
      tmev=T10/(1.1605d0)       !in MeV

c      xnuc_small = 30.97d0*rho10**(-0.75d0)*T10**1.125d0*ex(-6.096d0/T10)
c      xnuc_small = 22.16d0*rho10**(-0.75d0)*T10**1.125d0*ex(-8.209d0/T10)
c      xnuc_small = 26.2d0*rho10**(-0.75d0)*tmev**1.125d0*ex(-7.074d0/tmev)

c      xnuc = min(1.d0, max(xnuc_small,9.d-2))     !check1

      xnuc=xnuccalc(rhob,tk)

c      tmp(1)=xnuc_small**2
c      xnuc = tmp(1)/(1.d0+tmp(1))+xnuc_small*1.d0/(1.d0+tmp(1))

c     Pressure due to radiation, gas and degeneracy


      press(1)=Prad(rho10,T11)
      press(2)=Pgas(rho10,T11)
      press(3)=Pelposi(etae,T11) !check2

c      press(3)=Pdeg(rho10,T11)


      pressure0=press(1)+press(2)+press(3) !check2
c      pressure0=Prad(rho10,T11)+Pgas(rho10,T11)+Pelposi(etae,T11) !check2
c      pressure0=Prad(rho10,T11)+Pgas(rho10,T11)+Pdeg(rho10,T11) !check2


c     Neutrino pressure is more complicated since it requires solving
c     the neutrino transfer problem.  We do this iteratively.

c     subroutine muekT calculates the quantity eta_e = mu_e/kT and the
c     neutrino degeneracy parameter, as defined in Kohri & Mineshige
c     (2002)

      Hold=0.d0
      HoverR0=Hold/(2.95d5*em*r)

      call muekT (rho10,T11,etae,deg)

c     Using the estimate of deg, calculate the cooling rate in each of
c     the three neutrino species

c      qminusmu=qdotpair(rho10,T11,deg)+qdotbrem(rho10,T11)
      qminusmu=qdotpair(rho10,T11,deg)*0.7d0/4.8d0+qdotbrem(rho10,T11) !check1

      qminustau=qminusmu
      qminusel=qdotNe(rho10,T11,etae,deg)+qdotplasmon(rho10,T11,etae)
c     &     +qminusmu
     &     +qdotpair(rho10,T11,deg)*3.4d0/4.8d0+qdotbrem(rho10,T11) !check1

c     Next we need the neutrino opacities.  We start by assuming zero
c     disk height and zero opacities, and solve for these quantities
c     iteratively.

      tauael=1.d-30
      tauaemu=1.d-30
      tauatau=1.d-30
      tautel=tauael
      tautmu=tauamu
      tauttau=tauatau
      Hold=0.d0

      eps3=0.10d0
      do while (eps3.gt.1.d-11)

c      do 10 i=1,100000000

c     Calculate current guess of neutrino pressures and then estimate
c     disk height H.  With the new H, calculate the neutrino opacities.
c     Keep iterating until H does not change any more.

c     The neutrino pressure is (1/3) neutrino energy density.  The
c     latter is obtained via the Popham & Narayan (1995) approximation,
c     given the absorptive and scattering opacities.

      pnuel=2.21d29*T11**4*(0.5d0*tautel+0.5774d0)/(0.5d0*tautel+0.5774d0
     &        +1.d0/(3.d0*tauael))
      pnumu=2.21d29*T11**4*(0.5d0*tautmu+0.5774d0)/(0.5d0*tautmu+0.5774d0
     &        +1.d0/(3.d0*tauamu))
      pnutau=2.21d29*T11**4*(0.5d0*tauttau+0.5774d0)/(0.5d0*tauttau+0.5774d0
     &        +1.d0/(3.d0*tauatau))
      pressure=pressure0+pnuel+pnumu+pnutau

      press(4)=pnuel+pnumu+pnutau


      call heightvisc(rho10,pressure,H,visc)
      HoverR=H/(2.95d5*em*r)

      tauael=tauabs(qminusel,T11,H) !check
      tauamu=tauabs(qminusmu,T11,H)
      tauatau=tauabs(qminustau,T11,H)

      taus=tauscatt(rho10,T11,H)

      tautel=tauael+taus
      tautmu=tauamu+taus
      tauttau=tauatau+taus
      
      qminusel_ele=qdotNe_ele(rho10,T11,etae,deg)
     &     + ( qdotplasmon(rho10,T11,etae)
     &     +qdotpair(rho10,T11,deg)*3.4d0/4.8d0
     &     +qdotbrem(rho10,T11) ) /2.d0 !check /2.d0

      qminusel_antiele=qdotNe_antiele(rho10,T11,etae,deg)
     &     + (+qdotplasmon(rho10,T11,etae)
     &     +qdotpair(rho10,T11,deg)*3.4d0/4.8d0
     &     + qdotbrem(rho10,T11) ) /2.d0   !check /2.d0

      tauael_ele=tauabs(qminusel_ele,T11,H)*2.d0 !because denominater/2.d0
      tauael_antiele=tauabs(qminusel_antiele,T11,H)*2.d0 ! because denominater/2.d0

      eps3=abs(H-Hold)/Hold
c      if (abs(H-Hold).lt.1.d-15*H) go to 20

      Hold=H

      end do


c 10   continue

c      write (*,*) ' H convergence failed ! '

c 20   continue

c     Neutrino solution has converged.  Sum the cooling rates ffrom each
c     neutrino species, calculated via the Popham & Narayan (1995)
c     formula, to get the total cooling rate.  Calculate also the
c     heating rate, 3 G M Mdot/8 pi R^3, and the advection rate, taken
c     to be equal to Qp*(H/R)**2.

      Qmel=6.62d39*T11**4/(0.5d0*tautel+0.5774+1.d0/(3.d0*tauael))
      Qmmu=6.62d39*T11**4/(0.5d0*tautmu+0.5774+1.d0/(3.d0*tauamu))
      Qmtau=6.62d39*T11**4/(0.5d0*tauttau+0.5774+1.d0/(3.d0*tauatau))

      Qm=Qmel+Qmmu+Qmtau
      
      Qcool(1)=Qm

      Qp=1.223d42*emdot/(em**2*r**3)
      Qadv=Qp*HoverR**2

      Qcool(2)=Qadv

      kappar=0.4d0+0.2023d-5*rho10/T11**3.5
      tautot=kappar*1.d10*rho10*H
      Qrad=5.67d39*T11**4 / tautot
      Qrad=0.d0

      Qcool(4)=Qrad

c      vr=3.d0/2.d0*visc/(2.95d5*em*r)       !in cm/sec check2 
      vr=visc/(2.95d5*em*r)       !in cm/sec check2 
c      dxdr=abs(xnuc-xnucprev)/(2.95d5*em*abs(dr))  !in /cm check1
      dxdr=abs(xnuc-xnucprev)/(2.95d5*em*r)  !in /cm check1

cc      if (T11.gt.0.7737d0) then
c      if (T11.gt.0.1d0) then
c       Qphotodis=0.68d29*rho10*vr*dxdr*H*dim(1.d0-xnuc,0.d0)
                                                  !in erg/cm^2/sec check1
c      else 
c         Qphotodis=0.d0         !in erg/cm^3/sec check1
c      end if

      Qphotodis=0.d0            !in erg/cm^3/sec check1
      Qcool(3)=Qphotodis

      dQ=Qp-Qadv-Qm-Qphotodis-Qrad
c      write (*,*) ' Qm, Qp, Qadv, dQ = ',Qm,Qp,Qadv,dQ

c     Finally, estimate the mass accretion rate using the relation: 
c           emd = 6 pi nu rho H

      emd=18.85d0*visc*1.d10*rho10*H/1.989d33  !\dot{M}=6 \pi \nu \Sigma
c      write (*,*) ' T11, dQ = ',T11,dQ

      return
      end



C=====================================================================
      function Prad (rho10,T11) !in erg/cc, only for photons

c     Radiation pressure: (11/12)aT^4

      implicit double precision (a-h,o-z)

c      Prad=6.93d29*T11**4       !in erg/cc for photons and electrons 

      Prad=6.93d29*T11**4 *2.d0/5.5d0  !check2 in erg/cc only for photons

      return
      end


C=====================================================================
      function Pgas (rho10,T11) !in erg/cc, only for nucleons

c     Gas pressure: rho k T/m_p (fully dissociated nuclei)

      implicit double precision (a-h,o-z)

      Pgas=8.26d28*rho10*T11    !in erg/cc

      return
      end



C=====================================================================
      function Pelposi(etae,T11) !in erg/cc for relativistic e+e-

c     Degeneracy pressure: assumes molecular wt per electron = 2 (equal
c     numbers of free protons and neutrons)

      implicit double precision (a-h,o-z)

C=====================================================================
c      real*8 mev4toecc           !1 MeV^4 = (mev4toecc) erg/cc
c      parameter(mev4toecc=2.085d26)

C=====================================================================
      real*8 peint              !funcion
      real*8 pe              !funcion

C=====================================================================
      real*8 tmp(20)

C=====================================================================

      include 'const.dek'

      Pelposi=peint(T11,etae)*mev4toecc  !check

c      Pelposi=pe(T11,etae)*mev4toecc  !check2

      

      return
      end

C=====================================================================
      function Pdeg (rho10,T11)

c     Degeneracy pressure: assumes molecular wt per electron = 2 (equal
c     numbers of free protons and neutrons)

      implicit double precision (a-h,o-z)

      Pdeg=1.06d28*rho10**1.333333d0

      return
      end


C=====================================================================




      subroutine heightvisc(rho10,pressure,height,visc)

c     Given pressure and density, solves for sound speed, disk height
c     and viscosity coefficient visc by means of the alpha prescription.

      implicit double precision (a-h,o-z)
      common/params/alpha1,em,emdot,r
      common/guess/rho10guess,T11guess,ierrsolveT

      OmegaK=7.19d4/(em*r**1.5d0)    !in sec^-1
      csq=pressure/(1.d10*rho10)   !in (cm/sec)^2

      height=sqrt(csq)/OmegaK   !in cm

      visc=0.1d0*alpha1*csq/OmegaK !in cm^2/sec

      return
      end




C=====================================================================
      subroutine muekT (rho10,T11,amuekT,deg)

c     Solves for the quantity eta_e = amuekT = mu_e/kT

C=====================================================================
      implicit double precision (a-h,o-z)

C=====================================================================
      real*8 neint              !electrn number density in MeV^3

C=====================================================================
      common /check2/xnuc,xnucprev,dr

C=====================================================================
c      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d32) !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d0) !1 MeV^3 = (mev3tocc) /cc

C=====================================================================
      real*8 rnp                !function

C=====================================================================

      include 'const.dek'

c     QnpkT = (m_n-m_p)c^2/kT = 1.29MeV/kT

      QnpkT=0.150d0/T11

c     Start with widely separated initial guess values so as to bracket
c     the solution

c      amuekTmin=-1.d100   !check2
      amuekTmin=1.d-5
c      fmin=dexp(amuekTmin-QnpkT)+1.d0
      fmin=rnp(T11,amuekTmin) + 1.d0
     |     -rho10*xnuc*6.02d33/(neint(T11,amuekTmin)*mev3tocc) !check2
c     |     -rho10*6.02d33/(neint(T11,amuekTmin)*mev3tocc) !check2

c      amuekTmax=1.d100
      amuekTmax=1.d12   !check2

c      fmax=dexp(amuekTmax-QnpkT)+1.d0
      fmax=rnp(T11,amuekTmax) + 1.d0
     |     -rho10*xnuc*6.02d33/(neint(T11,amuekTmax)*mev3tocc) !check2
c     |     -rho10*6.02d33/(neint(T11,amuekTmax)*mev3tocc) !check2

c      fmax=min(fmax,1.d30)
c      write (*,*) amuekTmin,fmin,amuekTmax,fmax
      if (fmin*fmax.gt.0.d0) then
         write (*,*) ' solution not bracketed '
         write(*,108) xnuc, fmin, fmax
 108     format(20(1pe11.3,' '))
         amuekT=0.d0
         return
      endif

c     Solve for emuekT iteratively.  The equations we solve are:
c             n_p = n_e = (mu_e/hbar c)^3/(3 pi^2)
c             exp(mu_e-Q) = n_n/n_p

      eps4=0.1d0
      do while (eps4.gt.1.d-11)

         amuekT=sqrt(amuekTmin*amuekTmax)
         
c         f=dexp(amuekT-QnpkT)+1.d0
         f=rnp(T11,amuekT) + 1.d0
     |     - rho10*xnuc*6.02d33/(neint(T11,amuekT)*mev3tocc) !check2
c     |     - rho10*6.02d33/(neint(T11,amuekT)*mev3tocc) !check2

         if (f.gt.0.d0) then
            amuekTmax=amuekT
            fmax=f
         else
            amuekTmin=amuekT
            fmin=f
         endif

         eps4=(amuekTmax-amuekTmin)/abs(amuekTmin)
      end do


 20   continue

c     Solution for amuekT has converged.  Estimate degeneracy parameter
c     d = 0.178rho10/T11^3.  When d>1, we have a degenerate electron gas.

      amuekT=sqrt(amuekTmin*amuekTmax)

c      deg=0.178d0*rho10/T11**3 
      deg=abs(amuekT)
      
      return
      end



      function qdotNe (rho10,T11,etae,deg)

c     Cooling by neutronization, or URCA, reactions.  Calculate the
c     asymptotic nondegenerate and degenerate rates and interpolate for
c     the current value of deg.  This cooling applies only to electron
c     neutrinos.

      implicit double precision (a-h,o-z)

      common /check2/xnuc,xnucprev,dr
c      common /rnp1/rnp

      real*8 qpe2int, qne2int
      real*8 neint
      real*8 np, nn             !in /cc
      real*8 rnp                !function
      real*8 gammapeint, gammaneint

C=====================================================================
c      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C=====================================================================

      include 'const.dek'

cc      qdotNenondeg=9.2d33*T11**6*rho10
c      qdotNenondeg=9.2d33*T11**6*rho10 * xnuc  !check1
c      qdotNedeg=1.1d31*etae**9*T11**9
c      qdotNe=qdotNenondeg/(1.+deg*deg)+deg*deg*qdotNedeg/(1.d0+deg*deg)

      np=neint(T11,etae) * mev3tocc !in /cc
      nn=np*rnp(T11,etae)

      qdotNe=(qpe2int(T11,etae)*np + qne2int(T11,etae)*nn) !check2


      return
      end



      real*8 function qdotNe_ele(rho10,T11,etae,deg)

c     Cooling by neutronization, or URCA, reactions.  Calculate the
c     asymptotic nondegenerate and degenerate rates and interpolate for
c     the current value of deg.  This cooling applies only to electron
c     neutrinos.

      implicit double precision (a-h,o-z)

      common /check2/xnuc,xnucprev,dr

      real*8 qpe2int, qne2int
      real*8 neint
      real*8 np, nn             !in /cc
      real*8 rnp                !function
      real*8 gammapeint, gammaneint

C=====================================================================
c      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C=====================================================================

      include 'const.dek'


      np=neint(T11,etae) * mev3tocc !in /cc
      nn=np*rnp(T11,etae)

      qdotNe_ele=qpe2int(T11,etae)*np  !check2


      return
      end


      real*8 function qdotNe_antiele(rho10,T11,etae,deg)

c     Cooling by neutronization, or URCA, reactions.  Calculate the
c     asymptotic nondegenerate and degenerate rates and interpolate for
c     the current value of deg.  This cooling applies only to electron
c     neutrinos.

      implicit double precision (a-h,o-z)

      common /check2/xnuc,xnucprev,dr

      real*8 qpe2int, qne2int
      real*8 neint
      real*8 np, nn             !in /cc
      real*8 rnp                !function
      real*8 gammapeint, gammaneint

C=====================================================================
c      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C=====================================================================

      include 'const.dek'

      np=neint(T11,etae) * mev3tocc !in /cc
      nn=np*rnp(T11,etae)

      qdotNe_antiele=qne2int(T11,etae)*nn   !check2


      return
      end



      function qdotpair (rho10,T11,deg)
      implicit double precision (a-h,o-z)

c     Cooling by pair annihilation processes.  This cooling applies to
c     all three neutrino species.  According to an email from Kaz Kohri,
c     the rate given in his paper is the sum of the rates from the three
c     neutrino species, so each has only a third of the rate.

      qdotpair=4.8d33*T11**9/(1.d0+deg*deg)    !total
c      qdotpair=1.6d33*T11**9/(1.+deg*deg)   !check1

      return
      end



      function qdotbrem (rho10,T11)
      implicit double precision (a-h,o-z)

c     Cooling by neucleon-nucleon bremsstrahlung.  Applies to each
c     neutrino species..

      qdotbrem=1.5d27*rho10*rho10*T11**5.5d0

      return
      end



      function qdotplasmon (rho10,T11,etae)
      implicit double precision (a-h,o-z)

c     Cooling by plasmon decay.  Applies only to electron neutrinos.

      gammap=3.213d-2*sqrt(9.8696d0 + 3.d0*etae*etae)
      qdotplasmon=1.5d32*T11**9*gammap**6*exp(-gammap)
     &     *(2.d0+2.d0*gammap+gammap*gammap)

      return
      end



      function tauabs (qminus,T11,H)

c     Given the cooling rate for a particular neutrino species,
c     calculated the absorptive opacity by: tau_a = q- H/(7/2) sigma T^4

      implicit double precision (a-h,o-z)

      tauabs=qminus*H/(1.98d40*T11**4)

      return
      end



      function tauscatt (rho10,T11,H)

c     Scattering opacity for each neutrino species.  Taken from di
c     Matteo et al. (2002).  According to Tiziana, the rate given in the
c     paper is for each neutrino species.

      implicit double precision (a-h,o-z)

      tauscatt=2.7d-7*rho10*T11*T11*H

      return
      end

