C======================================================================
      real*8 function pressure_e(rho10,T11,xnuc,etae,npratiofree) !in erg/cc

C======================================================================
      implicit none

C======================================================================
c      common/tau/tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus
c      common/tau2/tauael_ele, tauael_antiele
c      common /variable1/ tdyn,hcm,rhob,tk

      include 'kazeos1.dek'
C======================================================================
c      real*8 tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus

C======================================================================
c      real*8 tauael_ele
c      real*8 tauael_antiele

C======================================================================
c      real*8 tdyn               !dynamical timescale in sec
c      real*8 hcm                !disk half thickness in cm
c      real*8 rhob               !baryon denstiy in g/cc
c      real*8 tk                 !temperature in K

C======================================================================
      real*8 rho10,T11,xnuc,etae,npratiofree
      real*8 yetotlocal
C======================================================================
      real*8 Pelposi            !function
      real*8 Pgas               !function

C======================================================================


      if ((etae.lt.1.d-6).and.(tk.lt.1.d6)) then
c     JCM : Kaz has a 4.0d0 instead of the correct 2.0d0
c rarely (if ever) used
c         pressure_e=Pgas(rho10,T11)*(xnuc/(1.d0+npratiofree)+(1.d0-xnuc)/2.0d0)
         call computeyetot(rhob,tk,npratiofree,xnuc,yetotlocal)
         pressure_e=Pgas(rho10,T11)*yetotlocal
      else
         pressure_e=Pelposi(etae,T11)
      end if

      return
      end





c     Below requires tau_calc() to be computed
c     Pressure
c      p_nu=pressure_nu(rho10,T11) !in erg/cc
c     Internal energy density
c      rho_nu = 3.d0 * p_nu      !in erg/cc
c     entropy density
c      s_nu=4.d0/3.d0*rho_nu/mev4toecc/tmev*mev3ergKcc/kerg !in erg/K/kb/cc
c
C======================================================================
c      real*8 function pressure_nu(rho10,T11) !in erg/cc
c
C======================================================================
c      implicit none
c
C======================================================================
c      include 'kazeos1.dek'
c      common/tau/tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus
c      common/tau2/tauael_ele, tauael_antiele
c      common /variable1/ tdyn,hcm,rhob,tk
c
C======================================================================
c      real*8 tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus

C======================================================================
c      real*8 tauael_ele
c      real*8 tauael_antiele
c
C======================================================================
c      real*8 tdyn               !dynamical timescale in sec
c      real*8 hcm                !disk half thickness in cm
c      real*8 rhob               !baryon denstiy in g/cc
c      real*8 tk                 !temperature in K
c
C======================================================================
c      real*8 rho10,T11
c
C======================================================================
C======================================================================
c      real*8 pnuel
c      real*8 pnumu
c      real*8 pnutau
c
C======================================================================
c      pnuel=2.21d29*T11**4*(0.5d0*tautel+0.5774d0)/(0.5d0*tautel+0.5774d0
c     &        +1.d0/(3.d0*tauael))
c      pnumu=2.21d29*T11**4*(0.5d0*tautmu+0.5774d0)/(0.5d0*tautmu+0.5774d0
c     &        +1.d0/(3.d0*tauamu))
c      pnutau=2.21d29*T11**4*(0.5d0*tauttau+0.5774d0)/(0.5d0*tauttau+0.5774d0
c     &        +1.d0/(3.d0*tauatau))
c
c      pressure_nu = pnuel + pnumu + pnutau
c
c      return
c      end




C=====================================================================
      real*8 function Prad (rho10,T11) !in erg/cc, only for photons

C=====================================================================
c     Radiation pressure: (11/12)aT^4

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 rho10,T11

C=====================================================================
c      Prad=6.93d29*T11**4       !in erg/cc for photons and electrons 

c      Prad=6.93d29*T11**4 *2.d0/5.5d0  !check2 in erg/cc only for photons
      Prad=2.52197d29*T11**4           !check2 in erg/cc only for photons

      return
      end


C=====================================================================
      real*8 function Pgas(rho10,T11) !in erg/cc, only for nucleons

C=====================================================================
c     Gas pressure: rho k T/m_p (fully dissociated nuclei)

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 rho10,T11

c      Pgas=8.26d28*rho10*T11    !in erg/cc
c   JCM (Pgas = kb/mb rhob T, where mb = (mn+mp)/2 assuming  npratiofree=1?)
c
c     
c
      Pgas=8.24902d28*rho10*T11    !in erg/cc



      return
      end


      real*8 function PN(rho10,T11,xnuc,etae,npratiofree) !in erg/cc, only for nucleons

C=====================================================================
c     Gas pressure with binding energy and species contribution included for alpha particles

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

c     Apparently in Fortran rhob,tk here are not globals if inside function rather than subroutine
c
      include 'kazeos1.dek'
c      common /variable1/ tdyn,hcm,rhob,tk
c      real*8 tdyn,hcm,rhob,tk  ! globals now
C=====================================================================
      real*8 Pgas ! function
      real*8 rho10,T11,xnuc,etae,npratiofree
c      real*8 Ebin,mN,Ebinrat,dissfactor,p_N
c      real*8 mn,mp,me,amu,AH,AHe,AO16,ASi28,AFe56
c      real*8 malpha
      real*8 temp,mutot,yefit
c      real*8 mutot1,mutot2,mutot3,mutot4,mutot5
c      real*8 tk12,tk23,tk34,tk45
c
c
      include 'const.dek'

c          p_N=Pgas(rho10,T11)*(1.d0+3.d0*xnuc)/4.d0   !in erg/cc
c         Modified by JCM using Kaz's new formula
c         Below 2 in MeV
c           Ebin = 28.3
c           mN = 938.919
c           Ebinrat = Ebin/(4.d0*mN)
c           dissfactor=(1.d0+(3.d0-Ebinrat)*xnuc)/(1.d0-Ebinrat)
c          PN=Pgas(rho10,T11)*dissfactor/4.d0   !in erg/cc
c          p_N=Pgas(rho10,T11)*(1.d0+3.d0*xnuc)/4.d0   !in erg/cc


c      call computemutotfit(rhob,tk,mutot)

c      JCM way:
c     Just sum up the special contribution to the equation:
c     utot = 3/2*kb*T*\Sum_i(rho_i/m_i)
c     pN = (gamma-1) utot
c     GODMARK:
c     Assumes A=1 if xnuc=1, so won't work for nucleon EOS in general
c     Need to make function of A and Z or use HELM abar/zbar calculation
c      pN=Pgas(rho10,T11)*(xnuc+(1.0-xnuc)*(mb/(kazabar*amu)))

c     Assume kazabar in terms of mb, not amu, and Pgas() is setup in terms of mb
      pN=Pgas(rho10,T11)*(1.0/kazabar)




      return
      end


C=====================================================================
      real*8 function Pelposi(etae,T11) !in erg/cc for relativistic e+e-

C=====================================================================
c     Degeneracy pressure: assumes molecular wt per electron = 2 (equal
c     numbers of free protons and neutrons)

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,etae

C=====================================================================
c      real*8 mev4toecc           !1 MeV^4 = (mev4toecc) erg/cc
c      parameter(mev4toecc=2.085d26)

C=====================================================================
      real*8 pe              !funcion
C=====================================================================
      real*8 tmp(20)

C=====================================================================

      include 'const.dek'

      Pelposi=pe(T11,etae)*mev4toecc  !check2 in erg/cc

      

      return
      end

