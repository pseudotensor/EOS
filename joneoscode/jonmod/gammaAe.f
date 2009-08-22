C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammaAe(T11,etae,etap,etan,etanu,abar,zbar,extrablock)

      implicit none

C----------PARAMTER.
C======================================================================
      real*8 T11,etae,etap,etan,etanu,abar,zbar,extrablock
C======================================================================
      real*8 dum(100)            !dummy parameters
C======================================================================
      integer ratepow
C======================================================================
      real*8 gammaAe_gen ! function
C======================================================================


      include 'const.dek'

      ratepow=2
      gammaAe = gammaAe_gen(ratepow,T11,etae,etap,etan,etanu,abar,zbar,extrablock)

      return   !check2
      end




C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammaAe_gen(ratepow,T11,etae,etap,etan,etanu,abar,zbar,extrablock)

C======================================================================
C     This function computes the reaction rate through A + electron --> n
C     + nue , as functions of temperature and the degeneracy parameter
C     of electrons 
C     /04/10/4/ K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist_gammaAe/etae_gAe,met_gAe,Qt_gAe,etanu_gAe,ratepow_gAe,extrablock_gAe

C=====================================================================
C----------PARAMTER.
      integer   iter
      PARAMETER (iter=50)          !Number of gaussian quads.
      real*8 xintd              !function
      real*8 distintegrate      !function
C======================================================================
      real*8 T11                !temperature in 10^11 K
      real*8 mue11              !electron chemical potential in 10^11 K
      real*8 muemev             !electron chemical potential in MeV

C======================================================================
      real*8 tmev               !temperature in MeV
      real*8 etae,etap,etan,etanu     !electron degeneracy parameter \mue/T
      real*8 abar,zbar
      real*8 extrablock
      integer ratepow
C======================================================================
c      real*8 ne                 !check1
C======================================================================
      real*8 Np,Nh
C======================================================================
c      real*8 gammape                !number deinsty of e- in MeV^3 !check2
c      real*8 npositron          !number density of e+ in MeV^3

C======================================================================
c      real*8 memev              !electron mass in MeV
c      parameter(memev=0.511d0)
c      parameter(memev=0.d0)  !check1
      real*8 met_gAe           !me/T

C======================================================================
c      real*8 Qmev                  !qvalue of weak interaction between N,P
c      parameter(Qmev=1.29d0)
      real*8 Qt_gAe

C======================================================================
c      real*8 mev2K              ! 1MeV = mev2K  K
c      parameter(mev2K=1.1605d10) ! 1MeV = mev2K  K
c      real*8 PI
c      parameter(PI=3.14159265d0)
c      real*8 mevtosec           !1 MeV = (mevtosec) /sec
c      parameter(mevtosec=1.519d21) !1 MeV = (mevtosec) /sec

C======================================================================
c      real*8 ga                 !axial vector coupling
c      parameter(ga=1.39d0)
c      real*8 GF                 !Fermi coupling constant
c      parameter(GF=1.1664d-11)  !in MeV^-2

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae_gAe,etanu_gAe   !electron degeneracy parameter \mue/T
      real*8 extrablock_gAe
      real*8 ratepow_gAe
      real*8 nbar
C======================================================================
c      real*8 np                 !proton number density in MeV^3

C======================================================================
c      real*8 neint                 !function of electron number density in MeV^3

C======================================================================
      external integgammaAe

C======================================================================
c      read(*,*) T11,etae        !check2

      include 'const.dek'


      etae_gAe=etae
      etanu_gAe=etanu
      ratepow_gAe=dble(ratepow)
      extrablock_gAe=extrablock

      tmev = T11*1.d11/mev2K
c      JCM: Qt_gAe term different than p+e^- -> n + \nu
c     From Bruenn (1985) equation C28
      Qt_gAe = dmax1(Deltamev/tmev + etan - etap,0.0d0)
c     use above so generally truncated
c      Qt_gAe = Deltamev/tmev + etan - etap
      met_gAe=memev/tmev

      dum(5) = distintegrate(integgammaAe)

c      JCM: ga term different than p+e^- -> n + \nu
      dum(10)=GF**2/(2.d0*PI**3)*((2.0/7.0)*ga**2)*tmev**5 !in MeV
      if(ratepow.eq.2) then
         dum(11) =  dum(10) * mevtosec ! in 1/sec
      else if(ratepow.eq.3) then
         dum(11) =  dum(10) * tmev * mev2toergs ! in erg/sec
      end if


c     Kawanaka & Mineshige (2007)
      if(zbar.lt.20.0) then
         Np = 0.0
      else if(zbar.ge.20.0 .AND. zbar.lt.28.0) then
         Np = zbar - 20.0
      else if(zbar.ge.28.0) then
         Np = 8.0
      end if

      nbar = abar-zbar
      if(nbar.lt.34.0) then
         Nh = 6.0
      else if(nbar.ge.34.0 .AND. nbar.lt.40.0) then
         Nh = 40.0 - nbar
      else if(nbar.ge.40.0) then
         Nh = 0.0
      end if
 
c     DEBUG:
c      write(*,*) 'nbar',abar,zbar,nbar,Np,Nh


      gammaAe_gen = dum(11) * dum(5) * Np*Nh


c      write(6,100) tmev, etae, gammape  !check2
 100  format(20(1pe11.3,' '))

c      stop
      return   !check2
      end

C======================================================================
      real*8 function integgammaAe(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist_gammape/etae_gAe,met_gAe,Qt_gAe,etanu_gAe,ratepow_gAe,extrablock_gAe

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae_gAe,etanu_gAe              !electron degeneracy parameter
      real*8 extrablock_gAe
      real*8 met_gAe           !me/T
      real*8 Qt_gAe             !Q/T
      real*8 ratepow_gAe
C======================================================================
c     real*8 integgammape       !integrand of electrons distribution function

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
C     JCM: Note x = (E_e - Q')/(k_b T) and Q'=Qt/(k_B T)
      tmp(12)=ddim( (x+Qt_gAe)**2-met_gAe**2, 0.d0 )
      tmp(1)= (x+Qt_gAe)*dsqrt(tmp(12)) * x**(ratepow_gAe)
      tmp(2)=dexp(x+Qt_gAe - etae_gAe) + 1.d0
c      tmp(3)=1.d0/(1.d0 + exp(-x))

c     Pauli blocking factor
      tmp(3) = 1.0d0 - 1.0d0*extrablock_gAe/(dexp(x-etanu_gAe)+1.0d0)

      integgammaAe=tmp(3) * tmp(1) / tmp(2)
c      integgammaAe=tmp(1) / tmp(2) *tmp(3) !check

      RETURN

      END


