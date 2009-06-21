C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammapnu(T11,etae)

C======================================================================
C     This function computes the reaction rate through p + electron --> n
C     + nue , as functions of temperature and the degeneracy parameter
C     of electrons 
C     /04/10/4/ K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist1/etae1,metilde,Qtilde

C=====================================================================
C----------PARAMTER.
      integer   iter
      PARAMETER (iter=50)          !Number of gaussian quads.
      real*8 xintd              !function

C======================================================================
      real*8 T11                !temperature in 10^11 K
      real*8 mue11              !electron chemical potential in 10^11 K
      real*8 muemev             !electron chemical potential in MeV

C======================================================================
      real*8 tmev               !temperature in MeV
      real*8 etae               !electron degeneracy parameter \mue/T

C======================================================================
      real*8 ne                 !check1

C======================================================================
c      real*8 gammapnu                !number deinsty of e- in MeV^3 !check2
      real*8 npositron          !number density of e+ in MeV^3

C======================================================================
      real*8 memev              !electron mass in MeV
      parameter(memev=0.511d0)
c      parameter(memev=0.d0)  !check1
      real*8 metilde           !me/T

C======================================================================
      real*8 Qmev                  !qvalue of weak interaction between N,P
      parameter(Qmev=1.29d0)
      real*8 Qtilde

C======================================================================
      real*8 mev2K              ! 1MeV = mev2K  K
      parameter(mev2K=1.1605d10) ! 1MeV = mev2K  K
      real*8 PI
      parameter(PI=3.14159265d0)
      real*8 mevtosec           !1 MeV = (mevtosec) /sec
      parameter(mevtosec=1.519d21) !1 MeV = (mevtosec) /sec

C======================================================================
      real*8 ga                 !axial vector coupling
      parameter(ga=1.39d0)
      real*8 GF                 !Fermi coupling constant
      parameter(GF=1.1664d-11)  !in MeV^-2

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae1              !electron degeneracy parameter \mue/T

C======================================================================
      real*8 np                 !proton number density in MeV^3

C======================================================================
      real*8 neint                 !function of electron number density in MeV^3

C======================================================================
      external integgammapnu

C======================================================================
c      read(*,*) T11,etae        !check2

      etae1=etae

      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde=memev/tmev

      dum(1)=xintd(0.d0,1.d0,integgammapnu,iter)
      dum(2)=xintd(1.d0,10.d0,integgammapnu,iter)
      dum(3)=xintd(10.d0,100.d0,integgammapnu,iter)
      dum(4)=xintd(100.d0,1000.d0,integgammapnu,iter)
      dum(5) = dum(1) + dum(2) + dum(3) + dum(4)


      dum(10)=GF**2/(2.d0*PI**3)*(1.d0+3.d0*ga**2)*tmev**5 !in MeV
      dum(11) =  dum(10) * mevtosec ! in 1/sec
      gammapnu= dum(11) * dum(5)

c      write(6,100) tmev, etae, gammapnu  !check2
 100  format(20(1pe11.3,' '))

c      stop
      return   !check2
      end

C======================================================================
      real*8 function integgammapnu(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist1/etae1,metilde,Qtilde

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae1              !electron degeneracy parameter
      real*8 metilde           !me/T
      real*8 Qtilde             !Q/T

C======================================================================
c     real*8 integgammapnu       !integrand of electrons distribution function

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
      tmp(1)= (x+metilde)*Sqrt(x*(x+2.d0*metilde))*(x+metilde+Qtilde)**2
      tmp(2)=exp(x+metilde+Qtilde) + 1.d0
      tmp(3)=1.d0/(  1.d0 + exp(-(x+metilde+etae1)))
c      integgammapnu=tmp(1) / tmp(2)
      integgammapnu=tmp(1) / tmp(2)*tmp(3)

      RETURN

      END


