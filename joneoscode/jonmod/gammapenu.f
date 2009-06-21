C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammapenu(T11,etae,etanu,extrablock)

      implicit none

C----------PARAMTER.
C======================================================================
      real*8 T11,etae,etanu,extrablock
C======================================================================
      real*8 dum(100)            !dummy parameters
C======================================================================
      integer ratepow
C======================================================================
      real*8 gammapenu_gen ! function
C======================================================================


      include 'const.dek'

      ratepow=2
      gammapenu = gammapenu_gen(ratepow,T11,etae,etanu,extrablock)

      return   !check2
      end



C======================================================================
c      Program reaction_rate_p_electron  !check2
c GODMARK: This rate doesn't account for neutrino energy
c
c
c
      real*8 function gammapenu_gen(ratepow,T11,etae,etanu,extrablock)

C======================================================================
C     This function computes the reaction rate through p + electron --> n
C     + nue , as functions of temperature and the degeneracy parameter
C     of electrons 
C     /04/10/4/ K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist1/etae1,metilde,Qtilde,etanu1,ratepow1,extrablock1

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
      real*8 etae,etanu               !electron degeneracy parameter \mue/T
      real*8 extrablock
      integer ratepow
C======================================================================
      real*8 ne                 !check1

C======================================================================
c      real*8 gammapenu                !number deinsty of e- in MeV^3 !check2
      real*8 npositron          !number density of e+ in MeV^3

C======================================================================
c      real*8 memev              !electron mass in MeV
c      parameter(memev=0.511d0)
c      parameter(memev=0.d0)  !check1
      real*8 metilde           !me/T

C======================================================================
c      real*8 Qmev                  !qvalue of weak interaction between N,P
c      parameter(Qmev=1.29d0)
      real*8 Qtilde
      real*8 ratepow1
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
      real*8 etae1,etanu1              !electron degeneracy parameter \mue/T
      real*8 extrablock1

C======================================================================
      real*8 np                 !proton number density in MeV^3

C======================================================================
      real*8 neint                 !function of electron number density in MeV^3

C======================================================================
      external integgammapenu

C======================================================================
c      read(*,*) T11,etae        !check2


      include 'const.dek'


      etae1=etae
      etanu1=etanu
      extrablock1=extrablock
      ratepow1=dble(ratepow)

      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde=memev/tmev

      dum(50)=Qtilde-metilde
c     JCM: Notice the limited range of integration
c
c JCM: Is the below sufficient?  GODMARK:
      dum(1)=xintd(0.d0,dum(50),integgammapenu,iter)
      dum(2)=0.d0
      dum(3)=0.d0
      dum(4)=0.d0
      dum(5) = dum(1) + dum(2) + dum(3) + dum(4)
c JCM (below incorrect since goes above range of integration)
c      dum(5) = distintegrate(integgammapenu)



      dum(10)=GF**2/(2.d0*PI**3)*(1.d0+3.d0*ga**2)*tmev**5 !in MeV
      if(ratepow.eq.2) then
         dum(11) =  dum(10) * mevtosec ! in 1/sec
      else if(ratepow.eq.3) then
         dum(11) =  dum(10) * tmev * mev2toergs ! in erg/sec
      end if

      gammapenu_gen = dum(11) * dum(5)

c      write(6,100) tmev, etae, gammapenu  !check2
 100  format(20(1pe11.3,' '))

c      stop
      return   !check2
      end



C======================================================================
      real*8 function integgammapenu(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist1/etae1,metilde,Qtilde,etanu1,ratepow1,extrablock1

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae1,etanu1              !electron degeneracy parameter
      real*8 extrablock1
      real*8 metilde           !me/T
      real*8 Qtilde             !Q/T
      real*8 ratepow1
C======================================================================
c     real*8 integgammapenu       !integrand of electrons distribution function

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
C     JCM: Note: x=pc/(k_b T)
      tmp(12) = ddim((x-Qtilde)**2 -metilde**2,0.0)
      tmp(1)= dsqrt( tmp(12)) * (Qtilde-x)* x**(ratepow1)
      tmp(2)=dexp( -(x-Qtilde+etae1) ) + 1.d0

c     Pauli Blocking Factor
      tmp(3)=1.0d0*extrablock1/(  1.d0 + dexp(x-etanu1) )

      integgammapenu=tmp(1) / tmp(2)*tmp(3)

      RETURN

      END


