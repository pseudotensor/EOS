C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammape(T11,etae,etanu,extrablock)

      implicit none

C----------PARAMTER.
C======================================================================
      real*8 T11,etae,etanu,extrablock
C======================================================================
      real*8 dum(100)            !dummy parameters
C======================================================================
      integer ratepow
C======================================================================
      real*8 gammape_gen ! function
C======================================================================


      include 'const.dek'

      ratepow=2
      gammape = gammape_gen(ratepow,T11,etae,etanu,extrablock)

      return   !check2
      end






C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammape_gen(ratepow,T11,etae,etanu,extrablock)

C======================================================================
C     This function computes the reaction rate through p + electron --> n
C     + nue , as functions of temperature and the degeneracy parameter
C     of electrons 
C     /04/10/4/ K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist_gammape/etae_gpe,met_gpe,Qt_gpe,etanu_gpe,ratepow_gpe,extrablock_gpe

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
c      real*8 gammape                !number deinsty of e- in MeV^3 !check2
      real*8 npositron          !number density of e+ in MeV^3

C======================================================================
c      real*8 memev              !electron mass in MeV
c      parameter(memev=0.511d0)
c      parameter(memev=0.d0)  !check1
      real*8 met_gpe           !me/T

C======================================================================
      real*8 Qt_gpe
C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae_gpe,etanu_gpe,ratepow_gpe              !electron degeneracy parameter \mue/T
      real*8 extrablock_gpe

C======================================================================
      real*8 np                 !proton number density in MeV^3

C======================================================================
      real*8 neint                 !function of electron number density in MeV^3

C======================================================================
      external integgammape

C======================================================================
c      read(*,*) T11,etae        !check2

      include 'const.dek'


      etae_gpe=etae
      etanu_gpe=etanu
      extrablock_gpe=extrablock
      ratepow_gpe=dble(ratepow)

      tmev = T11*1.d11/mev2K
      Qt_gpe=Qmev/tmev
      met_gpe=memev/tmev

      dum(5) = distintegrate(integgammape)

      dum(10)=GF**2/(2.d0*PI**3)*(1.d0+3.d0*ga**2)*tmev**5 !in MeV
      if(ratepow.eq.2) then
         dum(11) =  dum(10) * mevtosec ! in 1/sec
      else if(ratepow.eq.3) then
         dum(11) =  dum(10) * tmev * mev2toergs ! in erg/sec
      end if

      gammape_gen = dum(11) * dum(5)

c      write(6,100) tmev, etae, gammape  !check2
 100  format(20(1pe11.3,' '))

c      stop
      return   !check2
      end



C======================================================================
      real*8 function integgammape(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist_gammape/etae_gpe,met_gpe,Qt_gpe,etanu_gpe,ratepow_gpe,extrablock_gpe

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae_gpe,etanu_gpe,ratepow_gpe    !electron degeneracy parameter
      real*8 extrablock_gpe
      real*8 met_gpe           !me/T
      real*8 Qt_gpe             !Q/T

C======================================================================
c     real*8 integgammape       !integrand of electrons distribution function

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

c ratepow = 2 for number process and = 3 for energy process

C======================================================================
C----------------------------------------------------------------------
C     JCM: Note: x=pc/(k_b T)
      tmp(12)=ddim( (x+Qt_gpe)**2-met_gpe**2, 0.d0 )
      tmp(1)= (x+Qt_gpe)*dsqrt(tmp(12)) * x**(ratepow_gpe)
      tmp(2)=dexp(x+Qt_gpe - etae_gpe) + 1.d0
c      tmp(3)=1.d0/(1.d0 + exp(-x))

c     Pauli blocking factor
      tmp(3) = 1.0d0 - 1.0d0*extrablock_gpe/(dexp(x-etanu_gpe)+1.0d0)

      integgammape=tmp(3) * tmp(1) / tmp(2)
c      integgammape=tmp(1) / tmp(2) *tmp(3) !check

      RETURN

      END


