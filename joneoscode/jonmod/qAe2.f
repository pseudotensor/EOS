C======================================================================
c     Program cooling_rate_p_electron !check2
      real*8 function qAe2(T11,etae,etap,etan,etanu,abar,zbar,extrablock)  !in erg/sec
C======================================================================
      implicit none
C======================================================================
c     Heavy nuclei electron capture rate
c     Accurately ignores alphas and other lower A/Z's
C=====================================================================
C=====================================================================
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
c      read(*,*) T11,etae        !check2

      include 'const.dek'

      ratepow=3
      qAe2 = gammaAe_gen(ratepow,T11,etae,etap,etan,etanu,abar,zbar,extrablock)


      return   !check2
      end




