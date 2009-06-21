C======================================================================
c     Program cooling_rate_p_electron !check2
      real*8 function qne2(T11,etae,etanu,extrablock)  !in erg/sec

C======================================================================
C     This function computes the cooling rate through p + electron --> n
C     + nue, per proton number density as functions of temperature and
C     the degeneracy parameter of electrons 
C     /04/10/4/ K. Kohri
C======================================================================
      implicit none

C=====================================================================
C=====================================================================
C----------PARAMTER.
C======================================================================
      real*8 T11,etae,etanu,extrablock
C======================================================================
      real*8 dum(100)            !dummy parameters
C======================================================================
      integer ratepow
C======================================================================
      real*8 gammane_gen ! function

C======================================================================
c      read(*,*) T11,etae        !check2

      include 'const.dek'

      ratepow=3
c     etanu = etanuebar here
      qne2 = gammane_gen(ratepow,T11,etae,etanu,extrablock)


      return   !check2
      end


