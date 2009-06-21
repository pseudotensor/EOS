C======================================================================
c      Program electron_number_density   !check1
      real*8 function ne(T11,etae)

C======================================================================
C     This function computes electron number density, as functions of
C     temperature and the degeneracy parameter of electrons
C     /04/10/3/   K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist_num/etae_num,metilde_num

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
c      real*8 ne                 !check1

C======================================================================
      real*8 nelectron          !number deinsty of e- in MeV^3
      real*8 npositron          !number density of e+ in MeV^3

C======================================================================
      real*8 memev              !electron mass in MeV
      parameter(memev=0.511d0)
      real*8 metilde_num           !me/T

C======================================================================
      real*8 Qmev                  !qvalue of weak interaction between N,P
      parameter(Qmev=1.29d0)
      real*8 Qtilde

C======================================================================
      real*8 mev2K              ! 1MeV = mev2K  K
      parameter(mev2K=1.1605d10) ! 1MeV = mev2K  K
      real*8 PI
      parameter(PI=3.14159265d0)

C=====================================================================
      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae_num              !electron degeneracy parameter \mue/T

C======================================================================
      external distele_num
      external distposi_num

C======================================================================
c      read(*,*) T11,etae        !check1

      etae_num=etae

      metilde_num=memev/tmev
      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde_num=memev/tmev

      if (etae_num.le.metilde_num) then
         dum(31)=1.d0
      else
         dum(31)=Sqrt(etae_num**2-metilde_num**2)   ! dum(31) < 100.d0
      end if

      dum(1)=xintd(0.d0,1.d-3,distele_num,iter)   !check1
      dum(23)=xintd(1.d-3,1.d-2,distele_num,iter)   !check1
      dum(22)=xintd(1.d-2,1.d-1,distele_num,iter)   !check1
      dum(21)=xintd(1.d-1,1.d0,distele_num,iter)   !check1
c      dum(1)=xintd(0.d0,dum(31),distele_num,iter)
      dum(2)=xintd(1.d0,10.d0,distele_num,iter)
c      dum(2)=xintd(dum(31),100.d0,distele_num,iter)
      dum(3)=xintd(10.d0,100.d0,distele_num,iter)
      dum(4)=xintd(100.d0,1000.d0,distele_num,iter)
      dum(5) = dum(1)+dum(21)+dum(22)+dum(23) + dum(2) + dum(3) + dum(4)
      nelectron=tmev**3/PI**2 * dum(5)

      dum(11)=xintd(0.d0,1.d0,distposi_num,iter)
      dum(12)=xintd(1.d0,10.d0,distposi_num,iter)
      dum(13)=xintd(10.d0,100.d0,distposi_num,iter)
      dum(14)=xintd(100.d0,1000.d0,distposi_num,iter)
      dum(15) = dum(11) + dum(12) + dum(13) + dum(14)
      npositron=tmev**3/PI**2 * dum(15)

      ne = nelectron - npositron !in MeV^3

c      write(6,100) tmev, etae, ne/nelectron, ne*mev3tocc   !check1
 100  format(20(1pe11.3,' '))

c      stop                      !check1
      return   
      end


C======================================================================
      REAL*8 FUNCTION eval(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist_num/etae_num,metilde_num

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae_num              !electron degeneracy parameter

C======================================================================
      real*8 distele_num         !integrand of electrons distribution function
      real*8 distposi_num           !integrand of positrons distribution function


C======================================================================
      real*8 memev              !electron mass in MeV
      parameter(memev=0.511d0)
      real*8 metilde_num           !me/T

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
      ENTRY distele_num(x)
      
      tmp(1)= x**2
      tmp(2)=ex(Sqrt(x**2+metilde_num**2) - etae_num) + 1.d0
      distele_num=tmp(1) / tmp(2)

      RETURN

C----------------------------------------------------------------------
      ENTRY distposi_num(x)

      tmp(11)= x**2
      tmp(12)=ex(Sqrt(x**2+metilde_num**2) + etae_num) + 1.d0
      distposi_num=tmp(11) / tmp(12)

      RETURN

      END


C======================IDENTIFICATION DIVISION============================

      real*8 FUNCTION ex(x)

C----------LINKAGES.
C     CALLED BY - [subroutine] start, rate2, rate3, rate4, sol
C               - [function] eval
C     CALLS     - none

C----------REMARKS.
C     Exponential function with underflow precaution.

      implicit none

c-----------------local variable
      real*8 x

C=========================PROCEDURE DIVISION==============================

      IF (x.gt.88.0290D0) THEN        !In danger of overflow.
         ex = dexp(88.0290D0)
      ELSE
        IF (x.lt.-88.7220D0) THEN     !In danger of underflow.
          ex = 0.d0
        ELSE                       !Value of x in allowable range.
          ex = dexp(x)
        END IF
      END IF
      RETURN       

C----------NOTE-------------------------------------------------------------
C     The overflow limit for the VAX/VMS system is exp(88.029).
C     The underflow limit for the VAX/VMS system is exp(-88.722).

      END

