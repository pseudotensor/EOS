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
      common /dist1/etae1,metilde

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

C=====================================================================
      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae1              !electron degeneracy parameter \mue/T

C======================================================================
      external distele
      external distposi

C======================================================================
c      read(*,*) T11,etae        !check1

      etae1=etae

      metilde=memev/tmev
      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde=memev/tmev

      if (etae1.le.metilde) then
         dum(31)=1.d0
      else
         dum(31)=Sqrt(etae1**2-metilde**2)   ! dum(31) < 100.d0
      end if

c      dum(1)=xintd(0.d0,1.d0,distele,iter)   !check1
      dum(1)=xintd(0.d0,dum(31),distele,iter)
c      dum(2)=xintd(1.d0,10.d0,distele,iter)
      dum(2)=xintd(dum(31),100.d0,distele,iter)
c      dum(3)=xintd(10.d0,100.d0,distele,iter)
      dum(3)=0.d0
      dum(4)=xintd(100.d0,1000.d0,distele,iter)
      dum(5) = dum(1) + dum(2) + dum(3) + dum(4)
      nelectron=tmev**3/PI**2 * dum(5)

      dum(11)=xintd(0.d0,1.d0,distposi,iter)
      dum(12)=xintd(1.d0,10.d0,distposi,iter)
      dum(13)=xintd(10.d0,100.d0,distposi,iter)
      dum(14)=xintd(100.d0,1000.d0,distposi,iter)
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
      common /dist1/etae1,metilde

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae1              !electron degeneracy parameter

C======================================================================
      real*8 distele            !integrand of electrons distribution function
      real*8 distposi           !integrand of positrons distribution function


C======================================================================
      real*8 metilde           !me/T

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function
      real*8 ex1                 !function

C======================================================================
C----------------------------------------------------------------------
      ENTRY distele(x)
      
      tmp(1)= x**4/Sqrt(x**2+metilde**2)
      tmp(2)=ex1(Sqrt(x**2+metilde**2) - etae1) + 1.d0
      distele=tmp(1) / tmp(2)

      RETURN

C----------------------------------------------------------------------
      ENTRY distposi(x)

      tmp(11)= x**4/Sqrt(x**2+metilde**2)
      tmp(12)=ex1(Sqrt(x**2+metilde**2) + etae1) + 1.d0
      distposi=tmp(11) / tmp(12)

      RETURN

      END


C======================IDENTIFICATION DIVISION============================

      real*8 FUNCTION xintd (xlow,xhi,func,nq)

C----------LINKAGES.
C     CALLED BY - [subroutine] rate1, nudens
C     CALLS     - none

C----------REMARKS.
C     Computes the integral of the function "func".

      implicit none

C========================DECLARATION DIVISION=============================

C----------INPUT VARIABLES.
      REAL*8    xlow                 !Array of low limits.
      REAL*8    xhi                  !Array of high limits.
      INTEGER nq                   !Number of six point gaussian quads.

C----------COMPUTATION VARIABLES.
      REAL*8    dist                 !Size of quad interval.
      REAL*8    cent                 !Center of quad interval.
      REAL*8    x                    !Variables of integration.
      REAL*8    sum                  !Summation of terms.
      real*8    f

C----------COUNTERS.
      INTEGER nint                 !Interval number.
      INTEGER npnt                 !Point number.
      INTEGER np                   !Total number of points in interval.

C----------ABSCISSAS AND WEIGHT FACTORS.
      REAL*8    u(6)                 !Abscissas.
      REAL*8    w(6)                 !Weight factor.

C----------function type
      REAL*8    func
c      external    func

C============================DATA DIVISION================================

C----------ABSCISSAS AND WEIGHT FACTORS.
      DATA u/-.93246951420315,-.66120938646627,-.23861918608320,
     |        .23861918608320, .66120938646627, .93246951420315/  
      DATA w/.17132449237917,.36076157304814,.46791393457269,
     |       .46791393457269,.36076157304814,.17132449237917/        
      DATA np/6/              !6 point Gaussian integration.


C=========================PROCEDURE DIVISION==============================

C10--------DO INTEGRATION----------------------------------------------------

      sum   = 0.d0       
      dist  = (xhi-xlow)/float(nq) !Size of quad interval.
      DO nint = 1,nq
        cent = xlow+(float(nint)-0.5)*dist  !Center of interval.
        DO npnt = 1,np
          x   = cent+0.5*dist*u(npnt) !Integration point.
          f   = func(x)            !Evaluate function x(1).
          sum = sum+f*w(npnt)      !Add up sum.
        END DO
      END DO

C20--------GET INTEGRAL VALUE------------------------------------------------

      xintd = sum*dist*0.5         !Do integral.
      RETURN        

      END 

C======================IDENTIFICATION DIVISION============================

      real*8 FUNCTION ex1(x)

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

c      IF (x.gt.88.0290D0) THEN        !In danger of overflow.
       IF (x.gt.2.d3) THEN        !In danger of overflow.  !check2
c          ex1 = dexp(88.0290D0)  
          ex1 = dexp(2.d3)       !check2
      ELSE
        IF (x.lt.-88.7220D0) THEN     !In danger of underflow.
          ex1 = 0.d0
        ELSE                       !Value of x in allowable range.
          ex1 = dexp(x)
        END IF
      END IF
      RETURN       

C----------NOTE-------------------------------------------------------------
C     The overflow limit for the VAX/VMS system is exp(88.029).
C     The underflow limit for the VAX/VMS system is exp(-88.722).

      END

