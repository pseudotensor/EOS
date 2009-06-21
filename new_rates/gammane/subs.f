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
