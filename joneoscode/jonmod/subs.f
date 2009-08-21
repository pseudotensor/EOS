







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
c      DATA u/-.93246951420315,-.66120938646627,-.23861918608320,
c     |        .23861918608320, .66120938646627, .93246951420315/  
c      DATA w/.17132449237917,.36076157304814,.46791393457269,
c     |       .46791393457269,.36076157304814,.17132449237917/        
      DATA np/6/              !6 point Gaussian integration.

      DATA u/-0.9324695142031520278123015545,
     |     -0.661209386466264513661399595,
     |     -0.238619186083196908630501722,
     |     0.238619186083196908630501722,
     |     0.661209386466264513661399595,
     |     0.932469514203152027812301554/

      DATA w/0.17132449237917034504029614217,
     |     0.36076157304813860756983351384,
     |     0.46791393457269104738987034399,
     |     0.46791393457269104738987034399,
     |     0.36076157304813860756983351384,
     |     0.17132449237917034504029614217/

C=========================PROCEDURE DIVISION==============================

C10--------DO INTEGRATION----------------------------------------------------

      sum   = 0.d0       
c      dist  = (xhi-xlow)/float(nq) !Size of quad interval.
      dist  = (xhi-xlow)/dble(nq) !Size of quad interval.
      DO nint = 1,nq
c        cent = xlow+(float(nint)-0.5)*dist  !Center of interval.
        cent = xlow+(dble(nint)-0.5)*dist  !Center of interval.
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

      ex1 = dexp(x)
      return


      IF (x.gt.600.0d0) THEN        !In danger of overflow.
c       IF (x.gt.1.d3) THEN        !In danger of overflow.  !check2
          ex1 = dexp(600.0d0)  
c          ex1 = dexp(1.d3)       !check2
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

      ex = dexp(x)
      return


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

