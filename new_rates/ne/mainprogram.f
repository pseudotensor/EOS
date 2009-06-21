C======================================================================
      Program data_taking

C======================================================================
      implicit none

C======================================================================
      integer i,j,k

C======================================================================
      real*8 T11                !temperature in 10^11 K
      integer nT
      parameter(nT=701)
c      parameter(nT=71)
      real*8 intT
      real*8 Tmax
      parameter(Tmax=1.d2)
      real*8 Tmin
      parameter(Tmin=1.d-5)
      
C======================================================================
      real*8 etae               !electron degeneracy parameter \mue/T
      integer neta
      parameter(neta=501)
c      parameter(neta=101)
      real*8 inteta
      real*8 etamax
      parameter(etamax=1.d2)
      real*8 etamin
      parameter(etamin=1.d-3)
c      parameter(etamin=1.d-40)

C=====================================================================
      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C=====================================================================
      real*8 ne                 !function

C=====================================================================
      real sint(nT,neta)

C======================================================================
      real*8 tmp(20)

C======================================================================
      open(1,file="ne1.dat",status="unknown")
      open(2,file="imne",form='unformatted',status="unknown")

      intT=log10(Tmax/Tmin)/real(nT-1)
      inteta=log10(etamax/etamin)/real(neta-1)
      do i=1,nT
         T11=Tmin*10.d0**(intT*real(i-1))
         do j=1,neta
            etae=etamin*10.d0**(inteta*real(j-1))
            tmp(1)=ne(T11,etae)
            write(1,100) T11, etae, tmp(1)
            sint(i,j)=sngl( log10(tmp(1))+log10(mev3tocc) )
c            sint(i,j)=sngl( log10(tmp(1)) )
            
 100        format(20(1pe11.3,' '))
         end do
      end do

      write(2) nT,neta
      write(2) sint
      close(2)

      close(1)

      stop
      end


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

