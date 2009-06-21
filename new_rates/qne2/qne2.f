C======================================================================
c      Program cooling_rate_n_positron  !check2
      real*8 function qne2(T11,etae)  !in erg/sec 

C======================================================================
C     This function computes the cooling rate through p + electron --> n
C     + nue per neutron number density, as functions of temperature and
C     the degeneracy parameter of electrons 
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
c      real*8 qne2          !number deinsty of e- in MeV^3  !check2
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
c      real*8 mev5toergcms       !1 MeV^5 = (mev5toergcms) erg /cc/sec
c      parameter(mev5toergcms=3.1678d47) !1 MeV^5 = (mev5toergcms) erg /cc/sec
      real*8 mev2toergs         ! 1MeV^2 = (mev2toergs) erg/sec
      parameter(mev2toergs=2.4341d15)


C======================================================================
      real*8 ga                 !axial vector coupling
      parameter(ga=1.39d0)
      real*8 GF                 !Fermi coupling constant
      parameter(GF=1.1664d-11)  !in MeV^-2

C======================================================================
      real*8 dum(20)            !dummy parameters
      real*8 etae1              !electron degeneracy parameter \mue/T

C======================================================================
      real*8 nn                 !proton number density in MeV^3

C======================================================================
      real*8 neint                 !function of electron number density in MeV^3
      real*8 ex1 

C======================================================================
      external integqne2

C======================================================================
c      write(6,*) "temperature in 1.e11 K = ", "eta_e ="
c      read(*,*) T11,etae        !check1

      etae1=etae

      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde=memev/tmev

      dum(1)=xintd(0.d0,1.d0,integqne2,iter)
      dum(2)=xintd(1.d0,10.d0,integqne2,iter)
      dum(3)=xintd(10.d0,100.d0,integqne2,iter)
      dum(4)=xintd(100.d0,1000.d0,integqne2,iter)
      dum(5) = dum(1) + dum(2) + dum(3) + dum(4)

c      nn=neint(T11,etae)*ex1(etae-Qtilde)
      nn=1.d0                   !check check
c      nn=2.313e1                   !check1

      dum(10)=GF**2/(2.d0*PI**3)*(1.d0+3.d0*ga**2)*nn*tmev**6  !in MeV^5
c      dum(11) =  dum(10) * mev5toergcms ! in erg/cc/sec
      dum(11) =  dum(10) * mev2toergs ! in erg/sec
      qne2= dum(11) * dum(5)

c      write(6,100) tmev, etae, qne2
 100  format(20(1pe11.3,' '))

c      stop
      return                    !check2
      end

C======================================================================
      real*8 function integqne2(x)

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
c      real*8 integqne2            !integrand of electrons distribution function

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex1                 !function

C======================================================================
C----------------------------------------------------------------------
      tmp(1)= (x+metilde)*Sqrt((x+2.d0*metilde)*x) * (x+metilde+Qtilde)**3
      tmp(2)=ex1(x+Qtilde + etae1) + 1.d0
c      tmp(2)=exp(x+Qtilde + etae1) + 1.d0
      integqne2=tmp(1) / tmp(2)

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

