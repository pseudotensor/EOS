C======================================================================
c      Program electron_number_density   !check1 in MeV^4
      real*8 function pe(T11,etae) !in MeV^4

C======================================================================
C     This function computes electron pressure, as functions of
C     temperature and the degeneracy parameter of electrons /04/10/3/
C     K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist_pres/etae_pres,metilde_pres

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
      real*8 etae               !electron degeneracy parameter \mue/T

C======================================================================
c      real*8 pe                 !check1 in MeV^4

C======================================================================
      real*8 pelectron          !pressure of e- in MeV^4
      real*8 ppositron          !pressure of e+ in MeV^4

C======================================================================
c      real*8 memev              !electron mass in MeV
c      parameter(memev=0.511d0) 
c      parameter(memev=0.d0) 
      real*8 metilde_pres           !me/T

C======================================================================
c      real*8 Qmev                  !qvalue of weak interaction between N,P
c      parameter(Qmev=1.29d0)
c      real*8 Qtilde

C======================================================================
c      real*8 mev2K              ! 1MeV = mev2K  K
c      parameter(mev2K=1.1605d10) ! 1MeV = mev2K  K
c      real*8 PI
c      parameter(PI=3.14159265d0)

C=====================================================================
c      real*8 mev4toecc           !1 MeV^4 = (mev4toecc) erg/cc
c      parameter(mev4toecc=2.085d26)

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 lowerdum(100)       !extra dummy parameters JCM
      real*8 etae_pres              !electron degeneracy parameter \mue/T

C======================================================================
      external distele_pres
      external distposi_pres

C======================================================================
c      read(*,*) T11,etae        !check1

      include 'const.dek'



      etae_pres=etae

      tmev = T11*1.d11/mev2K

c      Qtilde=Qmev/tmev
      metilde_pres=memev/tmev

      if (etae_pres.le.metilde_pres) then
         dum(31)=1.d0
      else
         dum(31)=Sqrt(etae_pres**2-metilde_pres**2)  ! dum(31) < 100.d0
      end if


      dum(5) = distintegrate(distele_pres)

      pelectron=1.d0/(3.d0*PI**2)*tmev**4 * dum(5)

      dum(15) = distintegrate(distposi_pres)

      ppositron=1.d0/(3.d0*PI**2)*tmev**4 * dum(15)

      pe = pelectron + ppositron !in MeV^4

c      write(6,100) tmev, etae, pelectron, ppositron, pe*mev4toecc   !check1
 100  format(20(1pe11.3,' '))

c      stop                      !check1
      return   
      end

C======================================================================
      REAL*8 FUNCTION eval_pre(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist_pres/etae_pres,metilde_pres

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae_pres              !electron degeneracy parameter

C======================================================================
      real*8 distele_pres            !integrand of electrons distribution function
      real*8 distposi_pres           !integrand of positrons distribution function


C======================================================================
      real*8 metilde_pres           !me/T

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function
      real*8 ex1                 !function

C======================================================================
C----------------------------------------------------------------------
      ENTRY distele_pres(x)
      
      tmp(1)= x**4/Sqrt(x**2+metilde_pres**2)
      tmp(2)=ex1(Sqrt(x**2+metilde_pres**2) - etae_pres) + 1.d0
      distele_pres=tmp(1) / tmp(2)

      RETURN

C----------------------------------------------------------------------
      ENTRY distposi_pres(x)

      tmp(11)= x**4/Sqrt(x**2+metilde_pres**2)
      tmp(12)=ex1(Sqrt(x**2+metilde_pres**2) + etae_pres) + 1.d0
      distposi_pres=tmp(11) / tmp(12)

      RETURN

      END



c   function to perform segmented integration of distribution functions
c   used for all code -- not just this file
C   ======================================================================
      real*8 function distintegrate(funtype)
C======================================================================
      implicit none

C=====================================================================
c      common /variable1/ tdyn,hcm,rhob,tk
C======================================================================
      real*8 funtype
C======================================================================
c      real*8 tdyn               !dynamical timescale in sec
c      real*8 hcm                !disk half thickness in cm
c      real*8 rhob               !baryon denstiy in g/cc
c      real*8 tk                 !temperature in K
C======================================================================
C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 lowerdum(100)       !extra dummy parameters JCM
C======================================================================
      integer   iter,guessiter
      PARAMETER (guessiter=50)          !Number of gaussian quads.
      real*8 xintd              !function

c      if(tk<1E8) then
c         iter=200
c      else
      iter=guessiter
c      end if

c      iter=1000


c     JCM
      lowerdum(1) = xintd(1.d-9,1.d-8,funtype,iter)
      lowerdum(2) = xintd(1.d-8,1.d-7,funtype,iter)
      lowerdum(3) = xintd(1.d-7,1.d-6,funtype,iter)
      lowerdum(4) = xintd(1.d-6,1.d-5,funtype,iter)
      lowerdum(5) = xintd(1.d-5,1.d-4,funtype,iter)
      dum(1)=xintd(1.d-4,1.d-3,funtype,iter)
      dum(21)=xintd(1.d-3,1.d-2,funtype,iter)
      dum(22)=xintd(1.d-2,1.d-1,funtype,iter)
      dum(23)=xintd(1.d-1,1.d0,funtype,iter)
      dum(2)=xintd(1.d0,10.d0,funtype,iter)
      dum(3)=xintd(10.d0,100.d0,funtype,iter)
      dum(4)=xintd(100.d0,1000.d0,funtype,iter)
      dum(24)=xintd(1.d3,1.d4,funtype,iter)
      dum(25)=xintd(1.d4,1.d5,funtype,iter)
      dum(26)=xintd(1.d5,1.d6,funtype,iter)
c     JCM
      dum(27)=xintd(1.d6,1.d7,funtype,iter)
      dum(28)=xintd(1.d7,1.d8,funtype,iter)
      dum(29)=xintd(1.d8,1.d9,funtype,iter)

      dum(5)=dum(1)+dum(2)+dum(3)+dum(4)
     |             +dum(21)+dum(22)+dum(23)
     |             +dum(24)+dum(25)+dum(26)
     |             +dum(27)+dum(28)+dum(29)
     |             +lowerdum(1)+lowerdum(2)+lowerdum(3)
     |             +lowerdum(4)+lowerdum(5)

      distintegrate = dum(5)

c DEBUG:
c      write(*,*) dum(25),dum(26),dum(27)

      return
      end
