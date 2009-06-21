C======================================================================
c      Program electron_number_density   !check1
      real*8 function ne(T11,etae)  !MeV^3

      implicit none

      real*8 T11                !temperature in 10^11 K
      real*8 etae               !electron degeneracy parameter \mue/T

      real*8 ne_jon
      real*8 ne_kaz

      ne=ne_kaz(T11,etae)
c      ne=ne_jon(T11,etae)


      return
      end



C======================================================================
c      Program electron_number_density   !check1
c JCM: Apparently my function is not correct/working
      real*8 function ne_jon(T11,etae)  !MeV^3

C======================================================================
C     This function computes electron number density, as functions of
C     temperature and the degeneracy parameter of electrons
C     /04/10/3/   K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist_num_jon/etae_num,beta_num

C=====================================================================
C----------PARAMTER.
      integer   iter
      PARAMETER (iter=50)          !Number of gaussian quads.
      real*8 xintd              !function
      real*8 distintegrate
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

C=====================================================================
c      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C======================================================================
      real*8 tk
      real*8 prefactor
C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae_num            !electron degeneracy parameter \mue/T
      real*8 beta_num
C======================================================================
      external distele_num_jon
      external distposi_num_jon

C======================================================================
c      read(*,*) T11,etae        !check1

      include 'const.dek'

      tk=T11*1D11
      etae_num=etae
      beta_num=(kerg*tk)/(me*clight**2)
      prefactor = (sqrt(2.0)/hbarcgs**3/pi**2*(kerg*tk)**2*me/clight)/mev3tocc ! in MeV^3


      dum(5) = distintegrate(distele_num_jon)
      nelectron= prefactor * dum(5)

      dum(15) = distintegrate(distposi_num_jon)
      npositron=prefactor * dum(15)

c DEBUG      
c      write(6,*) 'etae,nele,npos'
c      write(6,100) etae_num,nelectron,npositron

      ne_jon = dim(nelectron - npositron, 0.d0) !in MeV^3

c      write(6,100) tmev, etae, ne/nelectron, ne*mev3tocc   !check1
 100  format(20(1pe11.3,' '))

c      stop                      !check1
      return   
      end

C======================================================================
      REAL*8 FUNCTION eval_num_jon(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist_num_jon/etae_num,beta_num

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae_num              !electron degeneracy parameter
      real*8 beta_num              !relativity parameter
C======================================================================
      real*8 distele_num_jon         !integrand of electrons distribution function
      real*8 distposi_num_jon        !integrand of positrons distribution function


C======================================================================

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
      ENTRY distele_num_jon(x)

      tmp(2)=dexp(dsqrt(x**2+beta_num**(-2.0)) - etae_num) + 1.d0
      tmp(1)= x/(1+beta_num**2*x**2)**(1/4)
      distele_num_jon=tmp(1) / tmp(2)
      

      RETURN

C----------------------------------------------------------------------
      ENTRY distposi_num_jon(x)


      tmp(12)=dexp(dsqrt(x**2+beta_num**(-2.0)) + etae_num) + 1.d0
      tmp(11)= x/(1+beta_num**2*x**2)**(1/4)
      distposi_num_jon=tmp(11) / tmp(12)

      RETURN

      END














C======================================================================
c      Program electron_number_density   !check1
      real*8 function ne_kaz(T11,etae)  !MeV^3

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
      real*8 distintegrate
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
c      real*8 memev              !electron mass in MeV
c      parameter(memev=0.511d0)
      real*8 metilde_num           !me/T

C======================================================================
c     real*8 Qmev                  !qvalue of weak interaction between N,P
c      parameter(Qmev=1.29d0)
      real*8 Qtilde

C======================================================================
c      real*8 mev2K              ! 1MeV = mev2K  K
c      parameter(mev2K=1.16041d10) ! 1MeV = mev2K  K
c      real*8 PI
c      parameter(PI=3.14159265d0)

C=====================================================================
c      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae_num              !electron degeneracy parameter \mue/T

C======================================================================
      external distele_num
      external distposi_num

C======================================================================
c      read(*,*) T11,etae        !check1

      include 'const.dek'

      etae_num=etae

      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde_num=memev/tmev

      if (etae_num.le.metilde_num) then
         dum(31)=1.d0
      else
         dum(31)=Sqrt(etae_num**2-metilde_num**2)   ! dum(31) < 100.d0
      end if

      dum(5) = distintegrate(distele_num)

      nelectron=tmev**3/PI**2 * dum(5)

      dum(15) = distintegrate(distposi_num)

      npositron=tmev**3/PI**2 * dum(15)

c DEBUG      
c      write(6,*) 'etae,nele,npos'
c      write(6,100) etae_num,nelectron,npositron

      ne_kaz = dim(nelectron - npositron, 0.d0) !in MeV^3

c      write(6,100) tmev, etae, ne/nelectron, ne*mev3tocc   !check1
 100  format(20(1pe11.3,' '))

c      stop                      !check1
      return   
      end

C======================================================================
      REAL*8 FUNCTION eval_num(x)

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
c      tmp(2)=ex(Sqrt(x**2+metilde_num**2) - etae_num) + 1.d0
c      tmp(2)=exp(Sqrt(x**2+metilde_num**2) - etae_num) + 1.d0
      tmp(2)=dexp(dsqrt(x**2+metilde_num**2) - etae_num) + 1.d0
      distele_num=tmp(1) / tmp(2)

      RETURN

C----------------------------------------------------------------------
      ENTRY distposi_num(x)

      tmp(11)= x**2
c      tmp(12)=ex(Sqrt(x**2+metilde_num**2) + etae_num) + 1.d0
c      tmp(12)=exp(Sqrt(x**2+metilde_num**2) + etae_num) + 1.d0
      tmp(12)=dexp(dsqrt(x**2+metilde_num**2) + etae_num) + 1.d0
      distposi_num=tmp(11) / tmp(12)

      RETURN

      END


