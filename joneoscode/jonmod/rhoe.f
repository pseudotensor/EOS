C======================================================================
      real*8 function rhoe(T11,etae) !in MeV^4

C======================================================================
C     This function computes electron energy density, as functions of
C     temperature and the degeneracy parameter of electrons
C     /04/10/3/   K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist_rhoe/etae_rhoe,metilde_rhoe

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
c      real*8 rhoe                 !check1

C======================================================================
      real*8 rho_ele          !number deinsty of e- in MeV^3
      real*8 rho_posi          !number density of e+ in MeV^3

C======================================================================
c      real*8 memev              !electron mass in MeV
c      parameter(memev=0.511d0)
      real*8 metilde_rhoe           !me/T

C======================================================================
c      real*8 Qmev                  !qvalue of weak interaction between N,P
c      parameter(Qmev=1.29d0)
      real*8 Qtilde

C======================================================================
c      real*8 mev2K              ! 1MeV = mev2K  K
c      parameter(mev2K=1.1605d10) ! 1MeV = mev2K  K
c      real*8 PI
c      parameter(PI=3.14159265d0)

C=====================================================================
c      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
c      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae_rhoe              !electron degeneracy parameter \mue/T

C======================================================================
      external distele_rhoe
      external distposi_rhoe

C======================================================================
c      read(*,*) T11,etae        !check1

      include 'const.dek'



      etae_rhoe=etae

      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde_rhoe=memev/tmev

      if (etae_rhoe.le.metilde_rhoe) then
         dum(31)=1.d0
      else
         dum(31)=Sqrt(etae_rhoe**2-metilde_rhoe**2)   ! dum(31) < 100.d0
      end if

      dum(5) = distintegrate(distele_rhoe)

      rho_ele=tmev**4/PI**2 * dum(5)

      dum(15) = distintegrate(distposi_rhoe)

      rho_posi=tmev**4/PI**2 * dum(15)

      rhoe = dim(rho_ele + rho_posi, 0.d0) !in MeV^3

 100  format(20(1pe11.3,' '))

c      stop                      !check1
      return   
      end

C======================================================================
      REAL*8 FUNCTION eval_rhoe(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist_rhoe/etae_rhoe,metilde_rhoe

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae_rhoe              !electron degeneracy parameter

C======================================================================
      real*8 distele_rhoe         !integrand of electrons distribution function
      real*8 distposi_rhoe           !integrand of positrons distribution function


C======================================================================
      real*8 memev              !electron mass in MeV
      parameter(memev=0.511d0)
      real*8 metilde_rhoe           !me/T

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
      ENTRY distele_rhoe(x)
      
      tmp(1)= x**2*Sqrt(x**2+metilde_rhoe**2)
      tmp(2)=exp(Sqrt(x**2+metilde_rhoe**2) - etae_rhoe) + 1.d0
      distele_rhoe=tmp(1) / tmp(2)

      RETURN

C----------------------------------------------------------------------
      ENTRY distposi_rhoe(x)

      tmp(11)= x**2*Sqrt(x**2+metilde_rhoe**2)
      tmp(12)=exp(Sqrt(x**2+metilde_rhoe**2) + etae_rhoe) + 1.d0
      distposi_rhoe=tmp(11) / tmp(12)

      RETURN

      END


