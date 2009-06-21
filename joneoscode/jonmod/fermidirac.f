C======================================================================
c     Compute Fermi-Dirac integral
c
c     Computes F_m(\eta) = \int_0^\infty dx x^m/(exp(x-\eta)+1)
c
      real*8 function fermidiraca(jj,eta)

C======================================================================
C======================================================================
      implicit none

C=====================================================================
      real*8 jj
      real*8 eta               !electron degeneracy parameter \mue/T
C======================================================================
      real*8 ex                 !function
      real*8 zeta3,zeta5,alpha,log2
C======================================================================
      include 'const.dek'

      zeta3=1.20205690315959
c     alpha = N[2*Pi^2/(9*Zeta[3]) - 1, 15]
      alpha=0.824577036826941
      log2=0.693147180559945
      zeta5=1.03692775514337
      
      if(dabs(jj-2.0)<1E-5) then
         if(eta>0.0) then
            fermidiraca =  (1.0/3.0)*(eta**3 + pi**2*eta) + (3.0/2.0)*zeta3*ex(-alpha*eta)
         else
            fermidiraca = (3.0/2.0)*zeta3*ex(eta)
         end if
      else if(dabs(jj-3.0)<1E-5) then
         if(eta>0.0) then
            fermidiraca = (1.0/4.0)*(eta**4 + 2.0*pi**2*eta**2 + 7.0*pi**4/15.0) - 7.0*pi**4/120.0*ex(-eta)
         else
            fermidiraca = 7.0*pi**4/120.0*ex(eta)
         end if
      else if(dabs(jj-4.0)<1E-5) then
         if(eta>0.0) then
            fermidiraca = (-248832.0 + ex(eta)*(1008207.0 + 625.0*ex(eta)*
     1           (-6335.0 + 32.0*ex(eta)*(1375.0 +
     1           27.0*ex(eta)*(-1485.0 + 2.0*ex(eta)*
     1           (720.0 + (-3.0 + 6.0*ex(eta))*eta**5 + 10.0*(-1.0 + 2.0*ex(eta))*eta**3*pi**2 +
     1           7.0*(-1.0 + 2.0*ex(eta))*eta*pi**4 + 30.0*eta**4*log2 + 270.0*eta**2*zeta3 +
     1           675.0*zeta5))))))/(3.24d7*ex(6.0*eta))
         else
            fermidiraca = -24.0*(-ex(eta) + ex(2.0*eta)/32.0 - ex(3.0*eta)/243.0 + ex(4.0*eta)/1024.0 - ex(5.0*eta)/3125.0 +
     1           ex(6.0*eta)/7776.0)
         end if
      else if(dabs(jj-5.0)<1E-5) then
         if(eta>0.0) then
            fermidiraca = (1.0/6.0)*(eta**6 + 5.0*pi**2*eta**4 + 7.0*pi**4*eta**2 + 31.0*pi**6/21.0) - 
     1           31.0*pi**6/252.0*ex(-eta)
         else
            fermidiraca = 31.0*pi**6/252.0*ex(eta)
         end if
      end if
      

      return
      end




C======================================================================
c     Compute very accurate Fermi-Dirac integral
c
c     Computes F_m(\eta) = \int_0^\infty dx x^m/(exp(x-\eta)+1)
c
      real*8 function fermidirac(jj,eta,Q)

C======================================================================
C======================================================================
      implicit none

C=====================================================================
      common /dist_fd/eta_fd,j_fd,Q_fd
      real*8 eta_fd,j_fd,Q_fd

C=====================================================================
C----------PARAMTER.
      integer   iter
      PARAMETER (iter=50)          !Number of gaussian quads.
      real*8 xintd              !function
      real*8 distintegrate      !function
C======================================================================
c      real*8 T11                !temperature in 10^11 K
c      real*8 mue11              !electron chemical potential in 10^11 K
c      real*8 muemev             !electron chemical potential in MeV

C======================================================================
c      real*8 tmev               !temperature in MeV
      real*8 jj
      real*8 eta               !electron degeneracy parameter \mue/T
      real*8 Q
C======================================================================
      external intfd
      real*8 fermidiraca ! function
      real*8 ex                 !function
      integer useapprox
C======================================================================


 

c     Whether to use approximate functions
      useapprox=0

      if(useapprox.eq.1 .AND. dabs(Q)<1E-5) then
         fermidirac = fermidiraca(jj,eta)

      else

         eta_fd = eta
         j_fd = jj
         Q_fd = Q
         fermidirac = distintegrate(intfd)

      end if

c      write(6,100) tmev, etae, gammape  !check2
 100  format(20(1pe11.3,' '))

      return
      end

C======================================================================
      real*8 function intfd(x)

C======================================================================
C======================================================================
      implicit none

C======================================================================
      common /dist_fd/eta_fd,j_fd,Q_fd
      real*8 eta_fd,j_fd,Q_fd

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
      tmp(1)= x**(j_fd)
      tmp(2)=dexp(x + Q_fd - eta_fd) + 1.d0
      intfd=tmp(1) / tmp(2)

      RETURN

      END


