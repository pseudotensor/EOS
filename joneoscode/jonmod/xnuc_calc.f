C======================================================================
      subroutine xnuc_calc(rho10,T11)
c           xnuc = xnuccalc(rho10,T11)      
C======================================================================
      implicit none
 
C=====================================================================
      include 'kazeos1.dek'
c      common/xnucstuff/xnuc
c      common /variable1/ tdyn,hcm,rhob,tk
C======================================================================
      real*8 xnuc_small,xnuc,T10,E
      real*8 rho10, T11min, T11, disc1, disc2
C======================================================================
c      real*8 tdyn               !dynamical timescale in sec
c      real*8 hcm                !disk half thickness in cm
c      real*8 rhob               !baryon denstiy in g/cc
c      real*8 tk                 !temperature in K
C======================================================================
C==========================  xnuc-calculation part=====================
C======================================================================
      rho10=rhob/1.d10
      T10=T11/1.d1

c      xnuc_small = 30.97d0*rho10**(-0.75d0)*T10**1.125d0*ex(-6.096d0/T10)
c      xnuc_small = 22.16d0*rho10**(-0.75d0)*T10**1.125d0*exp(-8.209d0/T10)
c     Below is what KAZ now things is good/corrected approximation
c      xnuc_small = 11.3226d0*rho10**(-0.75d0)*T10**1.125d0*exp(-8.20899d0/T10)
c      xnuc_small = 26.2d0*rho10**(-0.75d0)*tmev**1.125d0*ex(-7.074d0/tmev)

c     For T11<T11min just set xnuc to 0
      T11min=1.479833198237527*0.1
c     Below is full solution to quartic
c      if(T11.gt.T11min) then
c
         E = exp(1.0)
c
c
         disc1=(9.704093223738058e29*
     -                T11**13.5)/
     -              (E**(9.850785833272978/T11)*
     -                rho10**9) + 
     -             (5.319398204278261e37*T11**18)/
     -              (E**(13.134381111030637/T11)*
     -                rho10**12)
c
         if(disc1.lt.0.0) then
            xnuc = 0.0
            xnuc = min(1.d0, max(xnuc,9.d-10)) !check1
            return
            end if
c
c
c
c
c
          disc2 = (-2.619310999550215e9*T11**4.5)/
     -      (E**(3.283595277757659/T11)*rho10**3*
     -        ((7.293420462497868e18*T11**9)/
     -            (E**(6.567190555515318/T11)*
     -              rho10**6) + 
     -           sqrt(disc1))**
     -         0.3333333333333333) + 
     -     0.26456684199469993*
     -      ((7.293420462497868e18*T11**9)/
     -          (E**(6.567190555515318/T11)*
     -            rho10**6) + 
     -         sqrt(disc1))**0.3333333333333333
c
         if(disc2.lt.0.0) then
            xnuc = 0.0
            xnuc = min(1.d0, max(xnuc,9.d-10)) !check1
            return
            end if
c
c
c
         xnuc = -0.5*sqrt(disc2
     -     ) + 0.5*sqrt((2.619310999550215e9*
     -        T11**4.5)/
     -      (E**(3.283595277757659/T11)*rho10**3*
     -        ((7.293420462497868e18*T11**9)/
     -            (E**(6.567190555515318/T11)*
     -              rho10**6) + 
     -           sqrt(disc1))**
     -         0.3333333333333333) - 
     -     0.26456684199469993*
     -      ((7.293420462497868e18*T11**9)/
     -          (E**(6.567190555515318/T11)*
     -            rho10**6) + 
     -         sqrt(disc1))**0.3333333333333333
     -       + (1.0394742590294718e9*T11**4.5)/
     -      (E**(3.283595277757659/T11)*rho10**3*
     -        sqrt(disc2)))
c         else
c            xnuc = 0.0
c         end if
c
c
c
      xnuc = min(1.d0, max(xnuc,9.d-10))     !check1
c
c
c
      return
      end
