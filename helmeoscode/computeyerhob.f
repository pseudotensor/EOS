

      double precision function computeye(rhob)
C======================================================================
      implicit none
c     Passed variables
      double precision rhob  ! baryon mass density in g cm^(-3)
c--------------------------------------------------------------------
c electron fraction prescription by Liebendorfer (2005):
      double precision rho1, rho2
      double precision Y1, Y2, Yc
      double precision x_liebendorfer      ! Liebendorfer's x
c--------------------------------------------------------------------
c Calculating the Electron Fraction Ye=Ye(rho) using the analytical
c prescription by M. Liebendorfer (2005) ApJ, 633, 1042. --- G15
c
c----------------------------------------------------------------
c   PARAMETERS FOR LIEBENDORFER's Prescription of Y_e(rho) [in cgs units]
c
      rho1=3.d7
      rho2=2.d13
      Y1=0.5
      Y2=0.278
      Yc=0.035
c
      x_liebendorfer=(2.*dlog(rhob)-dlog(rho2)-dlog(rho1))/
     #        (dlog(rho2)-dlog(rho1))
      x_liebendorfer=max(-1.,min(1.,x_liebendorfer))
c     
      computeye=(Y1+Y2)/2. + x_liebendorfer*(Y2-Y1)/2.+
     #     Yc*(1.-dabs(x_liebendorfer)+4.*dabs(x_liebendorfer)*
     #     (dabs(x_liebendorfer)-0.5)*(dabs(x_liebendorfer)-1.))
c     
      return
      end



c determine whether degenerate (Kohri & Mineshige 2002 equation 10 and 12)
      real*8 function isdegen(mass,rhob,tk)
C======================================================================
      implicit none
c     Passed variables
      real*8 mass,rhob,tk,kb
c--------------------------------------------------------------------
      real*8 relfact,nden,degenrel,degennonrel
c--------------------------------------------------------------------
      include 'const.dek'

      kb=kerg

c     if 0, then non-rel, if 1, then say 
      relfact = kb*tk/(mass*clight**2)
      if(relfact.gt.1.0) then
         relfact=1.0
      end if
      if(relfact.lt.0.0) then
         relfact=0.0
      end if
c     So now relfact goes between 0 and 1 only (should be good enough for this estimation)

c     Only roughly correct for nucleons, otherwise not self-consistent for protons or electrons/positrons
      nden = rhob/mass

      degenrel = nden*pi**2/(1*(kb*tk/(hbarcgs*clight))**3)
      degennonrel = nden/(1*(mass*kb*tk/(2*pi*hbarcgs**2))**(3/2))

c     Estimate general degeneracy parameter
      isdegen=(degenrel*relfact + degennonrel*(1.0-relfact))

      return
      end
