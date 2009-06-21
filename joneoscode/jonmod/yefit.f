

c     pre-supernovae fit
      subroutine computeyefit_tkbased(rhob,tk,yefit)
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 yefit
      real*8 yeiron
C======================================================================
      include 'const.dek'

      yeiron=0.428493

c     From SM, compute Ye for presupernova-like abundances
      if(tk<10**9.532) then
         yefit = 0.5
      else if(tk>8.12d+09) then
         yefit = yeiron
      else
         yefit = (yeiron-0.5)/(8.12d+09-10**9.532)*(tk-10**9.532)+0.5
      end if
      
      return
      end




c     pre-supernovae fit
c     based upon density as should be
      subroutine computeyefit(rhob,tk,yefit)
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 yefit
      real*8 yenucleons,yene,yemg,yesi,yeiron
      real*8 rhobnucleons2ne, rhobne2mg,rhobmg2si,rhobsi2iron
C======================================================================
      include 'const.dek'

c modified FROM SM macro yefitrhob in grbmodel.m
c compute Ye for presupernova-like abundances
      yeiron=0.428493
      yenucleons=0.5
      yene=0.497957
      yemg=0.480252
      yesi=0.462273
      yeiron=0.428493
c
      rhobnucleons2ne=10**6.71
      rhobne2mg=10**7.103
      rhobmg2si=10**8.24
      rhobsi2iron=3.616d+09
c
      if(rhob<rhobnucleons2ne) then
         yefit=yenucleons
      else if(rhob<rhobne2mg) then
         yefit=yene
      else if(rhob<rhobmg2si) then
         yefit=yemg + (rhob-rhobne2mg)*(yesi-yemg)/(rhobmg2si-rhobne2mg)
      else if(rhob<rhobsi2iron) then
         yefit=yesi + (rhob-rhobmg2si)*(yeiron-yesi)/(rhobsi2iron-rhobmg2si)
      else
         yefit=yeiron
      end if
c

      return
      end



c     Compute Y_e [total] assuming alphas or free nucleons only
      subroutine computeyetot(rhob,tk,npratiofree,xnuc,yetot)
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 npratiofree
      real*8 xnuc
      real*8 yetot
C======================================================================
      real*8 npratiotot
c      real*8 malpha
      real*8 mpomalpha
C======================================================================

      include 'const.dek'

c      malpha   = AHe*amu
c     m_p/m_\alpha
      mpomalpha=mp/malpha

c     total neutron proton ratio
c     Only acccounts for He4 in bound state (iron has slightly smaller Ye)
      npratiotot = (npratiofree/(1.0+npratiofree)*xnuc
     1     + 2.0*mpomalpha*(1.0-xnuc))/(1.0/(1.0+npratiofree)*xnuc
     1     + 2.0*mpomalpha*(1.0-xnuc))

      yetot = 1.0/(1.0+npratiotot)


      return
      end




c     Compute Y_e [free,bound] assuming alphas or free nucleons only
c     Computed for input of Y_e[total]
      subroutine computeyefreebound(rhob,tk,tdynorye,xnuc,yefree,yebound)
C======================================================================
      implicit none
      
C======================================================================
c     Passed:
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 tdynorye !Y_e total
      real*8 xnuc
c     Returns:
      real*8 yefree,yebound
C======================================================================
c     Local:
      real*8 npratiobound,npratiofree
      real*8 mpomalpha
C======================================================================

      include 'const.dek'

c      malpha   = AHe*amu
c     m_p/m_\alpha
      mpomalpha=mp/malpha

c     total neutron proton ratio
c     Only acccounts for He4 in bound state (iron has slightly smaller Ye)

c     GODMARK: Notice the limitation of our definition means that
c     bound nuclei always have Y_e=0.5, which is not generally true
c     However, this definition of bound and free is not
c     *really* used when doing whichrnpmethod==1 since assumed using
c     reasonable nuclear EOS that sets Y_e even for bound nuclei
c     GODMARK: Note this is used to set rnp() for qdotNe() that uses xnuc to set rho.  Assume ok since for Xnuc~1 then correct, while for Xnuc~0 qdotNe is not important
      npratiobound = 1.0
      npratiofree  =
     1 xnuc/(xnuc*tdynorye-2.0*mpomalpha*(xnuc-1.0)*(2.0*tdynorye-1.0))-1.0

      if(npratiofree.lt.0.0) then
         npratiofree=0.0
      end if

c      write(*,*) 'inside=',xnuc,tdynorye,mpomalpha


      yefree = 1.0/(1.0+npratiofree)
      yebound = 1.0/(1.0+npratiobound)


      return
      end





c based upon temperature (GODMARK: maybe should be based upon density too)
c     mutot is output and is in mass per particle per amu, not per mb
      subroutine computemutotfit(rhob,xnuc,yefree,tk,mutot)
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
c     Passed:
      real*8 xnuc,yefree
c     Returned:
      real*8 mutot
c     Local:
      real*8 mutot1,mutot2,mutot3,mutot4,mutot5
      real*8 tk12,tk23,tk34,tk45
      real*8 mutotxnuc0
c
C======================================================================
      include 'const.dek'


      call computemutotfit_xnuc0(rhob,tk,mutotxnuc0)

c      write(*,*) 'mutotxnuc0',mutotxnuc0

c     correct for xnuc
      mutot=1.0/((1.0-xnuc)/mutotxnuc0 + xnuc*yefree/(mp/amu) + xnuc*(1.0-yefree)/(mn/amu)   )



      return
      end





c     Stellar Model motivated mutot = \bar{A} over *all* species
c     mutot is output and is in mass per particle per amu, not per mb
      subroutine computemutotfit_xnuc0(rhob,tk,mutotxnuc0)
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 mutotxnuc0
      real*8 mutot1,mutot2,mutot3,mutot4,mutot5
      real*8 tk12,tk23,tk34,tk45
      real*8 rho12,rho23,rho34,rho45
      real*8 tempbased
C======================================================================
      include 'const.dek'


c     Whether to use temp or density based.  Use density based for now since creates monotonic behavior as function of temperature
      tempbased=0
      

c     Compute (as in SM and in kaz-code) the presupernova-like abundances
c     Only used if xnuc~0, as near T=0 and not much of disk or jet
c     Molecular weights only cause minor non-monotonicity that is acceptable
c     Linearly interpolate between to obtain smoother results

      mutot1=AH
      mutot2=AHe
      mutot3=AO16
      mutot4=ASi28
      mutot5=1.0/(0.004779/AH + 0.007451/AHe + 2.44d-07/AMg24 + 5.696d-05/ASi28 + 0.9877/AFe56 )


      if(tempbased.eq.1) then
         tk12=10**7.50
         tk23=10**8.55
         tk34=10**9.45
         tk45=10**9.6

         if(tk>=tk45) then
            mutotxnuc0 =  mutot5
         else if(tk>=tk34) then
            mutotxnuc0 = mutot4 + (tk-tk34)*(mutot5-mutot4)/(tk45-tk34)
         else if(tk>=tk23) then
            mutotxnuc0 = mutot3 + (tk-tk23)*(mutot4-mutot3)/(tk34-tk23)
         else if(tk>=tk12) then
            mutotxnuc0 = mutot2 + (tk-tk12)*(mutot3-mutot2)/(tk23-tk12)
         else
            mutotxnuc0 = mutot1
         end if

      else

	rho12=10**0.0
        rho23=10**3.34
        rho34=10**6.64
        rho45=10**7.09

        if(rhob>=rho45) then
           mutotxnuc0 =  mutot5
        else if(rhob>=rho34) then
           mutotxnuc0 = mutot4 + (rhob-rho34)*(mutot5-mutot4)/(rho45-rho34)
        else if(rhob>=rho23) then
           mutotxnuc0 = mutot3 + (rhob-rho23)*(mutot4-mutot3)/(rho34-rho23)
        else if(rhob>=rho12) then
           mutotxnuc0 = mutot2 + (rhob-rho12)*(mutot3-mutot2)/(rho23-rho12)
        else
           mutotxnuc0 = mutot1
        end if

      end if



      return
      end






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
