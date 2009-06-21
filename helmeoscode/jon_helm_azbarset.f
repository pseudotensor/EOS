c..   ionmax  = number of isotopes in the network
c..   xmass   = mass fractions      (JCM: local rest-mass density fraction of this species)
c..   ymass   = molar fractions     (JCM: set internally by HELM code)
c..   aion    = number of nucleons  (JCM: set to be total rest-mass of species per unit amu -- and this is consistent with how used to define pressure, etc. in this code)
c..   zion    = number of protons   (JCM: = Y_e aion since HELM defines Y_e = zion/aion)



c     See appendix of Timmes & Swesty (2000)
c     Note that $Y_e$ = yein = electrons per baryon = zbar/abarnum
c     abarnum!=abar

c..
c..tests the helmholtz eos routine
c..
c..ionmax  = number of isotopes in the network
c..xmass   = mass fractions
c..ymass   = molar fractions
c..aion    = number of nucleons (mass per unit amu)
c..anumion = number of nucleons
c..zion    = number of protons

      subroutine azbar(xmass,aion,anumion,zion,ionmax,
     1                 xtot,ymass,abar,abarnum,zbar,startbound,abarbound)
      implicit none
      save

c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax

c..output:
c..molar abundances        = ymass(1:ionmax), 
c  JCM (abar must be mass of particle per unit amu given how used in pressure)
c..mean mass number of nucleons = abar
c..mean number of nucleons = abarnum
c..mean nucleon charge     = zbar
c..

c..declare
      integer          i,ionmax,startbound
      double precision xmass(ionmax),aion(ionmax),anumion(ionmax),zion(ionmax)
     1     ,ymass(ionmax),abar,abarnum,zbar,zbarxx,abarbound,ytot1
     2     ,ynum(ionmax),ytot2
     3     ,ytotbound,xtotbound
     4     ,xtot

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
c  JCM:
      ytot2   = 0.0d0
      xtot = 0.0d0
      ytotbound = 0.0d0
      xtotbound = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ynum(i) = xmass(i)/anumion(i)
       ytot1    = ytot1 + ymass(i)
c  JCM:
       ytot2    = ytot2 + ynum(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
       xtot = xtot + xmass(i)

       if(i>=startbound) then
          ytotbound = ytotbound + ymass(i)
          xtotbound = xtotbound + xmass(i) 
       end if

      enddo
c .. JCM: abar is exactly mutot computed in ~/sm/grbmodel.m
      abar   = xtot/ytot1
c    JCM: Now true Ye = zbar/abarnum
c    JCM: GODMARK: Not actually sure how HELM uses Ye.  Might be using a fake Ye that is based upon mass.  Real Ye should only be used when determining fraction of electrons or protons per unit baryon NUMBER.  JCM's Ye should be used when wanting mass fraction of protons per unit baryon mass
      abarnum   = xtot/ytot2
      zbar   = zbarxx * abar

      abarbound = xtotbound/ytotbound

      return
      end



c      subroutine azbarset1(rhob,tk,yein,xmass,aion,anumion,zion,ionmax)
      subroutine azbarset1(rhob,tk,yein,localabar,localabarnum,localzbar,localabarbound)
c     xnuc = xnuccalc(rho10,T11)      
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 localabar,localabarnum,localzbar,localabarbound
      real*8 yein,xnuc,xnuccalc
C======================================================================
      integer ionmax
      parameter        (ionmax=2)
      integer startbound
      parameter        (startbound=2)
      integer i
      real*8 mutotxnuc0,yefit1,yefit2
      double precision computeye
C======================================================================
      double precision xmass(ionmax),aion(ionmax),anumion(ionmax),zion(ionmax),ymass(ionmax),xtot
C======================================================================

      include 'const.dek'

c     Atomic mass (mass of neutrons, protons, electrons in ion)
c     aion = Atomic mass number (number of nucleons (protons and neutrons) in an ion nucleus
c     zion = Atomic number Z (number of protons)

c     Compute Xnuc
      xnuc=xnuccalc(rhob,tk)

c     compute Liebendorfer fit (to be used when xnuc=1)
      yefit2=computeye(rhob)

c     compute presupernovae fit (to be used when xnuc=0)
      call computeyefit(rhob,tk,yefit1)
      if(yein<0) then
c     force to use fit rather than any read-in data
         yein=yefit1
      end if

c     Compute \mu_e (mean molecular weight -- to be used when xnuc=0)
      call computemutotfit_xnuc0(rhob,tk,mutotxnuc0)




c..   set the mass fractions, z's and a's of the composition
      xmass(1) = xnuc
c     the average particle is defined to be half neutron half proton
c      aion(1)  = mb/amu
c     Now I redefined avo throughout code so aion is in terms of mb, not amu
      aion(1)  = 1.0
      anumion(1) = 1.0
c     zion(1)  = 1.0d0
c     1/2 charge per particle (i.e. ye=0.5)
c      zion(1)  = (yein*anumion(1))
      zion(1)  = (yefit2*anumion(1))



c     DEBUG
c     Apparently above fit causes non-monotonic behavior (i.e. as T increases, U,P,S sometimes decrease)
c     This means solution is multi-valued and can't invert properly
c     yein=0.5

c     write(*,*) yein,mutot,xnuc

c..   \alpha particle THROUGH Iron
      xmass(2) = 1.0-xnuc
      aion(2)  = mutotxnuc0
c     for now don't distinguish GODMARK
      anumion(2) = mutotxnuc0
      zion(2)  = (yein*anumion(2))

 

c..   get localabar, localzbar and a few other composition variables
      call azbar(xmass,aion,anumion,zion,ionmax,
     1     xtot,ymass,localabar,localabarnum,localzbar,startbound,localabarbound)





      return
      end









      subroutine azbarset2(rhob,tk,yein,passedionmax,xmass,localabar,localabarnum,localzbar,localabarbound)
c     xnuc = xnuccalc(rho10,T11)      
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 localabar,localabarnum,localzbar,localabarbound
      real*8 yein,xnuc,xnuccalc
      integer passedionmax
      double precision xmass(passedionmax)
C======================================================================
      include 'eosmemory.f'
      integer i
      integer localionmax
      parameter (localionmax=irowmax)
      integer startbound
      parameter        (startbound=2)
      real*8 yefit1,yefit2
      integer          ioni
c     function:
      double precision computeye
C======================================================================
      double precision aion(localionmax),anumion(localionmax),zion(localionmax),ymass(localionmax),xtot
C======================================================================

      include 'const.dek'

c     DEBUG:
c      write(*,*) 'passedionmax',passedionmax,irowmax


c     Compute Xnuc
      xnuc=xnuccalc(rhob,tk)

      if(xnuc.ge.1.0E-2) then
         write(6,04) 'xnuc =', xnuc
      end if
      
      if(xnuc.le.-1E-30) then
         write(6,04) 'xnuc =', xnuc
      end if

04    format(1x,a6,1pe26.15,a6,I3)


c     compute Liebendorfer fit (to be used when xnuc=1)
      yefit2=computeye(rhob)

c     compute presupernovae fit (to be used when xnuc=0)
      call computeyefit(rhob,tk,yefit1)
      if(yein<0) then
         yein=yefit1
      end if




c     Atomic mass (mass of neutrons, protons, electrons in ion)
c     aion = Atomic mass number (number of nucleons (protons and neutrons) in an ion nucleus
c     zion = Atomic number Z (number of protons)
      
c..   set the mass fractions, z's and a's of the composition
c     
c     Note that nucleon pressure is defined as
c     pion = N_A k_b T[K] rho[g/cc] (1/abar)
c     where abar = 1/ytot1
c     where ytot1 = Sum(i) xmass(i)/aion(i)
c     
c     So that aion should be in terms of amu's given N_A in front of pressure
c     
c     Including binding energy in the atomic mass
c     Set 0-point energy as free nucleon so that free nucleons have no additional energy
c     
c     
c     Free nucleons with equal numbers of neutrons and protons (npratio=1) GODMARK
c     This is a good approximation at t=0, but not for late time


c     Choose minimum for Ye between normal state of matter and deleptonized matter
c     Drop charge (zion) since mass roughly same but zion changes by alot per particle

      ioni=1
      xmass(ioni) = xnuc
c     average mass of free nucleon (as defined by Kaz) relative to the atomic mass unit
c      aion(ioni)  = (mn+mp)*0.5/amu
c     Changed definition of avo when used to define this, so now in terms of mb directly
      aion(ioni)  = 1.0
c     particle is half neutron half proton
      anumion(ioni)  = 1.0d0
c     zion(ioni)  = 1.0d0
c     1/2 charge per particle is default unless yefit2 alters this
      zion(ioni)  = 0.5d0
      if(zion(ioni)>yefit2*anumion(ioni)) then
         zion(ioni)  = (yefit2*anumion(ioni))
      end if
c     yein has no knowledge of rest-mass, only nucleon numbers
c     zion(ioni)  = (yein*aion(ioni))


c     Correct mass fractions

c     Hydrogen 1
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_H
      aion(ioni)  = AH
      anumion(ioni)  = 1.0d0
      zion(ioni)  = 1.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))
      
c..   Helium 2
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_He
      aion(ioni)  = AHe
      anumion(ioni)  = 4.0d0
      zion(ioni)  = 2.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))
      
c..   Carbon 12
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_C
      aion(ioni)  = AC12
      anumion(ioni)  = 12.0d0
      zion(ioni)  = 6.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))
      
c..   Oxygen 16
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_O
      aion(ioni)  = AO16
      anumion(ioni)  = 16.0d0
      zion(ioni)  = 8.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))

c..   Neon 20
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_N
      aion(ioni)  = ANe20
      anumion(ioni)  = 20.0d0
      zion(ioni)  = 10.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))
      
c..   Magnesium 22
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_Mg
      aion(ioni)  = AMg24
      anumion(ioni)  = 24.0d0
      zion(ioni)  = 12.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))
      
c..   Silicon 28
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_Si
      aion(ioni)  = ASi28
      anumion(ioni)  = 28.0d0
      zion(ioni)  = 14.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))
      
c..   Iron 56
      ioni=ioni+1
      xmass(ioni) = xmass(ioni)*(1.0-xnuc)
c     xmass(ioni) = X_Fe
      aion(ioni)  = AFe56
      anumion(ioni)  = 56.0d0
      zion(ioni)  = 26.0d0
      if(zion(ioni)>yein*anumion(ioni)) then
         zion(ioni)  = (yein*anumion(ioni))
      end if
c     zion(ioni)  = (yein*aion(ioni))



c..   get localabar, localzbar and a few other composition variables
      call azbar(xmass,aion,anumion,zion,passedionmax,
     1     xtot,ymass,localabar,localabarnum,localzbar,startbound,localabarbound)


      if(dabs(xtot-1.0)>1E-2) then
         write(*,*) 'Something wrong with xtot, which
     1        is presumed in this case to be 1.0 unless
     1        SM is messing up output, xtot=',xtot
         stop
      end if



      return
      end













      subroutine azbarset3(rhob,tk,yein,localabar,localabarnum,localzbar,localabarbound)
c     xnuc = xnuccalc(rho10,T11)      
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 localabar,localabarnum,localzbar,localabarbound
      real*8 yein,xnuc,xnuccalc
C======================================================================
      integer ionmax
      parameter        (ionmax=2)
      integer startbound
      parameter        (startbound=2)
      integer i
      real*8 yefit1,yefit2
      double precision computeye
C======================================================================
      double precision xmass(ionmax),aion(ionmax),anumion(ionmax),zion(ionmax),ymass(ionmax),xtot
C======================================================================

      include 'const.dek'

c     Atomic mass (mass of neutrons, protons, electrons in ion)
c     aion = Atomic mass number (number of nucleons (protons and neutrons) in an ion nucleus
c     zion = Atomic number Z (number of protons)

c     Compute Xnuc
      xnuc=xnuccalc(rhob,tk)


c     compute Liebendorfer fit (to be used when xnuc=1)
      yefit2=computeye(rhob)
      yein=yefit2
      
c..   set the mass fractions, z's and a's of the composition
      xmass(1) = xnuc
c     the particle is half neutron half proton
c      aion(1)  = (mn+mp)*0.5/amu
c     Redefined avo in terms of mb everywhere that needed this definition
      aion(1) = 1.0
      anumion(1) = aion(1)
c     zion(1)  = 1.0d0
c     1/2 charge per particle (i.e. ye=0.5)
c      zion(1)  = (yein*anumion(1))
      zion(1)  = (yefit2*anumion(1))



c     DEBUG
c     Apparently above fit causes non-monotonic behavior (i.e. as T increases, U,P,S sometimes decrease)
c     This means solution is multi-valued and can't invert properly
c     yein=0.5

c     write(*,*) yein,mutot,xnuc

c..   \alpha particle only
      xmass(2) = 1.0-xnuc
      aion(2)  = AHe
      anumion(2) = 4.0
      zion(2)  = (yein*anumion(2))



c..   get localabar, localzbar and a few other composition variables
      call azbar(xmass,aion,anumion,zion,ionmax,
     1     xtot,ymass,localabar,localabarnum,localzbar,startbound,localabarbound)



      return
      end







      subroutine azbarset4(rhob,tk,yein,localabar,localabarnum,localzbar,localabarbound)
c     xnuc = xnuccalc(rho10,T11)      
C======================================================================
      implicit none
      
C======================================================================
c     real*8 tdyn               !dynamical timescale in sec
c     real*8 hcm                !disk half thickness in cm
      real*8 rhob               !baryon denstiy in g/cc
      real*8 tk                 !temperature in K
      real*8 localabar,localabarnum,localzbar,localabarbound
      real*8 yein,xnuc,xnuccalc
C======================================================================
      integer ionmax
      parameter        (ionmax=1)
      integer startbound
      parameter        (startbound=2) ! none!
      integer i
      real*8 yefit1,yefit2
      double precision computeye
C======================================================================
      double precision xmass(ionmax),aion(ionmax),anumion(ionmax),zion(ionmax),ymass(ionmax),xtot
C======================================================================

      include 'const.dek'

c     Atomic mass (mass of neutrons, protons, electrons in ion)
c     aion = Atomic mass number (number of nucleons (protons and neutrons) in an ion nucleus
c     zion = Atomic number Z (number of protons)


c     Compute Xnuc
      xnuc=xnuccalc(rhob,tk)


c     compute Liebendorfer fit (to be used when xnuc=1)
      yefit2=computeye(rhob)
      yein=yefit2

c DEBUG:
c     Pick Y_e
c      yefit2=1.0

      
c..   set the mass fractions, z's and a's of the composition
      xmass(1) = 1.0
c     the particle is half neutron half proton
c     Below is consistent with what LSEOS gives for a certain T,rhob defined by 
      aion(1)  = 2.3929E+02
      anumion(1) = aion(1)
c     zion(1)  = 1.0d0
c     1/2 charge per particle (i.e. ye=0.5)
c      zion(1)  = (yein*anumion(1))
      zion(1)  = (yefit2*anumion(1))


c..   get localabar, localzbar and a few other composition variables
      call azbar(xmass,aion,anumion,zion,ionmax,
     1     xtot,ymass,localabar,localabarnum,localzbar,startbound,localabarbound)




      return
      end




















