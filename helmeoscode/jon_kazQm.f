      subroutine kazQm_calc(rhobinput,tkinput,xnucinput,etaeinput,Qtot,QNetot)
      implicit none
      save
      include 'kazeos.dek'
      include 'kaz_fakestate.dek'
      include 'const.dek'

c     Returned vars
      real*8 Qtot,QNetot

      real*8 rhobinput,tkinput,xnucinput,etaeinput

      real*8 T
      real*8 ndegrel,ndegnonreln,ndegnonrele
      real*8 ndege,ndegn
      real*8 neisdeg
      real*8 QNenondegval,QNenondeg,QNedeg
      real*8 Qepnondeg,Qepdeg,Qep
      real*8 Qbremdeg,Qbremnondeg,Qbrem
      real*8 gamp,Qplasmon





      T=tkinput
      rhobfree=xnucinput*rhobinput

ccccccccccccccccccccccccccccccccc
c     From kaz.m
c     

      ndegrel=2.0*(kerg*T/(hbarcgs*clight))**3.0/pi**3.0
      ndegnonreln=2.0*(mn*kerg*T/(2.0*pi*hbarcgs**2.0))**(3.0/2.0)
      ndegnonrele=2.0*(me*kerg*T/(2.0*pi*hbarcgs**2.0))**(3.0/2.0)
c
c degeneracy parameter
      if(kerg*T>me*clight**2) then
         ndege=ndegrel
      else
         ndege=ndegnonrele
      end if
      if(kerg*T>mn*clight**2) then
         ndegn=ndegrel
      else
         ndegn=ndegnonreln
      end if
c
      if(nefree>ndege) then
         neisdeg=1
      else
         neisdeg=0
      endif
c
      QNenondegval=9.2E33*(T/1.0E11)**6.0*(rhobfree/1E10)
      QNenondeg=QNenondegval*(1.0-neisdeg)
      QNedeg = 1.1E31*etaeinput**9.0*(T/1.0E11)**9.0*neisdeg
c QNedeg = QNenondeg
c QNe = QNedeg*neisdeg + QNenondeg*(1-neisdeg)
      QNetot = QNedeg + QNenondeg
      Qepnondeg=4.8E33*(T/1E11)**9.0
      Qepdeg = 0.0*T
      if(nefree>ndege) then
         Qep = Qepdeg
      else
         Qep = Qepnondeg
      end if
      Qbremdeg=3.4E33*(T/1E11)**8.0*(rhobfree/1E13)**(1.0/3.0)
      Qbremnondeg=1.5E33*(T/1E11)*5.5*(rhobfree/1E13)**2.0
c
      if(nnfree>ndegnonreln) then
         Qbrem=Qbremdeg
      else
         Qbrem=Qbremnondeg
      end if
      gamp=5.565E-2*sqrt( (pi**2.0+3.0*etaeinput**2.0)/3.0 )
      Qplasmon=1.5E32*(T/1E11)**9.0*gamp**6.0*exp(-gamp)*(1.0+gamp)*(2.0+gamp**2.0/(1.0+gamp))
c
c Qsimplest=QNenondeg+QNedeg+Qepnondeg+Qepdeg+Qbremdeg+Qbremnondeg+Qplasmon
      Qtot=QNetot+Qep+Qbrem+Qplasmon

      return
      end

c      add-in pair capture (can approximately treat as just a sum)
c       kazQmnondegen = 9.2E33*(tkinput/1E11)**6 * (rhobinput/1E10)
c      below very accurate when rhobinput=1E13 and T=10MeV/kerg
c      GODMARK: maybe look at ~/sm/kaz.m for better way to do this
c       kazQmdegen = 1.1E31*(3.5)**9.0 *(tkinput/1E11)**9













      real*8 function Ypcalc(rhobinput,tkinput)
      implicit none
      save

      real*8 rhobinput,tkinput

      real*8 lgrhob,lgT
      integer myuse1,myuse3,myuse4
      real*8 Ypfalse,Yp


      lgrhob=log10(rhobinput)
      lgT=log10(tkinput)
c
c lower left corner
      if((lgrhob<=7.0).AND.(lgT<9.45)) then
         myuse1=1
      else
         myuse1=0
      end if
c lower middle band
cmyuse2=((lgrhob<7.0).AND.(lgT<11.2).AND.(lgT>9.45)) ? 1 : 0
c right part
      if((lgT>9.45).AND.(lgrhob<2.88*(lgT-9.74)+7.3)) then
         myuse3=1
      else
         myuse3=0
      end if
c left upper part
      if((lgrhob>7.0).AND.(lgrhob>2.88*(lgT-9.74)+7.3)) then
         myuse4=1
      else
         myuse4=0
      end if
c
clgT=log10(T)
clgr=log10(rhobinput)
cslope=(11.9782-8.30982)/(10.2258-9.31292)
cfun=8.30982+slope*(lgT-9.31292)
cmyuse4=((rhobinput>10**8.3).AND.(lgr>fun)) ? 1 : 0
c
      if(myuse1.eq.1) then
         Ypfalse=1.0
      else
         if(myuse3.eq.1) then
            Ypfalse=0.5
         else
            if(myuse4.eq.1) then
               Ypfalse=1E-3
            else
               Ypfalse=1E30
            end if
         end if
      end if
c     
      Ypcalc=Ypfalse
c
      
c
      return
      end














      subroutine state_calc(rhobinput,tkinput,xnucinput)
      implicit none
      save
      include 'kazeos.dek'
      include 'kaz_fakestate.dek'
      real*8 rhobinput,tkinput,Ypcalc,xnucinput

c     executable code
      include 'const.dek'

c     get Yp
      Yp=Ypcalc(rhobinput,tkinput)

      rhobfree=xnucinput*rhobinput
      rhoalpha=(1-xnucinput)*rhobinput
c

c     overwrite Kaz answer
c     number of free neutrons to free protons
      npratiofree=(1-Yp)/Yp
c     
c     total number of protons
      rhop=Yp*rhobfree+rhoalpha/2
c     Yp = rhopfree/rhobfree
      rhopfree=Yp*rhobfree
      npfree=rhopfree/mp
      nptotal=1/mp*rhop
      netot=nptotal
c     fully ionized electrons
      nefree=netot
      rhoetot=me*netot
c     
      nnfree=npratiofree*npfree
c     
c     
      
      return
      end
