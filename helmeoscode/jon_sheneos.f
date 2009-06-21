




cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     shenbox() is like lsbox() for LSEOS, this returns LSEOS-like quantities for given cgs inputs
c
c     shenbox() inputs cgs values of Y_e, T[K], and \rho_b[g/cc]
c     shenbox() does:
c     1) reads the Shen table
c     2) looks-up where in the table the inputs are
c     3) interpolates the functions within the table
c     4) Stores shen single global values into LS single global values stored in eos_m4c.commononly.inc, el_eos.inc and vector_sneos.dek
c     5) Finally converts LSEOS-type quantities to cgs
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shenbox(yeinput, tkinput, rhobinput, didconverge)
      implicit none
      save
      
c..   does all the stuff needed for calling the lattimer-swesty eos
      
      
c..   bring in the SHEN data structures
      include 'vector_sheneos.dek'
c     contains single globals
      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
      include 'eosparms.f'
c     Constants
      include 'const.dek'

c     JCM: passed parameters
      double precision yeinput,tkinput,rhobinput

c     Local parameters
      double precision rhoi,tki,ypi

      double precision tmin,rhomin,tmax,rhomax, ypmin, ypmax

c     JCM: passed parameter to return to tell if converged or not
      integer didconverge

      integer totalbadinterp
      integer goodlookup

cccccccccccccccccccccccccccccccccccc
c
c     Read Shen table (first since doesn't depend upon anything and other things need result of this)
c
cccccccccccccccccccccccccccccccccccc



      if(whichtable.eq.0) then
c     Read-in Shen table (will only  happen once per entire code run)
         call read_shen_table()
      else
         call read_matlab_shen_table()
      end if



ccccccccccccccccccccccccccccccccccccc
c JCM:
c     Set up input state (only allows cgs input)
c     
c     Note that above is true Ye, while rest-mass version would be using abar_row(index)
c
c     y_inp, temp_cgs, and den_cgs are global and used in any_electron()
c
cccccccccccccccccccccccccccccccccccccc

c     set globals:
      ye_inp   = yeinput
      temp_cgs = tkinput
      den_cgs  = rhobinput
      

ccccccccccccccccccccccccccccccccccccc
c
c     Set limits where want to use table

c     tmin = 2.0*10**(ltkminin)*mev2K
c     Recall that Shen EOS has T=0 data that I interpolated to, so original limit on temperature is excluding T=0.  Hence I expect I can set tmin to whatever lower temperature I interpolated to
      tmin = 1.05*10**(ltkminout)*mev2K
      tmax = 0.95*10**(ltkmaxin)*mev2K
c     tmin = 0
c     Should be chosen to interpolation chooses Shen eos values rather than outside
c     Don't trust Shen inside original density
      rhomin = 1.05*10**(lrhobminin)
c     rhomin = 10**(5.2)
      rhomax = 0.99*10**(lrhobmaxin)

      ypmin=10**lypminout
      ypmax=10**lypmaxout


c     Further restrict and change rho,T limits
      call tableminmaxfixes(tmin,tmax,rhomin,rhomax, ypmin, ypmax)



ccccccccccccccccccccccccccccccccccccc
c JCM:
c     Set up input state (only allows cgs input)
c     
c     Note that above is true Ye, while rest-mass version would be using abar_row(index)
c
cccccccccccccccccccccccccccccccccccccc

c     Assume input in CGS units
      temp_nuc = temp_cgs/mev2K
      temp_nuc_lseos = temp_cgs_lseos * k2mev
      den_nuc  = den_cgs * avo * fm3
      den_nuc_lseos  = den_cgs_lseos * avo * fm3

cccccccccccccccccccccccccccccccccccc
c
c..Call lookup table
c
cccccccccccccccccccccccccccccccccccc



c here "didconverge" is meant to see if within lookup table
      call lookup_shen(den_cgs_lseos, temp_cgs_lseos, ye_inp, goodlookup, rhoi,tki,ypi)


c     DEBUG:
c      write(*,*) den_cgs_lseos,temp_cgs_lseos,ye_inp
c      write(*,*) goodlookup,rhoi,tki,ypi
c      write(*,*) 'tki',tki

      if(goodlookup.eq.1) then

         call interpolate_shen_all(rhoi,tki,ypi,totalbadinterp)
      
      end if



ccccccccccccccccccccccccccccccccccccccccccc
c
c     Define meaning of "convergence"
c
ccccccccccccccccccccccccccccccccccccccccccc
      if((goodlookup.eq.1).AND.(totalbadinterp.eq.0)) then
c         write(*,*) '1didconverge',didconverge
         didconverge=1
      else
c     DEBUG:
c         write(*,*) 'SHENCHECKS',didconverge,goodlookup,totalbadinterp
         didconverge=0
      end if



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Convert from shen???_sing to LSEOS-like values
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*) '2didconverge',didconverge,ye_inp
      call store_row_fromsheneos2lseos(den_cgs,temp_cgs,den_cgs_lseos,temp_cgs_lseos,ye_inp,didconverge)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Convert from LSEOS-like values to cgs values
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call store_row_fromlseos2cgs(didconverge)




      return
      end







      subroutine interpolate_shen_all(rhoi,tki,ypi,totalbadinterp)
      implicit none
      save
      
c..   does all the stuff needed for calling the lattimer-swesty eos
      
      
      include 'eosparms.f'
c..   bring in the SHEN data structures
      include 'vector_sheneos.dek'
c     contains single globals
      include 'vector_eos.single.dek'
c     DEBUG: below added for index
      include 'vector_eos.dek'
      include 'vector_sneos.dek'
c     Constants
      include 'const.dek'

c     Passed parameters
      double precision rhoi,tki,ypi

c     Local parameters
      integer badinterp
c     below are now global, but set here just before used
c      double precision lYpfloor,ltempfloor,Ypfloor,tempfloor,azmin,xmin,nbmin,mstarmin
      integer interptype
      integer whichfunction

c     To be returned parameters
      integer totalbadinterp

      



cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate Shen functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
         totalbadinterp=0
c     Then lookup value within table since within table
c     The 0 = interpolate linear function
c     The 1 = interpolate log10 of function
         interptype=0
         whichfunction=0
         call interpolate_shen(whichfunction, rhoi, tki, ypi,
     1   shenlrhob_row, shenlrhob_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=1
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shennb_row, shennb_sing,interptype,nbmin,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=2
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenlyp_row, shenlyp_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=3
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenyp_row, shenyp_sing,interptype,Ypfloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=4
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenf_row, shenf_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=5
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenebulk_row, shenebulk_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=6
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shensbulk_row, shensbulk_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

c         interptype=1
c         whichfunction=0
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenaheav_row, shenaheav_sing,interptype,azmin,badinterp)
         interptype=0
         whichfunction=7
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenaheav_row, shenaheav_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=8
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenzheav_row, shenzheav_sing,interptype,azmin,badinterp)
c         interptype=0
c         whichfunction=0
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenzheav_row, shenzheav_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=9
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenmstar_row, shenmstar_sing,interptype,mstarmin,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=10
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxneut_row, shenxneut_sing,interptype,xmin,badinterp)
c         interptype=0
c         whichfunction=0
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxneut_row, shenxneut_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=11
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxprot_row, shenxprot_sing,interptype,xmin,badinterp)
c         interptype=0
c         whichfunction=0
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxprot_row, shenxprot_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=12
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxalfa_row, shenxalfa_sing,interptype,xmin,badinterp)
c         interptype=0
c         whichfunction=0
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxalfa_row, shenxalfa_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=13
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxh_row, shenxh_sing,interptype,xmin,badinterp)
c         interptype=0
c         whichfunction=0
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenxh_row, shenxh_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=14
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenpbulk_row, shenpbulk_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=15
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenmunminusmup_row, shenmunminusmup_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=16
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenmup_row, shenmup_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=0
         whichfunction=17
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shenltemp_row, shenltemp_sing,interptype,fakefloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp

         interptype=1
         whichfunction=18
         call interpolate_shen(whichfunction,rhoi, tki, ypi,
     1    shentemp_row, shentemp_sing,interptype,tempfloor,badinterp)
         totalbadinterp=totalbadinterp+badinterp




         if(0) then
            write(*,*) lYpfloor,ltempfloor,Ypfloor,tempfloor,azmin,xmin

c     DEBUG:
            write(*,*) 'totalbadinterp',totalbadinterp
            write(*,*) 'S1',shenlrhob_sing, shennb_sing
            write(*,*) 'S2',shenlyp_sing, shenyp_sing
            write(*,*) 'S3',shenf_sing,shenebulk_sing,shensbulk_sing
            write(*,*) 'S4',shenaheav_sing, shenzheav_sing, shenmstar_sing
            write(*,*) 'S5',shenxneut_sing, shenxprot_sing, shenxalfa_sing, shenxh_sing
            write(*,*) 'S6',shenpbulk_sing,shenmunminusmup_sing,shenmup_sing
            write(*,*) 'S7',shenltemp_sing,shentemp_sing
         end if
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Enforce consistency of interpolation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(totalbadinterp.eq.0) then
            call enforce_consistency_sheneos()
         end if

         if(0) then
c     DEBUG:
            write(*,*) 'totalbadinterp',totalbadinterp
            write(*,*) 'F1',shenlrhob_sing, shennb_sing
            write(*,*) 'F2',shenlyp_sing, shenyp_sing
            write(*,*) 'F3',shenf_sing,shenebulk_sing,shensbulk_sing
            write(*,*) 'F4',shenaheav_sing, shenzheav_sing, shenmstar_sing
            write(*,*) 'F5',shenxneut_sing, shenxprot_sing, shenxalfa_sing, shenxh_sing
            write(*,*) 'F6',shenpbulk_sing,shenmunminusmup_sing,shenmup_sing
            write(*,*) 'F7',shenltemp_sing,shentemp_sing
         end if


c     DEBUG:
c         if(totalbadinterp.ne.0) then
c            write(*,*) 'index=',index
c         end if


         return
         end








ccccccccccccccccccccccccc
c
c Enforce conditions on interpolation since interpolation of a and z and x's separately is unconstrained
c
c Should be consistent with shen_interp.m Matlab code
c
ccccccccccccccccccccccccc

      subroutine enforce_consistency_sheneos()
      implicit none
      save

      include 'eosparms.f'
c     Below includes xnut,xprot, etc. used for ye computation
      include 'eos_m4c.commononly.inc'
c..   bring in the SHEN data structures
      include 'vector_sheneos.dek'
c     contains single globals (also has abar,abarnum,zbar)
c      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
c     contains global parameters and Y_e calculation limits
      include 'vector_eos.dek'
c     Constants
      include 'const.dek'

c     Local parameters
      logical doxcheck,doyecheck,dothermo1check,dothermo2check
c     Local variables
      double precision abarnum,abar,zbar,yelocal
c      double precision yelocal

      double precision ye_inpbackup
      double precision origye,origxneut,origxprot,origxalfa,origxh,origx,origa

      double precision yetrue,temptrue,rhobtrue,nbtrue

c     local variables
      integer xconsistent,yeconsistent


ccccccccccccccccccccccccccccccccccc
c
c Pick which Y_p to use as reference.  Could use interpolated value, but rough interpolation means Y_e rough and need electron EOS to be smooth
c
cccccccccccccccccccccccccccccccccccc
c         yetrue = shenyp_sing
c     ye_inp is LSEOS global in vector_sneos.dek set at shenbox()
      yetrue = ye_inp
      temptrue = temp_nuc_lseos
c      temptrue = shentemp_sing
      rhobtrue = den_cgs_lseos

      if(whichtable.eq.0) then
c     Below is consistent with Shen definition
         nbtrue   = rhobtrue/amu
      else
c     Now assuming Matlab table read-in with already-corrected density
         nbtrue   = rhobtrue/mb
c         nbtrue   = rhobtrue/amu
      end if
c      nbtrue   = shennb_sing


c     Whether to do certain checks
      doxcheck=.TRUE.
      doyecheck=.TRUE.
c      doyecheck=.FALSE.
      dothermo1check=.TRUE.
c      dothermo2check=.TRUE.
      dothermo2check=.FALSE. ! see comments below


 
c     Set globals so can use functions that use them
      ye_inpbackup=ye_inp
      origye=yetrue
      origxneut=shenxneut_sing
      origxprot=shenxprot_sing
      origxalfa=shenxalfa_sing
      origxh=shenxh_sing
      origx=shenzheav_sing/shenaheav_sing
      origa=shenaheav_sing

      ye_inp = origye
      xnut = origxneut
      xprot = origxprot
      xalfa = origxalfa
      xh = origxh
      a     = origa
      x = origx


ccccccccccccccccccccccccc
c     
c     Enforce minimums before corrections (otherwise corrections won't make sense since don't assume <0 for things before-hand)
c     
ccccccccccccccccccccccccc
         if(a<aheavtrust) then
            a=aheavtrust
         end if
         if(x<zheavtrust/aheavtrust) then
            x=zheavtrust/aheavtrust
         end if
         if(xnut<xmin) then
            xnut=xmin
         end if
         if(xprot<xmin) then
            xprot=xmin
         end if
         if(xalfa<xmin) then
            xalfa=xmin
         end if
         if(xh<xhtrust) then
            xh=xhtrust
         end if

c         write(*,*) 'ax1',a,x
c         write(*,*) 'trust',aheavtrust,zheavtrust

      
      if(doxcheck) then
         call enforce_x_consistency(ye_inp, xnut, xprot, xalfa, xh, a, x, xconsistent)
      end if

c      write(*,*) 'ax2',a,x


      if(doyecheck) then
         call compute_nuclear_azbar(xnut, xprot, xalfa, xh, a, x, abarnum,abar,zbar,yelocal)
         
         if(abs(yelocal-yetrue)/(abs(yelocal)+abs(yetrue))>yetolerance) then
c            write(*,*) 'Y1',
c     1           'Problem with yecheck, yelocal=',yelocal,'yp=',yetrue
            call enforce_ye_consistency(ye_inp, xnut, xprot, xalfa, xh, a, x, yeconsistent)
         end if
      end if

c         write(*,*) 'ax2',a,x


ccccccccccccccccccccccccc
c     
c     Enforce minimums after corrections (should only change negligibly small values)
c     
ccccccccccccccccccccccccc
         if(a<aheavtrust) then
            a=aheavtrust
         end if
         if(x<zheavtrust/aheavtrust) then
            x=zheavtrust/aheavtrust
         end if
         if(xnut<xmin) then
            xnut=xmin
         end if
         if(xprot<xmin) then
            xprot=xmin
         end if
         if(xalfa<xmin) then
            xalfa=xmin
         end if
         if(xh<xhtrust) then
            xh=xhtrust
         end if
    
c         write(*,*) 'ax3',a,x


ccccccccccccccccccccccccc
c     
c     Check Y_e value against desired value
c     
ccccccccccccccccccccccccc
         if(doyecheck) then
                  
         call compute_nuclear_azbar(xnut, xprot, xalfa, xh, a, x, abarnum,abar,zbar,yelocal)

c     Can't make right at yetolerance due to machine precision issues, or else would complain again
         if(abs(yelocal-yetrue)/(abs(yelocal)+abs(yetrue))>100.0*yetolerance) then
         write(*,*) 'Y1','STILL problem with yecheck, yelocal=',yelocal,'yp=',yetrue
         write(*,*) 'Y2',xnut,xprot,xalfa,xh,a,x
c         write(*,*) 'Y3',trustheavy,xhdominates
c         write(*,*) 'Y4',changeaheav,changezheav,changexprot,changexalfa,changexh
c         write(*,*) 'Y5',ifxh,ifnonelargest
c         write(*,*) 'Y6',shenxh_sing,shenxneut_sing,shenxprot_sing,shenxalfa_sing
c         write(*,*) 'Y7',shenaheav_sing,shenzheav_sing
c         write(*,*) 'Y8',aheavtrust,xhtrust
c         write(*,*) 'Y9',ifaheavlarger,ifxh,ifxneut,ifxprot,ifxalfa,ifnonelargest
c         write(*,*) 'Y0',origxneut,origxprot,origxalfa,origxh
         write(*,*) 'Yp',yetrue,shenyp_sing
         end if



      end if



cccccccccccccccccccccccccccccccccccccc
c
c     After xcheck and yecheck, restore ye_inp to backup and save other quantities
c
cccccccccccccccccccccccccccccccccccccc
      ye_inp=origye
      shenxneut_sing=xnut
      shenxprot_sing=xprot
      shenxalfa_sing=xalfa
      shenxh_sing=xh
      shenaheav_sing=a
      shenzheav_sing=a*x
c     Also restore backup of ye_inp in case was used outside function (which it presently is)
      ye_inp=ye_inpbackup





      if(dothermo1check) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c 
c
c Thermodynamic consistency relation #3 in Shen EOS guide.ps
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c assume error is in F
c below is $m_{u}$ in MeV from Shen's guide.tex
c     mu = 931.49432;
         shenf_sing = shenebulk_sing - temptrue*shensbulk_sing + mumev - shenmstar_sing



      end if
  

      if(dothermo2check) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 
c
c Thermodynamic consistency relation #1 in Shen EOS guide.ps
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c assume error is in \mu_n and \mu_p equally
c for now assume only in \mu_n for simplicity
c     Assume all error is in shenf_sing, since problems in Shen table at high T and low \rho
         shenmunminusmup_sing = (shenf_sing + shenpbulk_sing/nbtrue - shenmup_sing*yetrue)/(1.0-yetrue) - shenmup_sing


      end if






      return
      end













c Since reading in Shen EOS directly, rhob is defined as rhob = n_b amu instead of rest of code's way, so must correct
c     NOT USING ANYMORE
c     Decided to use Matlab for primary table generation/interpolation
       subroutine read_shen_table()
       implicit none
       save

c..bring in the data structure for HELM or TIMMES EOS
      include 'vector_sheneos.dek'
      include 'const.dek'

c     Input loop variables
      integer ii,jj,kk

c     Indicator of whether already read table (need "save" keyword for this function)
      integer didreadtable
      data didreadtable/0/


      if(didreadtable.eq.1) then
         return
      else
         didreadtable=1
c     And continue
      end if


c Shen EOS assumes 
c
c TM1 parameter set: e.g. Symmetry energy of 36.9MeV and 281MeV incompressibility
c
c
c Objects:
c
c  1) log10(rhob) [g/cm^3] : log10 Baryon rest-mass density
c  2) n_b [fm^{-3}]        : Baryon number density (rhob = n_b*amu)
c  3) log10(Yp)            : log10 Proton number
c  4) Y_p                  : Proton fraction
c                            Y_p = N_p/N_B in general, or Y_p = (n_p+2n_\alpha)/(n_n+n_p+4n_\alpha)
c  5) F [MeV]              : Free energy per baryon
c                            F = f/n_b - M  where f=total free energy and M=938MeV
c  6) E [MeV]              : Internal energy per baryon
c                            E = \ep/n_b - m_u where \ep is total energy
c                            m_u = 931.49432MeV
c  7) S [k_b]              : Entropy per baryon
c                            S = s/n_b where s=entropy density
c  8) aheav                : Mass number of heavy nucleus
c  9) zheav                : Charge number of heavy nucleus
c 10) M* [MeV]             : Effective mass
c 11) xneut                : Free neutron mass fraction
c 12) xprot                : Free proton mass fraction
c 13) xalfa                : Alpha mass fraction
c 14) xheav                : Heavy nucleus mass fraction
c 15) P [MeV/fm^3]         : Pressure
c 16) \tilde{\mu_n} [MeV]  : Chemical potential of neutrons relative to free nucleon mass M
c                            n_n = (1-Y_p)n_b
c 17) \tilde{\mu_p} [MeV]  : Chemical potential of protons relative to free nucleon mass M
c                            n_p = Y_p n_b
c
c Consistency requires:
c
c 1) F + pulk/n_b = \tilde{\mu_n} (1-Y_p) + \mu_p Y_p
c
c 2) xneut + xprot + xalfa + xheav = 1
c
c 3) F = E - TS + m_u - M  : original table has these conserved to 1 part in 1/0.001
c
c
c We give the resulting EOS by three tables, which are named

c (1) eos.tab (main EOS table, size: 54.5MB)
c \item temperature $T (MeV)$: 
c $ -1.0 \leq \log_{10}(T) \leq 2.0$
c mesh of $\log_{10}(T) \sim 0.1      $
c \item proton fraction $Y_p$:
c $ -2.00 \leq \log_{10}(Y_p) \leq -0.25 $
c mesh of $\log_{10}(Y_p) \sim 0.025 $
c \item baryon mass density $\rho_B (g/cm^3)$:
c $ 5.1 \leq \log_{10}(\rho_B) \leq 15.4$
c mesh of $\log_{10}(\rho_B) \sim 0.1 $

c (2) eos.t00 (EOS at $T=0$, same range of $Y_p$ and $\rho_B$, size: 1.76MB) 

c (3) eos.yp0 (EOS at $Y_p=0$, same range of $T$ and $\rho_B$, size: 0.79MB) \\
c We putted them in the website after compressed by 'gzip' \\
c http://www-server.rcnp.osaka-u.ac.jp/$^{\sim}$shen/table/  \\

c We write the table in the following order, first fix $T$ which is noted
c at the top of each block, second fix $Y_p$, third fix $\rho_B$.
c The blocks with different $T$ are divided by line 'ccccccc'.
c For each $T$, $Y_p$, and $\rho_B$, we tabulate all the  quantities  
c in one line in the order
c
c
c
c Interpretation:
c Data exist in Lines:
c
c StartTblock=4+i*N to EndTblock=7457+i*N with N=7458 for i=0,30  where line# starts at 1c
c
c Thus there are: 31 temperature points in blocks of size N lines with N-3 data lines
c
c This N corresponds to (rhob)*(nYp+1) where the +1 on nYp is because of a line break in each Yp block
c
c Within each temperature block there are nYp=71 Yp blocks with single line separations, accessed via:
c
c StartYpblock=4 + i*N + j*(nrhob+1)  EndYpblock=7457 + i*N + j*(nrhob+1) + nrhob
c 
c They note that the grid is not exactly uniform in log due to difficulties in computing the EOS for specific parameters
c 
c
c We store T=0 and Yp=0 tables in same array and perform linear interpolation on T or Yp if those values used
c Otherwise log-log-log interpolation is used
c
c So the grid is rectangular so can be read-in directly as a block into memory after first eliminating unecessary lines/comments via:
c
c egrep '^ [0-9]' eos.tab > sheneos.tab
c egrep '^ [0-9]' eos.t00 > sheneos.t00
c egrep '^ [0-9]' eos.yp0 > sheneos.yp0
c
c Recall order is: rhob fastest, Yp next fastest, T slowest
c


c         T11 = tk/1.d11         !in 10^{11}K
c         T10 = T11*10.d0        !in 10^{10}K
c         tmev=T10/(1.1605d0)    !in MeV
c         rho10=rhob/1.d10       !in 10^{10} g/cc




cccccccccccccccccccccccccccccccccccccccc
c
c     READ IN PRIMARY TABLE
c
cccccccccccccccccccccccccccccccccccccccc

      open (unit=3,file='sheneos.tab')

      do kk=2,numt
         do jj=2,numYp
            do ii=1,numrhob

               read (3,*)  shenlrhob_row(kk,jj,ii), shennb_row(kk,jj,ii),
     2              shenlyp_row(kk,jj,ii),   shenyp_row(kk,jj,ii),
     3              shenf_row(kk,jj,ii),     shenebulk_row(kk,jj,ii),
     3              shensbulk_row(kk,jj,ii),
     4              shenaheav_row(kk,jj,ii), shenzheav_row(kk,jj,ii),
     5              shenmstar_row(kk,jj,ii), shenxneut_row(kk,jj,ii),
     5              shenxprot_row(kk,jj,ii),
     6              shenxalfa_row(kk,jj,ii), shenxh_row(kk,jj,ii),
     7              shenpbulk_row(kk,jj,ii), shenmunminusmup_row(kk,jj,ii),
     8              shenmup_row(kk,jj,ii)
	
            end do
         end do
      end do

c     Close the input file
      close(3)




cccccccccccccccccccccccccccccccccccccccc
c
c     READ IN T=0 table
c
cccccccccccccccccccccccccccccccccccccccc

      open (unit=3,file='sheneos.t00')

      do kk=1,1
         do jj=2,numYp
            do ii=1,numrhob

               read (3,*)    shenlrhob_row(kk,jj,ii), shennb_row(kk,jj,ii),
     2              shenlyp_row(kk,jj,ii),   shenyp_row(kk,jj,ii),
     3              shenf_row(kk,jj,ii),     shenebulk_row(kk,jj,ii),
     3              shensbulk_row(kk,jj,ii),
     4              shenaheav_row(kk,jj,ii), shenzheav_row(kk,jj,ii),
     5              shenmstar_row(kk,jj,ii), shenxneut_row(kk,jj,ii),
     5              shenxprot_row(kk,jj,ii),
     6              shenxalfa_row(kk,jj,ii), shenxh_row(kk,jj,ii),
     7              shenpbulk_row(kk,jj,ii), shenmunminusmup_row(kk,jj,ii),
     8              shenmup_row(kk,jj,ii)

c     Set temperature to be 0
               shentemp_row(kk,jj,ii)=0.0
	

            end do
         end do
      end do

c     Close the input file
      close(3)


cccccccccccccccccccccccccccccccccccccccc
c
c     READ IN Yp=0 table
c
cccccccccccccccccccccccccccccccccccccccc

      open (unit=3,file='sheneos.yp0')

      do kk=2,numt
         do jj=1,1
            do ii=1,numrhob

               read (3,*)    shenlrhob_row(kk,jj,ii), shennb_row(kk,jj,ii),
     2              shenlyp_row(kk,jj,ii),   shenyp_row(kk,jj,ii),
     3              shenf_row(kk,jj,ii),     shenebulk_row(kk,jj,ii),
     3              shensbulk_row(kk,jj,ii),
     4              shenaheav_row(kk,jj,ii), shenzheav_row(kk,jj,ii),
     5              shenmstar_row(kk,jj,ii), shenxneut_row(kk,jj,ii),
     5              shenxprot_row(kk,jj,ii),
     6              shenxalfa_row(kk,jj,ii), shenxh_row(kk,jj,ii),
     7              shenpbulk_row(kk,jj,ii), shenmunminusmup_row(kk,jj,ii),
     8              shenmup_row(kk,jj,ii)


c     Set Yp to be 0
c     log10(Yp) shouldn't be used
               shenlyp_row(kk,jj,ii)=-1E50
               shenyp_row(kk,jj,ii)=0.0
	

            end do
         end do
      end do

c     Close the input file
      close(3)


cccccccccccccccccccccccccccccccccccccccc
c
c     Set corner for T=0 AND Yp=0 as average of T=0 and Yp=0 cases
c
cccccccccccccccccccccccccccccccccccccccc

      do kk=1,1
         do jj=1,1
            do ii=1,numrhob
               shenlrhob_row(kk,jj,ii) = 0.5*(shenlrhob_row(1,2,ii)+shenlrhob_row(2,1,ii))
               shennb_row(kk,jj,ii) = 0.5*(shennb_row(1,2,ii)+shennb_row(2,1,ii))
               shenlyp_row(kk,jj,ii) = 0.5*(shenlyp_row(1,2,ii)+shenlyp_row(2,1,ii))
               shenyp_row(kk,jj,ii) = 0.5*(shenyp_row(1,2,ii)+shenyp_row(2,1,ii))
               shenf_row(kk,jj,ii) = 0.5*(shenf_row(1,2,ii)+shenf_row(2,1,ii))
               shenebulk_row(kk,jj,ii) = 0.5*(shenebulk_row(1,2,ii)+shenebulk_row(2,1,ii))
               shensbulk_row(kk,jj,ii) = 0.5*(shensbulk_row(1,2,ii)+shensbulk_row(2,1,ii))
               shenaheav_row(kk,jj,ii) = 0.5*(shenaheav_row(1,2,ii)+shenaheav_row(2,1,ii))
               shenzheav_row(kk,jj,ii) = 0.5*(shenzheav_row(1,2,ii)+shenzheav_row(2,1,ii))
               shenmstar_row(kk,jj,ii) = 0.5*(shenmstar_row(1,2,ii)+shenmstar_row(2,1,ii))
               shenxneut_row(kk,jj,ii) = 0.5*(shenxneut_row(1,2,ii)+shenxneut_row(2,1,ii))
               shenxprot_row(kk,jj,ii) = 0.5*(shenxprot_row(1,2,ii)+shenxprot_row(2,1,ii))
               shenxalfa_row(kk,jj,ii) = 0.5*(shenxalfa_row(1,2,ii)+shenxalfa_row(2,1,ii))
               shenxh_row(kk,jj,ii) = 0.5*(shenxh_row(1,2,ii)+shenxh_row(2,1,ii))
               shenpbulk_row(kk,jj,ii) = 0.5*(shenpbulk_row(1,2,ii)+shenpbulk_row(2,1,ii))
               shenmunminusmup_row(kk,jj,ii) = 0.5*(shenmunminusmup_row(1,2,ii)+shenmunminusmup_row(2,1,ii))
               shenmup_row(kk,jj,ii) = 0.5*(shenmup_row(1,2,ii)+shenmup_row(2,1,ii))
            end do
         end do
      end do





cccccccccccccccccccccccccccccccccccccccc
c
c     Convert to cgs units
c
cccccccccccccccccccccccccccccccccccccccc

      do kk=1,numt
         do jj=1,numYp
            do ii=1,numrhob

               shenlrhob_row(kk,jj,ii)  = log10(10**(shenlrhob_row(kk,jj,ii))/amu*mb)
               shennb_row(kk,jj,ii)     = shennb_row(kk,jj,ii)     * ( 1.0/fm3  )
c               shenlyp_row(kk,jj,ii)    = shenlyp_row(kk,jj,ii)    * ( 1.0  )
c               shenyp_row(kk,jj,ii)     = shenyp_row(kk,jj,ii)     * ( 1.0  )
               shenf_row(kk,jj,ii)      = shenf_row(kk,jj,ii)      * ( mev2erg  )
               shenebulk_row(kk,jj,ii)  = shenebulk_row(kk,jj,ii)  * ( mev2erg  )
c     Assume their [k_b] means k_b in MeV/MeV/baryon and we want to convert to erg/K/baryon
c     So just multiply by k_b in cgs units
               shensbulk_row(kk,jj,ii)  = shensbulk_row(kk,jj,ii)  * ( kerg  )
               shenaheav_row(kk,jj,ii)  = shenaheav_row(kk,jj,ii)  * ( 1.0  )
               shenzheav_row(kk,jj,ii)  = shenzheav_row(kk,jj,ii)  * ( 1.0  )
               shenmstar_row(kk,jj,ii)  = shenmstar_row(kk,jj,ii)  * ( mev2erg  )
               shenxneut_row(kk,jj,ii)  = shenxneut_row(kk,jj,ii)  * ( 1.0  )
               shenxprot_row(kk,jj,ii)  = shenxprot_row(kk,jj,ii)  * ( 1.0  )
               shenxalfa_row(kk,jj,ii)  = shenxalfa_row(kk,jj,ii)  * ( 1.0  )
               shenxh_row(kk,jj,ii)     = shenxh_row(kk,jj,ii)     * ( 1.0  )
               shenpbulk_row(kk,jj,ii)  = shenpbulk_row(kk,jj,ii)  * ( mev2erg/fm3  )
               shenmunminusmup_row(kk,jj,ii)    = shenmunminusmup_row(kk,jj,ii)    * ( mev2erg  )
               shenmup_row(kk,jj,ii)    = shenmup_row(kk,jj,ii)    * ( mev2erg  )
	
            end do
         end do
      end do


cccccccccccccccccccccccccccccccccccccccc
c
c     Convert to energy densities for F, E, and S
c
c     So then $F$, $E$, $P$ are same units as $\rho_b c^2$ and $s/k_b=S/n_b/k_b$
c
cccccccccccccccccccccccccccccccccccccccc

      do kk=1,numt
         do jj=1,numYp
            do ii=1,numrhob

               shenf_row(kk,jj,ii)      = shenf_row(kk,jj,ii)      * ( shennb_row(kk,jj,ii) )
               shenebulk_row(kk,jj,ii)  = shenebulk_row(kk,jj,ii)  * ( shennb_row(kk,jj,ii) )
c     Entropy is now erg/K/cc instead of erg/K/baryon
               shensbulk_row(kk,jj,ii)  = shensbulk_row(kk,jj,ii)  * ( shennb_row(kk,jj,ii) )
	
            end do
         end do
      end do








      return
      end


















c Although "Shen EOS", rhob is defined in terms of n_b m_b rather than n_b amu if using Matlab table

cccccccccccccccccccccccc
c
c     USING THIS function to read-in table outputted from shen_interp.m matlab function
c
cccccccccccccccccccccccccccc
c     This matlab function co-locates \rho_b and T so perfect interpolation
c     Still linearly interpolate to get everything, but only Yp has linear error
      subroutine read_matlab_shen_table()
      implicit none
      save
      
c     next 2 are for loop parameters
      include 'kazeos.loopvars.dek'
      include 'kazeos.loopparms.dek'
c..   bring in the data structure for HELM or TIMMES EOS
      include 'vector_sheneos.dek'
      include 'const.dek'
       
c     Input loop variables
      integer ii,jj,kk

c     Choose which units to convert read-in table into
c     0 = nuclear (no change)
c     1 = CGS (not setup to use this, so don't use)
      integer whichunits
      data whichunits/0/

      integer didreadtable
      data didreadtable/0/

c     Local loop-check parameters
      integer index1,index2,index3


c     DEBUG:
c      write(*,*) 'didreadtable=',didreadtable




      if(didreadtable.eq.1) then
         return
      else
         write(*,*) "Reading Shen table"
         didreadtable=1
c     And continue
      end if


c     Same objects as normal Shen table but with 2 extra columns corresponding to:
c     log10(T[MeV]) T[MeV]
c
c
c     Recall order is: rhob fastest, Yp next fastest, T slowest
c
c         T11 = tk/1.d11         !in 10^{11}K
c         T10 = T11*10.d0        !in 10^{10}K
c         tmev=T10/(1.1605d0)    !in MeV
c         rho10=rhob/1.d10       !in 10^{10} g/cc




cccccccccccccccccccccccccccccccccccccccc
c
c     READ IN DATA HEADER (20 things)
c
cccccccccccccccccccccccccccccccccccccccc

      open (unit=3,file='sheneos.head')
      read (3,*) ncin,ncout,nrhobin,nypin,ntkin,nrhobout,nypout,ntkout
     1     ,lrhobminin,lrhobmaxin,lypminin,lypmaxin,ltkminin,ltkmaxin
     1     ,lrhobminout,lrhobmaxout,lypminout,lypmaxout
     1     ,ltkminout,ltkmaxout
      close (3)

      
      if(      (ncin.ne.ncinexpected)
     1     .OR.(ncout.ne.ncoutexpected)
     1     .OR.(nrhobin.ne.numrhoborig)
     1     .OR.(nypin.ne.numyporig)
     1     .OR.(ntkin.ne.numtorig)
     1
c     Check if table is consistent with expected matlab size
     1     .OR.(nrhobout.ne.numrhobmatlab)
c     Below not yet
c     1     .OR.(nypout.ne.numypmatlab)
     1     .OR.(ntkout.ne.numtmatlab)
     1
c     Check if table is consistent with expected array size
     1     .OR.(nrhobout.ne.numrhob)
c     Below not yet
c     1     .OR.(nypout.ne.numyp)
     1     .OR.(ntkout.ne.numt)
     1     ) then
      write(*,*) "Shen EOS header has incorrect information"
      stop
      end if


c     Check if table is consistent with Kaz loop ranges
      if(
     1     (nrhobout.ne.nrhob)
c     Below not yet
c     1     .OR.(nypout.ne.nyp)
     1     .OR.(ntkout.ne.ntk)
     1     ) then
      write(*,*) 'nrhobout=',nrhobout,'nrhob=',nrhob
      write(*,*) 'ntkout=',ntkout,'ntk=',ntk
      write(*,*) 'Shen EOS header is not setup as identical
     1     as Kaz size for desired density/temperature points'
      write(*,*) 'This is ok if size is different
     1     but points in Shen table match'
      write(*,*) 'Or this is ok if want to just use
     1     linear interpolation regardless of matching'
c      stop
      end if

c     DEBUG:
c      write(*,*) nrhobout,nrhob,ntkout,ntk

      

cccccccccccccccccccccccccccccccccccccccc
c
c     READ IN PRIMARY MATLAB TABLE (ncout=19 things)
c
cccccccccccccccccccccccccccccccccccccccc

      open (unit=3,file='sheneos.dat')

      do kk=1,numt
         do jj=1,numyp
            do ii=1,numrhob

               read (3,*)  index1,index2,index3,
     1              shenwithintable_row(kk,jj,ii),
     1              shenlrhob_row(kk,jj,ii), shennb_row(kk,jj,ii),
     2              shenlyp_row(kk,jj,ii),   shenyp_row(kk,jj,ii),
     3              shenf_row(kk,jj,ii),     shenebulk_row(kk,jj,ii),
     3              shensbulk_row(kk,jj,ii),
     4              shenaheav_row(kk,jj,ii), shenzheav_row(kk,jj,ii),
     5              shenmstar_row(kk,jj,ii), shenxneut_row(kk,jj,ii),
     5              shenxprot_row(kk,jj,ii),
     6              shenxalfa_row(kk,jj,ii), shenxh_row(kk,jj,ii),
     7              shenpbulk_row(kk,jj,ii), shenmunminusmup_row(kk,jj,ii),
     8              shenmup_row(kk,jj,ii),
     9              shenltemp_row(kk,jj,ii), shentemp_row(kk,jj,ii)

c     Check consistency of table indicies
               if( (index3==kk).AND.(index2==jj).AND.(index1==ii)) then
c     then good
                  else
                     write(*,*) 'Table not consistent'
     1                    ,index1,index2,index3,kk,jj,ii
                     stop
                  end if

            end do
         end do
      end do

c     Close the input file
      close(3)




      if(whichunits.eq.0) then
cccccccccccccccccccccccccccccccccccccccc
c     
c     Keep nuclear units and convert to LSEOS units
c     
cccccccccccccccccccccccccccccccccccccccc

         do kk=1,numt
            do jj=1,numyp
               do ii=1,numrhob

                  shenlrhob_row(kk,jj,ii)  = shenlrhob_row(kk,jj,ii)  * ( 1.0 )
                  shennb_row(kk,jj,ii)     = shennb_row(kk,jj,ii)     * ( 1.0  )
c     shenlyp_row(kk,jj,ii)    = shenlyp_row(kk,jj,ii)    * ( 1.0  )
c     shenyp_row(kk,jj,ii)     = shenyp_row(kk,jj,ii)     * ( 1.0  )
                  shenf_row(kk,jj,ii)      = shenf_row(kk,jj,ii)      * ( 1.0  )
                  shenebulk_row(kk,jj,ii)  = shenebulk_row(kk,jj,ii)  * ( 1.0  )
                  shensbulk_row(kk,jj,ii)  = shensbulk_row(kk,jj,ii)  * ( 1.0  )
c                  shenaheav_row(kk,jj,ii)  = shenaheav_row(kk,jj,ii)  * ( 1.0  )
c                  shenzheav_row(kk,jj,ii)  = shenzheav_row(kk,jj,ii)  * ( 1.0  )
                  shenmstar_row(kk,jj,ii)  = shenmstar_row(kk,jj,ii)  * ( 1.0  )
c                  shenxneut_row(kk,jj,ii)  = shenxneut_row(kk,jj,ii)  * ( 1.0  )
c                  shenxprot_row(kk,jj,ii)  = shenxprot_row(kk,jj,ii)  * ( 1.0  )
c                  shenxalfa_row(kk,jj,ii)  = shenxalfa_row(kk,jj,ii)  * ( 1.0  )
c                  shenxh_row(kk,jj,ii)     = shenxh_row(kk,jj,ii)     * ( 1.0  )
                  shenpbulk_row(kk,jj,ii)  = shenpbulk_row(kk,jj,ii)  * ( 1.0  )
                  shenmunminusmup_row(kk,jj,ii)    = shenmunminusmup_row(kk,jj,ii)    * ( 1.0  )
                  shenmup_row(kk,jj,ii)    = shenmup_row(kk,jj,ii)    * ( 1.0  )
                  shenltemp_row(kk,jj,ii)  = log10(10**shenltemp_row(kk,jj,ii)  * ( 1.0 ))
                  shentemp_row(kk,jj,ii)   = shentemp_row(kk,jj,ii)   * ( 1.0 )
                  
               end do
            end do
         end do

      end if



      if(whichunits.eq.1) then
cccccccccccccccccccccccccccccccccccccccc
c     
c     Convert to cgs units (If done, not correct for rest of code)
c     
cccccccccccccccccccccccccccccccccccccccc

         do kk=1,numt
            do jj=1,numyp
               do ii=1,numrhob

                  shenlrhob_row(kk,jj,ii)  = shenlrhob_row(kk,jj,ii)  * ( 1.0 )
                  shennb_row(kk,jj,ii)     = shennb_row(kk,jj,ii)     * ( 1.0/fm3  )
c     shenlyp_row(kk,jj,ii)    = shenlyp_row(kk,jj,ii)    * ( 1.0  )
c     shenyp_row(kk,jj,ii)     = shenyp_row(kk,jj,ii)     * ( 1.0  )
                  shenf_row(kk,jj,ii)      = shenf_row(kk,jj,ii)      * ( mev2erg  )
                  shenebulk_row(kk,jj,ii)  = shenebulk_row(kk,jj,ii)  * ( mev2erg  )
c     Assume their [k_b] means k_b in MeV/K and we want to convert to erg/K
                  shensbulk_row(kk,jj,ii)  = shensbulk_row(kk,jj,ii)  * ( mev2erg  )
                  shenaheav_row(kk,jj,ii)  = shenaheav_row(kk,jj,ii)  * ( 1.0  )
                  shenzheav_row(kk,jj,ii)  = shenzheav_row(kk,jj,ii)  * ( 1.0  )
                  shenmstar_row(kk,jj,ii)  = shenmstar_row(kk,jj,ii)  * ( mev2erg  )
                  shenxneut_row(kk,jj,ii)  = shenxneut_row(kk,jj,ii)  * ( 1.0  )
                  shenxprot_row(kk,jj,ii)  = shenxprot_row(kk,jj,ii)  * ( 1.0  )
                  shenxalfa_row(kk,jj,ii)  = shenxalfa_row(kk,jj,ii)  * ( 1.0  )
                  shenxh_row(kk,jj,ii)     = shenxh_row(kk,jj,ii)     * ( 1.0  )
                  shenpbulk_row(kk,jj,ii)  = shenpbulk_row(kk,jj,ii)  * ( mev2erg/fm3  )
                  shenmunminusmup_row(kk,jj,ii)    = shenmunminusmup_row(kk,jj,ii)    * ( mev2erg  )
                  shenmup_row(kk,jj,ii)    = shenmup_row(kk,jj,ii)    * ( mev2erg  )
                  shenltemp_row(kk,jj,ii)  = log10(10**shenltemp_row(kk,jj,ii)  * ( mev2K ))
                  shentemp_row(kk,jj,ii)   = shentemp_row(kk,jj,ii)   * ( mev2K )
                  
               end do
            end do
         end do

      end if








cccccccccccccccccccccccccccccccccccccccc
c
c     DO NOT do for now since writing in LSEOS form and these will later be converted to densities
c
c     Convert to energy densities for F, E, and S
c
c     So then $F$, $E$, $P$ are same units as $\rho_b c^2$ and $S/k_b$
c
cccccccccccccccccccccccccccccccccccccccc

c      do kk=1,numt
c         do jj=1,numYp
c            do ii=1,numrhob
c
c               shenf_row(kk,jj,ii)      = shenf_row(kk,jj,ii)      * ( shennb_row(kk,jj,ii) )
c               shenebulk_row(kk,jj,ii)  = shenebulk_row(kk,jj,ii)  * ( shennb_row(kk,jj,ii) )
c               shensbulk_row(kk,jj,ii)  = shensbulk_row(kk,jj,ii)  * ( shennb_row(kk,jj,ii) )
c	
c            end do
c         end do
c      end do








      return
      end












cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Lookup floating point position in Matlab version of Shen table 
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lookup_shen(den_cgs_lseos, temp_cgs_lseos, ye_inp, didconverge, rhobi, tki, ypi)
      implicit none

c     Passed parameters
c     double precision temp_cgs,ye_inp,den_cgs
c     didconverge is SET HERE and tells rest of code if solution found within lookup table
      integer didconverge

c     Passed parameters
      double precision den_cgs_lseos,temp_cgs_lseos,ye_inp

c     Parameters to return
      double precision rhobi,tki,ypi

c     Local parameters
      double precision rhob,yp,tk,lrhob,lyp,ltk
c     Nulcear unit versions
      double precision tin,yein,dnsin

      include 'const.dek'

c     Bring in Shen table parameters
      include 'vector_sheneos.dek'


      

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Convert CGS and nuclear units into table units
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Convert input values to nuclear units for lookup
      tin = temp_cgs_lseos*k2mev
      yein=ye_inp
      dnsin=den_cgs_lseos*(avo*fm3)


c     Note that lookup table was formed in mixed units
c     Density in CGS units in g/cc
c     Temperature in nuclear units in MeV
      rhob=den_cgs_lseos ! g/cc
      yp=ye_inp
      tk=tin   ! MeV

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Given rhob,yp,tk then look up kk,jj,ii to use
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     old simple log10 way:
c     i = (x - x0)*[(N-1)/(x1-x0)] such that if x=x0, then i=0, if x=x1, then i=N-1
c     Fortran form starting with 1 through nrhob,ntk,nyp

      lrhob = log10(rhob)
      rhobi = 1 + (lrhob - lrhobminout)*(nrhobout-1)/(lrhobmaxout-lrhobminout)

      lyp = log10(yp)
      ypi = 1 + (lyp - lypminout)*(nypout-1)/(lypmaxout-lypminout)
 
      ltk = log10(tk)
      tki = 1 + (ltk - ltkminout)*(ntkout-1)/(ltkmaxout-ltkminout)

c     DEBUG:
c      write(*,*) rhob,lrhob,lrhobminout,nrhobout,lrhobmaxout,rhobi
c      write(*,*) yp,lyp,lypminout,nypout,lypmaxout,ypi
c      write(*,*) tk,ltk,ltkminout,ntkout,ltkmaxout,tki
      

ccccccccccccccccc
c
c     Don't allow answer if beyond original table or beyond Matlab version of table
c
ccccccccccccccccc
      if(
     1     ((lrhob>lrhobmaxin).OR.(lrhob<lrhobminin))
     1     .OR.((lrhob>lrhobmaxout).OR.(lrhob<lrhobminout))
     1     .OR.((lyp>lypmaxin).OR.(yp<0.0))
c     Original table min is 0 for yp
c.OR.(lyp<lypminin)
     1     .OR.((lyp>lypmaxout).OR.(lyp<lypminout))
     1     .OR.((ltk>ltkmaxin))
c     Original talbe min is 0 for Tk
c.OR.(ltk<ltkminin)
     1     .OR.((ltk>ltkmaxout).OR.(ltk<ltkminout))
     1     ) then
c     
c
c
c     Note that want to not use Shen table if
c     a) Temperature is larger than in Shen table
c     b) Density is below Shen table
c     c) Yp is larger than Shen table
      didconverge=0
      else
         didconverge=1
      end if




      





      return
      end




















cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate from Shen table to desired values
c
c     This function based upon get_eos_fromlookup() in HARM code
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interpolate_shen(whichfunction, rhobi, tki, ypi, shen_row, shen_sing, funinterptype, funfloor, badinterp)
      implicit none

c     Get Shen table size for ???_row array below
      include 'vector_sheneos.dek'


c     Passed parameters
      integer whichfunction
      integer funinterptype
      double precision rhobi,tki,ypi
      double precision shen_row(numt,numyp,numrhob)
      double precision funfloor

c     To be returned parameters
      double precision shen_sing
      integer badinterp

      integer interptype

c     0 = nearest neighbord
c     1 = linear per dimension
c     2 = quadratic per dimension (not yet)
      interptype=1



      if(interptype.eq.0) then
         call interpolate_nearest_shen(whichfunction, rhobi, tki, ypi, shen_row, shen_sing, funinterptype, funfloor, badinterp)
      else if(interptype.eq.1) then
         call interpolate_linear_shen(whichfunction, rhobi, tki, ypi, shen_row, shen_sing, funinterptype, funfloor, badinterp)
      else if(interptype.eq.2) then
      end if


      return
      end















cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nearest Neighbord Interpolate from Shen table to desired values
c
c     This function based upon get_eos_fromlookup() in HARM code
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interpolate_nearest_shen(whichfunction, rhobi, tki, ypi, shen_row, shen_sing, funinterptype, funfloor, badinterp)
      implicit none

c     Get Shen table size for ???_row array below
      include 'vector_sheneos.dek'

c     Passed parameters
      integer whichfunction
      integer funinterptype
      double precision rhobi,tki,ypi
      double precision shen_row(numt,numyp,numrhob)
      double precision funfloor

c     To be returned parameters
      double precision shen_sing
      integer badinterp

      integer iii,jjj,kkk
      integer ii,jj,kk

      double precision funvalue
      integer waswithintable


      ii=idint(rhobi) ! fastest index in _row
      jj=idint(tki)   ! slowest index in _row
      kk=idint(ypi)   ! middle index in _row
	


      waswithintable=0
c     iii,jjj,kkk are offsets from real values, so start at 0 even in Fortran by subtracting 1 below
      do iii=1,2
         do jjj=1,2
            do kkk=1,2

               if(shenwithintable_row(jj+jjj-1,kk+kkk-1,ii+iii-1).eq.1) then

                  waswithintable=1

                  if(funinterptype.eq.0) then
                     funvalue = shen_row(jj+jjj-1,kk+kkk-1,ii+iii-1)
                  else if(funinterptype.eq.1) then
                     if(shen_row(jj+jjj-1,kk+kkk-1,ii+iii-1)<funfloor) then
                        funvalue = log10(funfloor);
                     else
                        funvalue = log10(shen_row(jj+jjj-1,kk+kkk-1,ii+iii-1)+funfloor)
                     end if
                  end if

                  goto 888
                  
               else
c     Not in valid table range
               end if

            end do
         end do
      end do
c     end loop over dimensions


 888  continue



cccccccccccccccccccccccccccccccccccc
      if(waswithintable.eq.0) then

         badinterp=1
         write(*,*) 'Never found good value',rhobi, tki, ypi

      else
         badinterp=0
         shen_sing = funvalue

         if(funinterptype.eq.1) then
            shen_sing=10**(shen_sing)
         end if
      end if





      return
      end





cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Linearly Interpolate from Shen table to desired values
c
c     This function based upon get_eos_fromlookup() in HARM code
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interpolate_linear_shen(whichfunction, rhobi, tki, ypi, shen_row, shen_sing, funinterptype, funfloor, badinterp)
      implicit none

c     Get Shen table size for ???_row array below
      include 'vector_sheneos.dek'

c     Passed parameters
      integer whichfunction
      integer funinterptype
      double precision rhobi,tki,ypi
      double precision shen_row(numt,numyp,numrhob)
      double precision funfloor

c     To be returned parameters
      double precision shen_sing
      integer badinterp

      integer ii,jj,kk
      double precision di(2),dj(2),dk(2)
      double precision dist(2,2,2),f(2,2,2)
      double precision totaldist,totalf
      integer iii,jjj,kkk

c     Below used to detect issues with data
      double precision negative1,mynan,myinf

c     Used to indicate if got good values
      integer totalsummed

c     Note that Matlab and Fortran agree on text appearance of NaN, so reading in NaN is used to indicate when table is out of range and don't want to use that data (like shifting away from bad regions)
c      negative1=-1.0
c      mynan=sqrt(negative1)

      negative1=0.0
      mynan=negative1/negative1
      myinf=1.0/0.0


      ii=idint(rhobi) ! fastest index in _row
      jj=idint(tki)   ! slowest index in _row
      kk=idint(ypi)   ! middle index in _row
	
c     distance from ii to rhobi varying from 0..1
      di(2)=rhobi-dble(ii)
      di(1)=1.0-di(2)
      dj(2)=tki-dble(jj)
      dj(1)=1.0-dj(2)
      dk(2)=ypi-dble(kk)
      dk(1)=1.0-dk(2)

    
c     Loop over nearby table values and determine tri-linearly interpolated value
c     3-D table means 2^3=8 positions
      totaldist=0.0
      totalf=0.0
      totalsummed=0

c     DEBUG:
c      if(whichfunction.eq.16) then
c     write(*,*) 'fdist',f(iii,jjj,kkk),totalf,dist(iii,jjj,kkk),totaldist
c         write(*,*) 'di',di(1),di(2),dj(1),dj(2),dk(1),dk(2)
c      end if

c     DEBUG:
c      write(*,*) ' \n'

c     iii,jjj,kkk are offsets from real values, so start at 0 even in Fortran by subtracting 1 below
      do iii=1,2
         do jjj=1,2
            do kkk=1,2

               if(shenwithintable_row(jj+jjj-1,kk+kkk-1,ii+iii-1).eq.1) then
c     Get distance and function value
                  dist(iii,jjj,kkk) = di(iii)*dj(jjj)*dk(kkk)

                  if(dist(iii,jjj,kkk).lt.0.0) then
                     write(*,*) 'dist was negative',ii,jj,kk,iii,jjj,kkk,dist(iii,jjj,kkk)
                  end if

                  if(funinterptype.eq.0) then
                     f(iii,jjj,kkk) = shen_row(jj+jjj-1,kk+kkk-1,ii+iii-1)
                  else if(funinterptype.eq.1) then
                     if(shen_row(jj+jjj-1,kk+kkk-1,ii+iii-1)<funfloor) then
                        f(iii,jjj,kkk) = log10(funfloor);
                     else
                        f(iii,jjj,kkk) = log10(shen_row(jj+jjj-1,kk+kkk-1,ii+iii-1)+funfloor)
c     DEBUG:
c     write(*,*) jj,kk,ii,jjj,kkk,iii,shen_row(jj+jjj-1,kk+kkk-1,ii+iii-1),funfloor
                     end if
                  end if
                  
c     Then do add
                  totalsummed=totalsummed+1
c     Here dist is really a weight for that value.
c     For the value to be important dist=1 that implies true distance is actually 0 from that point
                  totaldist = totaldist + dist(iii,jjj,kkk)
                  totalf = totalf + f(iii,jjj,kkk)*dist(iii,jjj,kkk)

c     DEBUG:
c                  if(whichfunction.eq.16) then
c                     write(*,*) 'fdist',f(iii,jjj,kkk),dist(iii,jjj,kkk)
c                  end if
                  
                  if (f(iii,jjj,kkk)/=f(iii,jjj,kkk)) then
                     write(*,*) 'Trapped NaN'
                  end if
                  
c     DEBUG:
c     write(*,*) 'add','totalf=',totalf,'f=',f(iii,jjj,kkk)
               else
c     Don't add then
c     DEBUG:
c                  write(*,*) 'noadd','totalf=',totalf,'f=',f(iii,jjj,kkk)
               end if

            end do
         end do
      end do
c     end loop over dimensions

c 0.2381685551976161E+00007   0.1122667773510816E+00013
cccccccccccccccccccccccccccccccccccc
c
c     finally normalize and check solution type
c
cccccccccccccccccccccccccccccccccccc
      if(totalsummed.eq.0) then
c     just use nearest neighbor if no valid inversion (NO!)
c     shen_sing = shen_row(jj,kk,ii)

         badinterp=1
         write(*,*) 'Never found good value',rhobi, tki, ypi
         

      else
c     Good interpolation, so normalize
         badinterp=0
         shen_sing =totalf/totaldist

         if(funinterptype.eq.1) then
            shen_sing=10**(shen_sing)
         end if
      end if


c     DEBUG:
c      if(whichfunction.eq.16) then
c         write(*,*) 'fdisttot',totalf,totaldist
c      end if
                  




 100  format(40(E27.16E5,' '))

      return
      end



















cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Convert from Shen to LS so can use LS storage routines for rest of quantities
c
ccccccc NOT:
c ccccccCompare this function with LS version in jon_lsbox.f called:
c ccccccstore_row_fromlseos2cgs()
c
c Instead of making like cgsnuc quantities, just call store_row_fromlseos2cgs() after making like *directly* LSEOS quantities
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine store_row_fromsheneos2lseos(density_cgs,
     1     temperature_cgs, density_cgs_lseos, temperature_cgs_lseos, yp, didconverge) 
      implicit none

c     Passed quantity
c     Local index
      integer didconverge

c     Passed quantities
      double precision density_cgs,temperature_cgs,yp
      double precision density_cgs_lseos,temperature_cgs_lseos


c     Below defined in vector_sneos.dek
      double precision tin,yein,dnsin

      include 'const.dek'
c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
      include 'el_eos.inc'
c  JCM:
c  contains single globals like sound, dse, etc.
      include 'vector_eos.single.dek'
      include 'vector_eos.dek'
c Below contains _cgs and _nuc globals
      include 'vector_sneos.dek'

c     Bring in Shen table parameters
      include 'vector_sheneos.dek'


      
c      write(*,*) 'shenmunp',shenmun_sing,shenmup_sing


      if(didconverge.eq.1) then

cccccccccccccccccccccccccccccccccccccc
c
c     Convert to nuclear units before calling any_electron()
c
c     Inside any_electron() tin and dnsin are used or not depending upon usediffinputs
c
cccccccccccccccccccccccccccccccccccccc
         tin = temperature_cgs_lseos*k2mev
         yein=yp
         dnsin=density_cgs_lseos*(avo*fm3)


ccccccccccccccccccccccccccccccccc
c     
c     Convert from Shen ???_sing to LS EOS single global structures that are in vector_sneos.dek and eos_m4c.commononly.inc and maybe el_eos.inc
c     
cccccccccccccccccccccccccccccccc

c     shenlrhob_sing
c     shennb_sing
c     shenlyp_sing
c     shenyp_sing
c     These are in nuclear units
         bftot=shenf_sing
         bu=shenebulk_sing
         bs=shensbulk_sing
         bpress=shenpbulk_sing

c     These are dimensionless
         a=shenaheav_sing
         x=shenzheav_sing/shenaheav_sing
c     Below not the same
         muhat=shenmstar_sing
         xnut=shenxneut_sing
         xprot=shenxprot_sing
         xalfa=shenxalfa_sing
         xh=shenxh_sing


c         write(*,*) 'azheav',a,x,xnut,xprot,xalfa,xh

c     These are in nuclear units (GODMARK: what to do with these?)
c         shenmun_sing
c         shenmup_sing
c     shenltemp_sing
c     shentemp_sing

c     These are actually relative to M, so are really \hat{\eta}
         etapls=shenmup_sing/tin
c     shenmunminusmup_sing = \mu_n - \mu_p, so to get \eta_n add back in \mu_p
c     Did this for accuracy of interpolation for \eta_n-\eta_p used often
         etanls=shenmunminusmup_sing/tin+etapls


c         write(*,*) 'rhotemp',density_cgs_lseos,temperature_cgs_lseos
c         write(*,*) 'etanls=',etanls,etapls
         
c     Calls any electron EOS and puts electron EOS quantities into nuclear form into direct LS EOS type quantities
c     abar,zbar,zbarnum will be recomputed in any_electron()
         call any_electron(tin,yein,dnsin)

c     any_electron() generates:
c     1) abarnum,abar,zbar     : species
c     2) nusbe,neplus,musube   : output vector
c     3) epress, eu, es, fsube : electrons
c     4) ppress,pu,ps,pf       : photons
c     5) demudt,demudn,demudy  : chemical potential derivatives
c     6) depdt,depdn,depdy,deudt,deudn,deudy,desdt,desdn,desdy : electron derivatives
c     7) dppdt, dppdn, dppdy, dpudt, dpudn, dpudy, dpsdt, dpsdn, dpsdy : photon derivatives

c     LSEOS processes these and generates:
         ftot = bftot  + fsube  + pf
         utot = bu     + eu     + pu
         stot = bs     + es     + ps
         ptot = bpress + epress + ppress
         
c     GODMARK:
c     Shen doesn't have derivatives, so just set as 0 for now and add up LSEOS type terms
         dbsdt = 0.0
         dsdt  = dbsdt + desdt + dpsdt
         dbsdn = 0.0
         dsdn  = dbsdn + desdn + dpsdn
         dbsdy = 0.0
         dsdy  = dbsdy + desdy + dpsdy

         dbudt = 0.0
         dudt  = dbudt + deudt + dpudt
         dbudn = 0.0
         dudn  = dbudn + deudn + dpudn
         dbudy = 0.0
         dudy  = dbudy + deudy + dpudy

         dbpdt = 0.0
         dpdt  = dbpdt + depdt + dppdt
         dbpdn = 0.0
         dpdn  = dbpdn + depdn + dppdn
         dbpdy = 0.0
         dpdy  = dbpdy + depdy + dppdy
         
c         dbmudt = 0.0
c         dmudt = dbmudt + ye*demudt
c         dbmudn = 0.0
c         dmudn = dbmudn + ye*demudn
c         dbmudy = 0.0
c         dmudy = dbmudy + ye*demudy



      else

c     Set things to 0 as indicator that did not converge
         call zero_lseos_quantities()

      end if





      return
      end







cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Zero-out LSEOS quantities
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zero_lseos_quantities()
      implicit none

      include 'const.dek'
c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
      include 'el_eos.inc'
c  JCM:
c  contains single globals like sound, dse, etc.
      include 'vector_eos.single.dek'
      include 'vector_eos.dek'
c Below contains _cgs and _nuc globals
      include 'vector_sneos.dek'



c     Pressure

      ppress=0.0
      epress=0.0
      bpress=0.0
      ptot=0.0
      dpdt=0.0
      dpdn=0.0
      dpdy=0.0

c     Entropy per baryon

      ps=0.0
      es=0.0
      bs=0.0
      stot=0.0
      dsdt=0.0
      dsdn=0.0
      dsdy=0.0


c     energy per baryon

      pu=0.0
      eu=0.0
      bu=0.0
      utot=0.0
      dudt=0.0
      dudn=0.0
      dudy=0.0


c     Free energy

      pf=0.0
      fsube=0.0
      bftot=0.0
      ftot=0.0


      etanls=0.0
      etapls=0.0


      return
      end





