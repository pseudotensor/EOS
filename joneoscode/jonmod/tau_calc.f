c TODO: Something wrong with:
c     Qmphoton too small compared to before
c
c
c     Assumes species computed via kaz_species_etae() or externally
c
c     Globally uses several quantities, see preparecall2kazeos() in jon_lsbox.f
c     Things globally need include (e.g.):
c     rhob,tk,hcm,tdynorye,tdynorynu
c     npfree,nnfree,npheav,nnheav,kazaheav,kazzheav
c     yetot,abarbound


c     Notice that \nu_e and \bar{\nu}_e are treated differently.
c     Notice that \nu_\mu and \nu_\tau are treated the same
c     Notice that \bar{\nu}_\mu and \bar{\nu}_\tau are treated the same

c     Note that when finding rho,T,Y_e,Y_\nu table, only things computed here change when changing Y_\nu.

c     ynumethod=0 : quasi-thermalization (must have hcm dep)
c     ynumethod=1 : Y_nu dependence (so no hcm dep)
c     ynumethod=2 : Y_nu thermal (so no hcm dep)
c     ynumethod=3 : like ynumethod==0, except eta_nue = - eta_nuebar and E_nu's are determined like ynumethod.eq.1


C======================================================================
      subroutine tau_calc(ynumethod,rho10,T11,xnuc,etae,etap,etan,etanu)
      
C======================================================================
      implicit none

      include 'kazeos1.dek'


c     Passed:
      integer ynumethod
      integer tauiter
      real*8 rho10,T11,etae,etap,etan,etanu
      real*8 xnuc

c     Locals:
      integer tempynumethod
      real*8 etanuthermal
      real*8 SMALL
      integer NUMTAUITER
      real*8 ETAERROR,ETASUPERSMALL,YNUERROR
      real*8 MINETA,MINYNU
      real*8 etanuold,Ynu0old,Ynuold,error1,error2,error3
      real*8 eta_nuemin,eta_nuemax,ferr,fmin,fmax
      real*8 eps1,eps2
c     Locals: (now use save in sub function)
c      real*8 eta_nue,eta_nuebar
c      real*8 qtaut_nue,qtaut_nuebar,qtautmu,qtauttau
c      real*8 ntaut_nue,ntaut_nuebar,ntautmu,ntauttau


c     A small number to avoid division by 0
      SMALL = 1D-150





 
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) 'getting thermal'


c     First get completely thermalized Ynu

c     Below is fake way to set extrablock->1
      Ynuthermal=tdynorynu
c     \beta-equilibrium
      etanuthermal = etae + etap - etan
      etanu = etanuthermal
c     Must tell tau_calc_sub() that not ymumethod==0 where etanu depends upon H and so would get reset near end of its call.  Doesn't happen to currently affect this calculation, but affects what etanu is fed into rest of calculations.
      tempynumethod=2
c     Do not iterate, just solve and force tauiter==0 so pure thermalization
      tauiter=0
      call tau_calc_sub(tauiter,tempynumethod
     1     ,rho10,T11,xnuc,etae,etap,etan,etanu
     1     )
c     Assign from one global to another
c     Ynuthermal for method that includes H in table
c     AND Ynuthermal is used to set blocking factor as relative to Ynu0 *always*.  Never use optical depth version.
      Ynuthermal0=Ynu0
      Ynuthermal=Ynu
c     Need 0 version of number densities to compute Ynuthermal in HARM and below recomputation of Ynuthermal
      nnueth0=nnue0
      nnuebarth0=nnuebar0


      
c      write(*,*) 'Ynuthermal',rho10,T11,etae,etap,etan,etanu,Ynuthermal



cccccccccccccccccccccccccccc
c
c     Assume complete thermalization (if ynumethod.eq.2, then first get to Ynuthermal is sufficient)
c
c     ynumethod.eq.2
c
cccccccccccccccccccccccccccc

      if(ynumethod.eq.2) then
c     Uses optical depth dupressed Ynuthermal
         Ynu = Ynuthermal
         tdynorynu = Ynuthermal
      end if


cccccccccccccccccccccccccccc
c
c     Iterate to obtain converged solution for \eta_\nu for neutrinos
c
c     ynumethod.eq.0  .OR. ynumethod.eq.3
c
cccccccccccccccccccccccccccc
      if(ynumethod.eq.0 .OR. ynumethod.eq.3) then

c     Set number of iterations (if ==1, then input etanu only used)
c     NUMTAUITER=1
c     Don't need many iterations to get reasonable result
         NUMTAUITER=10
c     Generally not more accurate than 10% anyways, so ok to seek error of 1%
         ETAERROR=1D-2
         ETASUPERSMALL=SMALL
c     Don't worry as much about relative error if quantities are below the following values:
         MINETA=1D-5
         MINYNU=SMALL
c     Start out with thermal as guess
         etanu=etanuthermal
         Ynu0=Ynuthermal0
         Ynu=Ynuthermal
c     Iterate
c     NOTE that iteration really iterates over not only etanu but also Ynu,Ynu0 because of blocking factor.
         do tauiter=1,NUMTAUITER
            etanuold=etanu
            Ynu0old=Ynu0
            Ynuold=Ynu
            call tau_calc_sub(tauiter,ynumethod
     1           ,rho10,T11,xnuc,etae,etap,etan,etanu
     1           )


c     DEBUG:
c            write(*,*) 'iter etanu',tauiter,rho10,T11,etae,etap,etan,etanu,etanuold,Ynu,Ynu0,Ynuiter


c     Check if converged if gone through at least one full iteration:
            if(tauiter.gt.1) then
               error1=dabs(etanu-etanuold)/(MINETA + dabs(etanu)+dabs(etanuold)+ETASUPERSMALL)
               error2=dabs(Ynu0-Ynu0old)/(MINYNU + dabs(Ynu0)+dabs(Ynu0old)+ETASUPERSMALL)
               error3=dabs(Ynu-Ynuold)/(MINYNU + dabs(Ynu)+dabs(Ynuold)+ETASUPERSMALL)
               if(error1<ETAERROR .AND. error2<ETAERROR .AND. error3<ETAERROR) then
                  goto 777
               end if
            end if

            if(tauiter.eq.NUMTAUITER) then
               write(*,*) 'Never found solution to etanu',rho10,T11,etae,etap,etan,etanu
            end if



         end do
 777     continue

c     Assign final converged Ynu to overwrite "guess" for what is generally an independent variable
         tdynorynu =  Ynuiter

      end if



      
cccccccccccccccccccccccccccc
c
c     Bisect to obtain converged solution for \eta_\nu for given fixed Y_\nu for neutrinos
c     Notice that tau_calc_sub() sets Ynu globally
c
c     ynumethod.eq.1
c
cccccccccccccccccccccccccccc

      if(ynumethod.eq.1) then

c     Set number of iterations (if ==1, then input etanu only used)
c     NUMTAUITER=1
c     Don't need many iterations to get reasonable result
         NUMTAUITER=30
c     Generally not more accurate than 10% anyways, so ok to seek error of 1%
         ETAERROR=1D-2
         YNUERROR=1D-5


c     Set limits of bisection
         eta_nuemin=1D-20
         eta_nuemax=1D20
            

c     Determine fmin and fmax (inside tau_calc_sub())

c         write(*,*) 'getting fmin'

ccccccccccccccccccccccccccccccccccccc
c     Determine fmin
         etanu = eta_nuemin
         tauiter=1
         call tau_calc_sub(tauiter,ynumethod
     1        ,rho10,T11,xnuc,etae,etap,etan,etanu
     1        )

         fmin = (Ynuiter - tdynorynu)/(dabs(tdynorynu)+SMALL)


       if(fmin.gt.0.0) then
            write(*,*) 'fmin is >0.0',fmin,eta_nuemin
            write(*,*) 'Ynuiter,tdynorynu',Ynuiter,tdynorynu
            write(*,*) 'fminmax',eta_nuemin,fmin,eta_nuemax,fmax

c     Try a negative \eta_nu for this case
            etanu = -(etae + etap - etan)*0.95d0
            tauiter=1
            call tau_calc_sub(tauiter,ynumethod
     1           ,rho10,T11,xnuc,etae,etap,etan,etanu
     1           )
            fmin = (Ynuiter - tdynorynu)/(dabs(tdynorynu)+SMALL)
            if(fmin.gt.0.0) then
               write(*,*) 'fmin is >0.0',fmin,eta_nuemin
               write(*,*) 'Ynuiter,tdynorynu',Ynuiter,tdynorynu
               write(*,*) 'Could not fix fmin1'
c     So just use the smallest anti-thermal value obtained rather than stopping -- user should check that generated Ynuiter is not too large
               return ! since done and will just use this anti-thermal etanu and resulting Ynuiter
            end if
         end if



         
c         write(*,*) 'getting fmax'
            
ccccccccccccccccccccccccccccccccccccc
c     Determine fmax
         etanu = eta_nuemax
         tauiter=1
         call tau_calc_sub(tauiter,ynumethod
     1        ,rho10,T11,xnuc,etae,etap,etan,etanu
     1        )
         fmax = (Ynuiter - tdynorynu)/(dabs(tdynorynu)+SMALL)

cfmin is >0.0   121177.009637179       9.999999999999999E-021
c Ynuiter,tdynorynu   15680962086.9804       9.999999939225290E-009
c fminmax  9.999999999999999E-021   121177.009637179       1.000000000000000E+020  1.568096218228100E+018


         if(fmax.lt.0.0) then
            write(*,*) 'fmax is <0.0',fmax,eta_nuemax
            write(*,*) 'Ynuiter,tdynorynu',Ynuiter,tdynorynu
            write(*,*) 'fminmax',eta_nuemin,fmin,eta_nuemax,fmax
c     Rather than stopping, use the Ynuiter from this largest \eta_nu -- user should check Ynuiter at end
            return
         end if
         
c         write(*,*) 'fminmax',eta_nuemin,fmin,eta_nuemax,fmax

c          fmin is >0.0   1.00000000000000       9.999999999999999E-021
c     61  fminmax  9.999999999999999E-021   1.00000000000000       1.000000000000000E+020
c     62  -4.394342758031380E+025



ccccccccccccccccccccccccccccccccccccc
c     Now bisect
         do tauiter=1,NUMTAUITER

c            write(*,*) 'tauiter',tauiter

            etanu = sqrt(eta_nuemin*eta_nuemax)

            call tau_calc_sub(tauiter,ynumethod
     1           ,rho10,T11,xnuc,etae,etap,etan,etanu
     1           )


c     error function
            ferr = (Ynuiter - tdynorynu)/(dabs(tdynorynu)+SMALL)

            if (ferr.gt.0.d0) then
               eta_nuemax=etanu
               fmax=ferr
            else
               eta_nuemin=etanu
               fmin=ferr
            endif

c            write(*,*) 'ferr',ferr,eta_nuemin,etanu,eta_nuemax
c            write(*,*) 'Ynuiter',Ynuiter,tdynorynu

            if(dabs(tdynorynu).gt.SMALL) then
c     Normally want Ynuiter to be close to desired answer
               eps1=dabs(ferr)
               eps2=dabs(log10(eta_nuemax/eta_nuemin))/dabs(log10(eta_nuemin))

               if(eps1<YNUERROR .AND. eps2<ETAERROR) then
                  goto 778
               end if
            else
c     if ynu~0, then seek convergence on \eta_\nu
               eps1=0.0
               eps2=dabs(log10(eta_nuemax/eta_nuemin))/dabs(log10(eta_nuemin))

               if(eps2<ETAERROR) then
                  goto 778
               end if

            end if
            
c     Check if converged if gone through at least one full iteration:

            if(tauiter.eq.NUMTAUITER) then
               write(*,*) 'Never found solution to ynu',rho10,T11,etae,etap,etan,etanu,ynu
            end if


         end do
 778     continue
      end if









c     DEBUG:
c      write(*,*) 'eta=',etae,etap,etan,etanu


      return
      end




C======================================================================
      subroutine tau_calc_sub(tauiter,ynumethod
     1     ,rho10,T11,xnuc,etae,etap,etan,etanu
     1        )
      
C======================================================================
      implicit none
      save ! needed for things only need to compute once
 
C=====================================================================

      include 'kazeos1.dek'
      include 'kazeos.parms.dek'

c      common/tau/tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus
c      common/tau2/tauael_ele, tauael_antiele
c      common /variable1/ tdyn,hcm,rhob,tk


C======================================================================
c      real*8 tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus

C======================================================================
c      real*8 tdyn               !dynamical timescale in sec
c      real*8 hcm                !disk half thickness in cm
c      real*8 rhob               !baryon denstiy in g/cc
c      real*8 tk                 !temperature in K

C======================================================================

      real*8 rate_2stream,density_2stream,tdifffromlambda ! functions

c     Passed:
      integer tauiter,ynumethod
      real*8 rho10,T11,etae,etap,etan,etanu
      real*8 xnuc
c     Below now passed and returned
c      real*8 eta_nue,eta_nuebar

c     Locals:
      real*8 T10,rho10xnuc,nb
      
      real*8 KRrhoscatt,tauphotonscattohcm,tauphotonscatt,tauphotonabs
      real*8 rhoblocal,tklocal,KRrho,tauphoton
      real*8 rhoA10
      real*8 Qphoton0
      real*8 u_photon0,n_photon0


C======================================================================
c     Functions:
      real*8 qdotne,qdotn,qdotpe, qdotAe, qdotpair_ele, qdotpair_mutau,qdotbrem
      real*8 qdotplasmonele, qdotplasmonmutau
      real*8 ndotne,ndotn,ndotpe, ndotAe, ndotpair_ele, ndotpair_mutau,ndotbrem
      real*8 ndotplasmonele, ndotplasmonmutau
C======================================================================
c     Locals:
      real*8 qminusne,qminuspe,qminusAe,qminusn,qminusplasmon,qminuspair,qminusbrem
      real*8 qminusHe,qminusALPHAe
      real*8 qminus_nue,qminus_nuebar
      real*8 nminusne,nminuspe,nminusAe,nminusn,nminusplasmon,nminuspair,nminusbrem
      real*8 nminusHe,nminusALPHAe
      real*8 nminus_nue,nminus_nuebar
c======================================================================
c     Locals:
      real*8 npluspnu,npluspenu,nplusnnu
c======================================================================
c     Functions:
      real*8 tauabs, tauscattNP, tauscattA, tauscattenue,tauscattenuebar,tauscattenu_mutau
      real*8 tausN,tausemu,tausetau
C======================================================================
c     Functions:
      real*8 fermidirac,fermidiraca
C======================================================================
c     Locals:
      real*8 H
      real*8 tmev               !temperature in MeV
C======================================================================
c     Locals:
      real*8 u_nue0,u_nuebar0,u_nuel0,u_numu0,u_nutau0
      real*8 n_nue0,n_nuebar0,n_nuel0,n_numu0,n_nutau0
      real*8 u_nue,u_nuebar,u_nuel,u_numu,u_nutau
      real*8 n_nue,n_nuebar,n_nuel,n_numu,n_nutau
C======================================================================
c     Locals:
      real*8 jj,Q
      real*8 eta_nue_eq,eta_nuebar_eq

c     Below now passed/returned
      real*8 eta_nue,eta_nuebar
      real*8 fd5_nue,fd4_nue,fd3_nue,fd2_nue
      real*8 fd5_nuebar,fd4_nuebar,fd3_nuebar,fd2_nuebar
      real*8 fd50,fd40,fd30,fd20

      real*8 qtaua_nue,qtaua_nuebar,qtauael,qtauamu,qtauatau
      real*8 qtausN_nue,qtausN_nuebar,qtausNnumutau,qtaus
      real*8 qtaus_enue,qtaus_enuebar,qtauseel,qtausemu,qtausetau
      real*8 qtaustot_nue,qtaustot_nuebar,qtaustotel,qtaustotmu,qtaustottau
      real*8 qtautel
c     just outside function now (No: now use save in sub function)
      real*8 qtaut_nue,qtaut_nuebar,qtautmu,qtauttau
      real*8 qtaut_numutau,ntaut_numutau ! temp vars

c     ,ntauael ! global now
      real*8 ntaua_nue,ntaua_nuebar,ntauamu,ntauatau
      real*8 ntausN_nue,ntausN_nuebar,ntausNnumutau,ntaus
      real*8 ntaus_enue,ntaus_enuebar,ntauseel,ntausemu,ntausetau
c     ntaustotel ! global now
      real*8 ntaustot_nue,ntaustot_nuebar,ntaustotmu,ntaustottau
      real*8 ntautel,ntautinel
c     just outside function now
      real*8 ntaut_nue,ntaut_nuebar,ntautmu,ntauttau
      real*8 ntautin_nue,ntautin_nuebar,ntautinmu,ntautintau ! inelastic totals

c     Locals:
      real*8 Elocal
      real*8 qtausNP_nue,qtausNH_nue,qtausNALPHA_nue
      real*8 qtausNP_nuebar,qtausNH_nuebar,qtausNALPHA_nuebar
      real*8 qtausNPnumutau,qtausNHnumutau,qtausNALPHAnumutau

      real*8 ntausNP_nue,ntausNH_nue,ntausNALPHA_nue
      real*8 ntausNP_nuebar,ntausNH_nuebar,ntausNALPHA_nuebar
      real*8 ntausNPnumutau,ntausNHnumutau,ntausNALPHAnumutau

c     Functions:
      real*8 gammannu,gammapnu,gammapenu,gammane,gamman,gammape,gammaAe
      real*8 gammaHe,gammaALPHAe

c     Locals:
      real*8 Qm_nue,Qm_nuebar
      real*8 Nm_nue,Nm_nuebar
c     Locals:
      real*8 npfreenondeg,nnfreenondeg,nA,nAnondeg,nH,nHnondeg
      real*8 nbfree,nbfreenondeg
      real*8 npheavnondeg,nnheavnondeg,npboundnondeg,nnboundnondeg
      real*8 nALPHA,nALPHAnondeg
      real*8 lambdain_nue,lambdain_nuebar,lambdainmu,lambdaintau
      real*8 ilambdaintot
c     ,lambdaintot ! global now

      real*8 gammap2nthin,gamman2pthin

      real*8 taulocal_nue,taulocal_nuebar
      real*8 thermalnpratio
      real*8 SMALL
      real*8 kbtk

      real*8 qtaulocal_nue,qtaulocal_nuebar,qtaulocal_numutau
      real*8 ntaulocal_nue,ntaulocal_nuebar,ntaulocal_numutau

      real*8 qEsq_nue,qEsq_nuebar,qE_nue,qE_nuebar,qEsq_numutau,qE_numutau
      real*8 nEsq_nue,nEsq_nuebar,nE_nue,nE_nuebar,nEsq_numutau,nE_numutau

      real*8 qEsqth_nue,qEsqth_nuebar,qEth_nue,qEth_nuebar,qEsqth_numutau,qEth_numutau
      real*8 nEsqth_nue,nEsqth_nuebar,nEth_nue,nEth_nuebar,nEsqth_numutau,nEth_numutau

      real*8 qEsqnt_nue,qEsqnt_nuebar,qEnt_nue,qEnt_nuebar,qEsqnt_numutau,qEnt_numutau
      real*8 nEsqnt_nue,nEsqnt_nuebar,nEnt_nue,nEnt_nuebar,nEsqnt_numutau,nEnt_numutau

      integer noepos,noeneg
      integer nonue,nonuebar
      real*8 fd2_eneg,fd2_epos,fd3_eneg,fd3_epos,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
      real*8 block_nnu,block_pnu,block_pe,block_ne
      real*8 block_ee_nue,block_ee_nuebar,block_ee_numutau
      real*8 block_ee_nuetot,block_ee_numutautot

      real*8 gammap,block_ge_nue,block_ge_nuebar,block_ge_nuetot,block_ge_numutau,block_ge_numutautot
      real*8 qtau_nnu,qtau_pnu,ntau_nnu,ntau_pnu

      real*8 Ebaryonfactor_nue,Ebaryonfactor_nuebar,Ebaryonfactor_mutau
      real*8 Eelectronfactor_nue,Eelectronfactor_nuebar,Eelectronfactor_mutau

      real*8 abaralpha,zbaralpha,yealpha
      real*8 zbarboundlocal

      real*8 eta_nue_factor,eta_nuebar_factor,VoCthermal_nue,VoCthermal_nuebar,VoCthermal_numutau


      real*8 nlambda_nue,nlambda_nuebar,nlambda_numu,nlambda_nutau
      real*8 qlambda_nue,qlambda_nuebar,qlambda_numu,qlambda_nutau
      real*8 nTdifftot,inlambdatot,nlambdatot
      real*8 nTdiff_nue,nTdiff_nuebar,nTdiff_numu,nTdiff_nutau
      real*8 qTdiff_nue,qTdiff_nuebar,qTdiff_numu,qTdiff_nutau
      real*8 qTdifftot,iqlambdatot,qlambdatot

      real*8 nTlosstot,nTloss_nue,nTloss_nuebar,nTloss_numu,nTloss_nutau
      real*8 qTlosstot,qTloss_nue,qTloss_nuebar,qTloss_numu,qTloss_nutau

      real*8 RufNm_nue,RufNm_nuebar,RufNm_numu,RufNm_nutau
      real*8 RufQm_nue,RufQm_nuebar,RufQm_numu,RufQm_nutau

      real*8 Ynulocal,Ynu0local
      real*8 etanbruenn,etapbruenn
      real*8 etanruffert,etapruffert

      real*8 Ynpfree,Ypnfree
      real*8 Ynpheav,Ypnheav
      real*8 Ynpbound,Ypnbound
      real*8 rho10free,rho10freenondeg

      real*8 Xpnondeg,Xnnondeg

      real*8 Elocalnue,Elocalnuebar,Elocalnumutau

      real*8 extrablock
      real*8 etapwithm,etanwithm,etabwithm,etaewithm


      include 'const.dek'




C======================================================================
C============================= tau-calculation part=====================
C======================================================================

      SMALL = 1D-150
      H = hcm
c      T10=T11/1.d1
      T10=T11*10.0
c      free baryons
      rho10xnuc=rho10*xnuc
      rhoblocal = (rho10*1.0D10)
      nb = rhoblocal/mb
      tklocal = (T10*1.0D10)
c      tmp(1)=xnuc_small**2
c      xnuc = tmp(1)/(1.d0+tmp(1))+xnuc_small*1.d0/(1.d0+tmp(1))
      tmev = T11*1.d11/mev2K
      kbtk=(kerg*tklocal)



      if(tauiter.eq.0) then
c     Perfect thermalization calculation
         extrablock = 1.0
      else
c     Setup Fermi blocking correction to estimate blocking factor better
c     in regimes where neutrinos are (and are not) part of collisional system
         extrablock = dmax1(0.0d0,dmin1(1.0d0,Ynu0/Ynuthermal0))
      end if



cccccccccccccccccccccccccccc
c
c     Set initial neutrino chemical potential
c     
c     If first-time call, then use input etanu, which is probably 0, which is initial guess
c
c     if ynumethod.eq.1, then always reset eta_nue and eta_nuebar to input values
c
c     These iterated-like quantities are saved using the "save" Fortran keyword for the current function call
c
cccccccccccccccccccccccccccc
      if(tauiter.eq.0 .OR. tauiter.eq.1 .OR. ynumethod.eq.1 .OR. ynumethod.eq.2) then
         eta_nue=etanu
c     GODMARK: below requires fast reactions once thermalized (used by ynumethod.eq.1 and 2 to determine \eta_\nu)
         eta_nuebar=-eta_nue
         qtaut_nue=SMALL
         qtaut_nuebar=SMALL
         qtautmu=SMALL
         qtauttau=SMALL
         ntaut_nue=SMALL
         ntaut_nuebar=SMALL
         ntautmu=SMALL
         ntauttau=SMALL
      end if

      if(ynumethod.eq.3) then
         eta_nue=etanu
c     GODMARK: below requires fast reactions once thermalized (used by ynumethod.eq.1 and 2 to determine \eta_\nu)
         eta_nuebar=-eta_nue
      end if




cccccccccccccccccccccccccccc
c
c     Get Fermi-Dirac integrals
c
cccccccccccccccccccccccccccc

c     Test of integrals at high eta
c      eta_nue=.66701606536848172545D+00007
c      eta_nuebar=-eta_nue

c     DEBUG:
c      write(*,*) 'etanus',tauiter,eta_nue,eta_nuebar



c     we often form the ratios fd5/fd3 or fd4/fd2 or use fd2 and fd3 multiplicatively for number/energy densities
      jj=2.0
      Q=0.0
      fd2_nue=fermidirac(jj,eta_nue,Q)
      fd2_nuebar=fermidirac(jj,eta_nuebar,Q)

      jj=3.0
      Q=0.0
      fd3_nue=fermidirac(jj,eta_nue,Q)
      fd3_nuebar=fermidirac(jj,eta_nuebar,Q)

      jj=4.0
      Q=0.0
      fd4_nue=fermidirac(jj,eta_nue,Q)
      fd4_nuebar=fermidirac(jj,eta_nuebar,Q)

      jj=5.0
      Q=0.0
      fd5_nue=fermidirac(jj,eta_nue,Q)
      fd5_nuebar=fermidirac(jj,eta_nuebar,Q)


c     If Fermi_Dirac integral is 0, then energy/number density is 0 and must avoid division by 0.
c     This can happen for very large negative eta_nuebar, so correct those terms
c     For example, number density of anti-neutrinos can be 0 if large enough \eta_nu, so avoid but keep negligible
      if(1.eq.1) then


         nonue=0
         nonuebar=0
         if( (fd2_nue.lt.SMALL).OR.(fd3_nue.lt.SMALL)
     1        .OR.(fd4_nue.lt.SMALL).OR.(fd5_nue.lt.SMALL) ) then
c     Assume no electron neutrinos
            nonue = 1
            fd2_nue=SMALL
            fd3_nue=SMALL
            fd4_nue=SMALL*1D-10
            fd5_nue=SMALL*1D-10
         end if

         if( (fd2_nuebar.lt.SMALL).OR.(fd3_nuebar.lt.SMALL)
     1        .OR.(fd4_nuebar.lt.SMALL).OR.(fd5_nuebar.lt.SMALL) ) then
c     Assume no anti-electron neutrinos
            nonuebar = 1
            fd2_nuebar=SMALL
            fd3_nuebar=SMALL
            fd4_nuebar=SMALL*1D-10
            fd5_nuebar=SMALL*1D-10
         end if
      else
         if(fd2_nuebar.lt.SMALL) then
            fd2_nuebar=SMALL
            fd4_nuebar=0.0d0
         end if
         if(fd3_nuebar.lt.SMALL) then
            fd3_nuebar=SMALL
            fd5_nuebar=0.0d0
         end if

         if(fd2_nue.lt.SMALL) then
            fd2_nue=SMALL
            fd4_nue=0.0d0
         end if
         if(fd3_nue.lt.SMALL) then
            fd3_nue=SMALL
            fd5_nue=0.0d0
         end if
      end if

c     FD for eta=0 for mu and tau neutrinos and anti-neutrinos
      fd20=1.80309
      fd30=5.6822
      fd40=23.3309
      fd50=118.266


c      write(*,*) 'etanue,etanuebar',eta_nue,eta_nuebar
c      write(*,*) 'nonue,nonuebar',nonue,nonuebar


      if(tauiter.eq.0 .OR. tauiter.eq.1) then
c     Only  need to compute this once
c     Fermi-Dirac integrals for electrons for neutrino-related blocking factors
         jj=2.0
         Q=0.0
         fd2_eneg=fermidirac(jj,etae,Q)
         fd2_epos=fermidirac(jj,-etae,Q)

         jj=3.0
         Q=0.0
         fd3_eneg=fermidirac(jj,etae,Q)
         fd3_epos=fermidirac(jj,-etae,Q)

         jj=4.0
         Q=0.0
         fd4_eneg=fermidirac(jj,etae,Q)
         fd4_epos=fermidirac(jj,-etae,Q)

         jj=5.0
         Q=0.0
         fd5_eneg=fermidirac(jj,etae,Q)
         fd5_epos=fermidirac(jj,-etae,Q)


c     Check for abscence of positrons (and electrons) to force correct behavior in that limit
c     Note that energy per neutrino estimated as $1.0$ if those SMALL checks are activated
c     But should imply very small actual rates since number rates don't divide by these SMALL's
         noepos=0
         noeneg=0
         if(1.eq.1) then


            if( (fd2_epos.lt.SMALL).OR.(fd3_epos.lt.SMALL)
     1           .OR.(fd4_epos.lt.SMALL).OR.(fd5_epos.lt.SMALL) ) then
c     Assume no positrons
               noepos = 1
               fd2_epos=SMALL
               fd3_epos=SMALL
               fd4_epos=SMALL*1D-10
               fd5_epos=SMALL*1D-10
            end if

            if( (fd2_eneg.lt.SMALL).OR.(fd3_eneg.lt.SMALL)
     1           .OR.(fd4_eneg.lt.SMALL).OR.(fd5_eneg.lt.SMALL) ) then
c     Assume no electrons
               noeneg = 1
               fd2_eneg=SMALL
               fd3_eneg=SMALL
               fd4_eneg=SMALL*1D-10
               fd5_eneg=SMALL*1D-10
            end if
         end if

      end if


ccccccccccccccccccccccccccccccc
c
c     Compute average THERMAL neutrino energy for energy and number processes
c     These are energies per k_b T
c     Properly accounts for <E_\nu^2> vs. <E_\nu>^2
c
ccccccccccccccccccccccccccccc

c     These are assuming the neutrinos are thermalized
c     KM07 appendix


c     energy processes
c     <E_\nu^2>/(kb T)^2
      qEsqth_nue = fd5_nue/fd3_nue
      qEsqth_nuebar = fd5_nuebar/fd3_nuebar

c     <E_\n>/(kb T)
      qEth_nue = fd4_nue/fd3_nue
      qEth_nuebar = fd4_nuebar/fd3_nuebar

c     <E_\nu^2>/(kb T)^2
      qEsqth_numutau = fd50/fd30
c     <E_\nu>/(kb T)
      qEth_numutau = fd40/fd30


c     number processes
c     <E_\nu^2>/(kb T)^2
      nEsqth_nue = fd4_nue/fd2_nue
      nEsqth_nuebar = fd4_nuebar/fd2_nuebar

c     <E_\nu>/(kb T)
      nEth_nue = fd3_nue/fd2_nue
      nEth_nuebar = fd3_nuebar/fd2_nuebar

c     <E_\nu^2>/(kb T)^2
      nEsqth_numutau = fd40/fd20
c     <E_\nu>/(kb T)
      nEth_numutau = fd30/fd20


c      write(*,*) 'fds_nue',fd4_nue,fd2_nue
c      write(*,*) 'fds_nuebar',fd4_nuebar,fd2_nuebar
cfds  1.000000000000000E-060  1.000000000000000E-050

ccccccccccccccccccccccccccccccc
c
c     Compute neutrino energy and number densities (before optical depth corrections)
c
ccccccccccccccccccccccccccccc


      u_nue0 = 5.82607D-16*tk**4*fd3_nue
      n_nue0 = 4.21965*tk**3*fd2_nue

      u_nuebar0 = 5.82607D-16*tk**4*fd3_nuebar
      n_nuebar0 = 4.21965*tk**3*fd2_nuebar

c     NOTE THAT THE BELOW DOES NOT DEPEND UPON H, so Y_\nu^0
c     NOTE that this is the perfectly thermalized value that can be supressed to 0 or even be made negative!
      Ynu0local = (n_nue0-n_nuebar0)/nb
     
c      write(*,*) 'n_',n_nue0,n_nuebar0,Ynu0local



      u_nuel0=u_nue0 + u_nuebar0 ! diag
      n_nuel0=n_nue0 + n_nuebar0 ! diag


c     for mu (for normal AND anti particles both, hence factor of 2.0)
      u_numu0 = 5.82607D-16*tk**4*fd30*2.0
      n_numu0 = 4.21965*tk**3*fd20*2.0

c     for tau (both normal AND anti particles both)
      u_nutau0 = u_numu0
      n_nutau0 = n_numu0

c     Here factor of 2 used above is to obtain both polarizations of photons and of course photons have 0 chemical potential due to fast electromagnetic interactions
      u_photon0 = u_numu0
      n_photon0 = n_numu0
      

c      write(*,*) 'etanu=',etanu
c      write(*,*) 'fd=',fd2,fd3
c      write(*,*) 'nu=',unu0,nnu0,unu0/nnu0/mev2erg

c     DEBUG:
c      write(*,*) 'eta_nue=',eta_nue
c      write(*,*) 'fd=',fd2_nue,fd3_nue,fd4_nue,fd5_nue
c      write(*,*) 'fd=',fd2_nuebar,fd3_nuebar,fd4_nuebar,fd5_nuebar



ccccccccccccccccccccccccccccccc
c
c     Some number densities -- accounting for nucleon degeneracy
c
ccccccccccccccccccccccccccccc


c     I see that in Ruffert et al. (1996) they express the number density term as \rho Y_{np,pn}, so apparently correct final KM07 expressions
c     Note that Ruffert and Bruenn's definitions are that etan and etap are without their respective rest-masses, so have to remove those from our etan and etap that removed the same mass for both

c      write(*,*) 'free',nnfree,npfree

c     From Bruenn (1985) equation 3.5, his definition of \hat{\mu} is \mu_n - \mu_p without rest-mass subtracted, or at most the same mass subtracted from both
      etanbruenn = etan
      etapbruenn = etap

      etanruffert = etan + (amu-mn)*clight**2/kbtk
      etapruffert = etap + (amu-mp)*clight**2/kbtk

      
c      call compute_old_Ynp_Ypn(yefree,etanbruenn,etapbruenn,Ynpfree,Ypnfree)
c      npfreenondeg = (nnfree+npfree)*Ypnfree
c      nnfreenondeg = (nnfree+npfree)*Ynpfree

c     New method that is easier to fix-up
      call compute_new_Ynp_Ypn(yefree,etanbruenn,etapbruenn,Ynpfree,Ypnfree)
      npfreenondeg = npfree*Ypnfree
      nnfreenondeg = nnfree*Ynpfree

ccccccccccccccccccccccccc
c     GODMARK: Need Sumi on whether below is correct -- how to define npfree from Xnut,Xprot?
ccccccccccccccccccccccccc
c      nbfree=dmax1(1D-150,(npfree*mp+nnfree*mn)/mb)
c See Shen guide just below equation 18 for definition of $n_b$
      nbfree=dmax1(1D-150,(npfree+nnfree))
c     Note that below is a *definition* for the mass-weighting used
c     Presumes below is used really for number of free baryons but in "mass form" so that don't really need true total mass-energy
      rho10free=mb*nbfree/1.0D10
c      nbfreenondeg=dmax1(1D-150,(mp*npfreenondeg+mn*nnfreenondeg)/mb)
      nbfreenondeg=dmax1(1D-150,(npfreenondeg+nnfreenondeg))
      rho10freenondeg=mb*nbfreenondeg/1.0D10

c     Free (and non-degen) neutron fraction
c     npfreenondeg=nfreenondeg=0 is possible when Ynpfree=Ypnfree=0, so limit so Xpnondeg,Xnnondeg non-NaN.  In end will multiply by nbfreenondeg anyways, so sets (e.g.) qdotbrem to 0
      Xpnondeg=npfreenondeg/nbfreenondeg
      Xnnondeg=nnfreenondeg/nbfreenondeg

      if(1.eq.1) then
         npheavnondeg = npheav
         nnheavnondeg = nnheav

         npboundnondeg = npbound
         nnboundnondeg = nnbound
      else
         npheavnondeg = 0.0
         nnheavnondeg = 0.0

         npboundnondeg = 0.0
         nnboundnondeg = 0.0
      end if


c     Apparently this degeneracy calculation only applies to free nucleons (Bruenn 1985)
c      call compute_old_Ynp_Ypn(yeheav,etanbruenn,etapbruenn,Ynpheav,Ypnheav)
c      npheavnondeg = (nnheav+npheav)*Ypnheav
c      nnheavnondeg = (nnheav+npheav)*Ynpheav

c      call compute_old_Ynp_Ypn(yebound,etanbruenn,etapbruenn,Ynpbound,Ypnbound)
c      npboundnondeg = (nnbound+npbound)*Ypnbound
c      nnboundnondeg = (nnbound+npbound)*Ynpbound

c     Above compute_old_Ynp_Ypn() takes care of the below now
c      npfreenondeg = dmax1((nnfree-npfree)/(dexp(etanbruenn-etapbruenn)-1.0d0),0.0d0)
c      nnfreenondeg = dmax1((npfree-nnfree)/(dexp(etapbruenn-etanbruenn)-1.0d0),0.0d0)
c      npheavnondeg = dmax1((nnheav-npheav)/(dexp(etanbruenn-etapbruenn)-1.0d0),0.0d0)
c      nnheavnondeg = dmax1((npheav-nnheav)/(dexp(etapbruenn-etanbruenn)-1.0d0),0.0d0)
c      npboundnondeg = dmax1((nnbound-npbound)/(dexp(etanbruenn-etapbruenn)-1.0d0),0.0d0)
c      nnboundnondeg = dmax1((npbound-nnbound)/(dexp(etapbruenn-etanbruenn)-1.0d0),0.0d0)

c     Above compute_old_Ynp_Ypn() takes care of the below now
c     Limit number density
c     In case etan-etap is smaller than expected
c     Noticed that at low density, high temperature then mun-mup~0 which is why this is needed (essentially to fix numerical error in mup and mun)
c      npfreenondeg = dmin1(npfreenondeg,dmax1(npfree,nnfree))
c      nnfreenondeg = dmin1(nnfreenondeg,dmax1(npfree,nnfree))
c      npheavnondeg = dmin1(npheavnondeg,dmax1(npheav,nnheav))
c      nnheavnondeg = dmin1(nnheavnondeg,dmax1(npheav,nnheav))
c      npboundnondeg =dmin1(npboundnondeg,dmax1(npbound,nnbound))
c      nnboundnondeg =dmin1(nnboundnondeg,dmax1(npbound,nnbound))

c     Heavy number density (A>4)
      if(kazaheav.gt.0.5) then
         nH = ((npheav+nnheav)/kazaheav)
         nHnondeg = ((npheavnondeg+nnheavnondeg)/kazaheav)
      else
         nHnondeg=0.0
      end if

c     DEBUG:
c      write(*,*) 'nH',npheav,nnheav,kazaheav
c      write(*,*) 'nH2',nH,nHnondeg
c      write(*,*) 'nbound',npbound,nnbound



c     alpha particles only A=4
      nALPHA = (rhoblocal/mb)*kazxalfa
c     Apparently alphas can't be degenerately (Fermi block) suppressed?
c      nALPHAnondeg = nALPHA/dmax1(dexp(etanbruenn-etapbruenn),1.0d0)
      nALPHAnondeg = nALPHA
      yealpha = 0.5
      abaralpha = 4.0
      zbaralpha = abaralpha*yealpha

c     All bound (A>1) number density
      if(abarbound.gt.0.5) then
         nA = ((npbound+nnbound)/abarbound)
         nAnondeg = ((npboundnondeg+nnboundnondeg)/abarbound)
      else
         nAnondeg=0.0
      end if










cccccccccccccccccccccccccccc
c
c     Compute photon optical depth and emission rate for Qm
c     Does not require iteration, so only done on tauiter.eq.1 and saved
c     Note that should compute photon pressure, etc. from this too via same method
c     So below only should presently be applied for optically thick case
c
cccccccccccccccccccccccccccc

      if(tauiter.eq.0 .OR. tauiter.eq.1) then
c     KM02 equation 34-36
c     Shapiro & Teukolsky, equation 14.5.55 and Appendix I.44
c     See also equation 14.5.28
c         KRrho = 1.0/(1.0/(4.0D9*yetot*rho10)
c     1        + 1.0/(2.04D4 * yetot*2*kazzbar*rho10**2*T11**(-3.5))
c     1        + 1.0/(1.37243D7 * xnuc*rho10**2*T11**(-3.5))
c     1        )


c     Photon scattering opacity (Thompson)
         KRrhoscatt=4.0D9*yetot*rho10

c     Photon absorption opacity (free-free and bound-free absorption)
c         KRrhoabs = 1.0/(
c     1        + 1.0/(2.04D4 * yetot*2*kazzbar*rho10**2*T11**(-3.5))
c     1        + 1.0/(1.37243D7 * xnuc*rho10**2*T11**(-3.5))
c     1        )

c     Shapiro & Teukolsky equation I.46 computed in ruffert_emissionrates.nb
         rhoA10 = nAnondeg*kazabar*mb/1.0D10
         Qphoton0 = 1.6D46 * rhoA10 * rho10 * T11**(0.5)*yetot*kazzbar**2/kazabar

         tauphotonscattohcm = KRrhoscatt
c         tauphotonabsohcm = KRrhoabs
         tauphotonabsohcm = (1.0d0/H)*dmax1(tauabs(Qphoton0,u_photon0,H),SMALL)

         tauphotonscatt=tauphotonscattohcm*H
         tauphotonabs=tauphotonabsohcm*H

c     total photon optical depth
         tauphoton = tauphotonscatt+tauphotonabs
         tauphotonohcm = tauphoton/hcm
         
c     New version modelled after 2-stream approximation like neutrinos
c     energy rates
         Qphoton=rate_2stream(clight,hcm,u_photon0,tauphoton,tauphotonabs)

c         write(*,*) 'DEBUG1',rhoA10,Qphoton0,kazabar,u_photon0,hcm,tauphoton,tauphotonabs,Qphoton
c         write(*,*) 'DEBUG2',rhoA10,nAnondeg,abarbound,npboundnondeg,nnboundnondeg,etan,etap




c     Don't need number rate or number density of photons

c     Old version:
c     H-independent quantity is $Q_{photon} H^2$
c         Qphoton = (5.6712D-5 * tklocal**4 / tauphotonohcm) /H**2

c     Assume expression (based upon Kirkoff's law) holds in optically thin and thick regimes
c         if(tauphoton.lt.1.0) then
c            write(*,*) 'Oops, optical depth to photons less than unity'
c     stop
c     For now don't compute emissivity
c            Qphoton = 0.0
c         end if
      end if










ccccccccccccccccccccccccccccccc
c
c     Some neutrino phase space blocking factors
c
ccccccccccccccccccccccccccccc

      if(1.eq.0) then
c     These only apply for hot matter and Ruffert form of rates
c     See Ruffert et al. (1996) equations A11-A12 and A15-A16 (for absoroption \beta-process)
c     [not currently used since don't use explicit absorption]
c     e^- + p <-- n + \nu_e
         block_nnu = 1.0/(1.0 + dexp(- (fd5_nue/fd4_nue - etae)))
c     e^+ + n <-- p + \bar{\nu}_e
         block_pnu = 1.0/(1.0 + dexp(- (fd5_nuebar/fd4_nuebar + etae)))

c     See Ruffert et al. (1996) equations B3-B4 (for \beta-process)
c     e^- + p --> n + \nu_e
         block_pe = 1.0/(1.0 + dexp(- (fd5_eneg/fd4_eneg - eta_nue)))
c     e^+ + n --> p + \bar{\nu}_e
         block_ne = 1.0/(1.0 + dexp(- (fd5_epos/fd4_epos - eta_nuebar)))
      else
c     Now using Kaz version of blocking factor directly during integration
         block_nnu = 1.0
         block_pnu = 1.0
         block_pe = 1.0
         block_ne = 1.0
      end if

c     See Ruffert et al. (1996) equations B9
c     e^- + e^+ -> \nu_i + \bar{\nu}_i for each species
      block_ee_nue = 1.0/(1.0 + dexp(- (0.5*fd4_eneg/fd3_eneg + 0.5*fd4_epos/fd3_epos - eta_nue)))
      block_ee_nuebar = 1.0/(1.0 + dexp(- (0.5*fd4_eneg/fd3_eneg + 0.5*fd4_epos/fd3_epos - eta_nuebar)))
      block_ee_nuetot = block_ee_nue*block_ee_nuebar

c     etamu,etatau=0 for both normal and anti-particles (each)
      block_ee_numutau = 1.0/(1.0 + dexp(- (0.5*fd4_eneg/fd3_eneg + 0.5*fd4_epos/fd3_epos - 0.0)))
      block_ee_numutautot = block_ee_numutau**2.0

c     See Ruffert et al. (1996) equations B13
c     \tilde{\gamma} -> \nu_i + \bar{\nu}_i for each species
      
c     See plasmon calculation too since \gamma_p defined there as well
      gammap=3.21295d-2*sqrt(pi**2 + 3.d0*etae*etae)
      block_ge_nue = 1.0/(1.0 + dexp(- (1.0+0.5*gammap**2/(1.0+gammap) - eta_nue)))
      block_ge_nuebar = 1.0/(1.0 + dexp(- (1.0+0.5*gammap**2/(1.0+gammap) - eta_nuebar)))
      block_ge_nuetot = block_ge_nue*block_ge_nuebar

      block_ge_numutau = 1.0/(1.0 + dexp(- (1.0+0.5*gammap**2/(1.0+gammap) - 0.0)))
      block_ge_numutautot = block_ge_numutau**2.0




ccccccccccccccccccccccccccccccc
c
c     Absorption opacities from Ruffert
c     Not to be mixed with absorption opacity using Kirkoff's law
c     Only applies for hot (relativistic) matter, otherwise need to modify Kaz gammapnu,gammannu,gammapenu
c
ccccccccccccccccccccccccccccc

c     Could use Chen & Beloborodov (2007) B2-B3 instead, but then have to integrate
c     Could use Burrows & Thompson (2002), section 3.1 equations 9-10-11
c     Ruffert et al. (2006) equations A11-A12
c     energy rates type
      qtau_nnu = 7.11725D-42*nnfreenondeg*T11**2*fd5_nue/fd3_nue*block_nnu*H
      qtau_pnu = 7.11725D-42*npfreenondeg*T11**2*fd5_nuebar/fd3_nuebar*block_pnu*H
c     number rates type
      ntau_nnu = 7.11725D-42*nnfreenondeg*T11**2*fd4_nue/fd2_nue*block_nnu*H
      ntau_pnu = 7.11725D-42*npfreenondeg*T11**2*fd4_nuebar/fd2_nuebar*block_pnu*H




ccccccccccccccccccccccccccccccc
c
c     Energy emission rates
c
ccccccccccccccccccccccccccccc

c     terms related to gammap2n
      qminuspe=block_pe*qdotpe(T11,etae,etap,etan,eta_nue,extrablock,npfreenondeg)
      qminusHe=block_pe*qdotAe(T11,etae,etap,etan,eta_nue,kazaheav,kazzheav,extrablock,nHnondeg)
      qminusALPHAe=0.0 ! is 0 generally
      qminusAe = qminusHe + qminusALPHAe

c      write(*,*) 'DEBUGQ',tauiter,qminuspe,qminusHe,qminusALPHAe,T11,etae,etap,etan,eta_nue,npfreenondeg,nHnondeg,kazaheav,kazzheav
      
c     Assume blocking factor similar for ne and n
      qminusn=block_ne*qdotn(T11,etae,etap,etan,eta_nuebar,extrablock,nnfreenondeg) ! Should this process be included?  Asked Kaz. GODMARK
      qminusne=block_ne*qdotne(T11,etae,etap,etan,eta_nuebar,extrablock,nnfreenondeg)

      qminusplasmon=block_ge_nuetot*qdotplasmonele(T11,etae)
      qminuspair=block_ee_nuetot*qdotpair_ele(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )
      qminusbrem=qdotbrem(rho10free,rho10freenondeg,T11,Xpnondeg,Xnnondeg)

c     GODMARK: Don't yet have gammannu, gammapnu,gammapenu versions for energy rates

c     creation of electron neutrino
c     0.5's are because those processes involve creating both anti and normal neutrinos
      qminus_nue   =qminuspe+qminusAe + 0.5*qminusplasmon+ 0.5*qminuspair+0.5*qminusbrem
c     creation of electron anti-neutrino
      qminus_nuebar=qminusne+qminusn  + 0.5*qminusplasmon+ 0.5*qminuspair+0.5*qminusbrem

c     for diagnostics
      qminusel = qminus_nue+qminus_nuebar

c     DEBUG:
      if((qminus_nue.lt.0.0).OR.(qminus_nuebar.lt.0.0)) then
         write(*,*) 'fds',fd2_nue,fd2_nuebar,fd3_nue,fd3_nuebar,fd4_nue,fd4_nuebar,fd5_nue,fd5_nuebar
         write(*,*) 'un',u_nue0,u_nuebar0
         write(*,*) 'Ypnfree',nnfree,npfree,Ypnfree
         write(*,*) 'Ynpfree',nnfree,npfree,Ynpfree
         write(*,*) 'ndeg',npfreenondeg,nnfreenondeg,nHnondeg,etan,etap,kazaheav,kazzheav,block_pe
         write(*,*) 'qminus',qminuspe,qminusAe,qminusn,qminusne,qminusplasmon,qminuspair,qminusbrem
         write(*,*) 'qminusel=',qminusel,qminus_nue,qminus_nuebar
      end if



c     For mu and tau particles have different rates
      qminusplasmon=block_ge_numutautot*qdotplasmonmutau(T11,etae)
c     See ruffert_emissionrates.nb for Reenue2Reenux
      qminuspair=block_ee_numutautot*qdotpair_mutau(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )

c     \mu neutrino AND anti-neutrino
      qminusmu=qminusplasmon+qminuspair+qminusbrem

c     \tau neutrino AND anti-neutrino
      qminustau=qminusmu



ccccccccccccccccccccccccccccccc
c
c     Number emission rates (In reality the energy and number rates are simply related by choosing
c     which power index in the Fermi-Dirac integral, but at the moment have as separate type calculations)
c
ccccccccccccccccccccccccccccc

         
c     terms related to gammap2n
      nminuspe = block_pe*npfreenondeg*ddim(gammape(T11,etae,eta_nue,extrablock),SMALL)
c     JCM: below takes into account fact that gammaAe() returns rate per nuclei with atomic mass A
c     gammaAe() n_A = +dn_n/dt contribution and -dn_p/dt contribution
c     Write as per proton in nuclei since each unit converts 1 proton/sec to 1 neutron/sec
      nminusHe   = block_pe*nHnondeg*ddim(gammaAe(T11,etae,etap,etan,eta_nue,kazaheav,kazzheav,extrablock),SMALL)
      nminusALPHAe = 0.0 ! generally 0
      nminusAe = nminusHe + nminusALPHAe

      npluspnu = block_pnu*npfreenondeg*ddim(gammapnu(T11,etae,eta_nuebar,extrablock),SMALL)
c     Assume blocking factor similar for pnu and penu
      npluspenu = block_pnu*npfreenondeg*ddim(gammapenu(T11,etae,eta_nuebar,extrablock),SMALL)

c     Terms related to gamman2p
c     Assume blocking factor similar for ne and n
      nminusn   =  block_ne*nnfreenondeg*ddim(gamman(T11,etae,eta_nuebar,extrablock),SMALL) ! Kaz says this is very small
      nminusne  =  block_ne*nnfreenondeg*ddim(gammane(T11,etae,eta_nuebar,extrablock),SMALL)
      nplusnnu  =  block_nnu*nnfreenondeg*ddim(gammannu(T11,etae,eta_nue,extrablock),SMALL)

      nminusplasmon = block_ge_nuetot*ndotplasmonele(T11,etae)
      nminuspair=block_ee_nuetot*ndotpair_ele(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )

      if(ynumethod.eq.1 .OR. ynumethod.eq.2 .OR. ynumethod.eq.3) then
c     ynumethod.eq.3 here so result as close to ynumethod.1 as possible
         Elocalnue = qEth_nue*kbtk
         Elocalnuebar = qEth_nuebar*kbtk
         Elocalnumutau = qEth_numutau*kbtk

      else
         Elocalnue = (3.15*kbtk)
         Elocalnuebar = Elocalnue
         Elocalnumutau = Elocalnue
      end if


      nminusbrem=ndotbrem(rho10free,rho10freenondeg,T11,etap,etan,Elocalnue,Elocalnuebar,Elocalnumutau,Xpnondeg,Xnnondeg)

c     creation of electron neutrino
      nminus_nue   =nminuspe+nminusAe + 0.5*nminusplasmon+ 0.5*nminuspair+0.5*nminusbrem
c     creation of electron anti-neutrino
      nminus_nuebar=nminusne+nminusn  + 0.5*nminusplasmon+ 0.5*nminuspair+0.5*nminusbrem


c     for diagnostics
      nminusel = nminus_nue + nminus_nuebar

      nminusplasmon=block_ge_numutautot*ndotplasmonmutau(T11,etae)
      nminuspair=block_ee_numutautot*ndotpair_mutau(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )

c     \mu neutrino AND anti-neutrino
      nminusmu=nminusplasmon+nminuspair+nminusbrem

c     \tau neutrino AND anti-neutrino
      nminustau=nminusmu








      if(ynumethod.eq.0) then
ccccccccccccccccccccccccccccccc
c
c     Compute average NON-THERMAL neutrino energy for energy and number processes
c     These are energies per k_b T
c     In reality <E_\nu^2> should be computed differently
c     In reality shouldn't energy and number processes be computed differently?
c
ccccccccccccccccccccccccccccc


c     energy processes
c     <E_\nu>
      qEnt_nue = qminus_nue/nminus_nue/kbtk
      qEnt_nuebar = qminus_nuebar/nminus_nuebar/kbtk

c     <E_\nu^2>
      qEsqnt_nue = (qEnt_nue)**2
      qEsqnt_nuebar = (qEnt_nuebar)**2


c     <E_\nu>
      qEnt_numutau = qminusmu/nminusmu/kbtk
c     <E_\nu^2>
      qEsqnt_numutau = (qEnt_numutau)**2


c     number processes
c     <E_\nu>
      nEnt_nue = qEnt_nue
      nEnt_nuebar = qEnt_nuebar

c     <E_\nu^2>
      nEsqnt_nue = qEsqnt_nue
      nEsqnt_nuebar = qEsqnt_nuebar

c     <E_\nu>
      nEnt_numutau = qEnt_numutau

c     <E_\nu^2>
      nEsqnt_numutau = qEsqnt_numutau




ccccccccccccccccccccccccccccccc
c     
c     Compute average THERMAL - NON-THERMAL neutrino energy for energy and number processes
c     *** These are energies per k_b T ***
c     Compute \mu and \tau even though optical depths can be different (average optical depth per species)
c     
ccccccccccccccccccccccccccccc

c     As in Chen & Beloborodov (2007), consider neutrino energy as linear interpolation between
c     thermal in optically thick case and non-thermal in optically thin case

         
c     These vary from 0 = thin to 1 = thick

c     Tdiff = 3\tau H/c  -> Vdiff = H/Tdiff  -> Vdiff = c/(3\tau)  VoCdiff = 1/(3\tau)

c     new way
         VoCthermal_nue = 1.0/(3.0*qtaut_nue)
         VoCthermal_nue = dmin1(VoCthermal_nue,1.0d0)
         qtaulocal_nue = (1.0-dexp(-qtaut_nue))*(1.0-Vocthermal_nue)

         VoCthermal_nuebar = 1.0/(3.0*qtaut_nuebar)
         VoCthermal_nuebar = dmin1(VoCthermal_nuebar,1.0d0)
         qtaulocal_nuebar = (1.0-dexp(-qtaut_nuebar))*(1.0-Vocthermal_nuebar)

         qtaut_numutau=0.5*(qtautmu+qtauttau)
         VoCthermal_numutau = 1.0/(3.0*qtaut_numutau)
         VoCthermal_numutau = dmin1(VoCthermal_numutau,1.0d0)
         qtaulocal_numutau = (1.0-dexp(-qtaut_numutau))*(1.0-Vocthermal_numutau)


         VoCthermal_nue = 1.0/(3.0*ntaut_nue)
         VoCthermal_nue = dmin1(VoCthermal_nue,1.0d0)
         ntaulocal_nue = (1.0-dexp(-ntaut_nue))*(1.0-Vocthermal_nue)

         VoCthermal_nuebar = 1.0/(3.0*ntaut_nuebar)
         VoCthermal_nuebar = dmin1(VoCthermal_nuebar,1.0d0)
         ntaulocal_nuebar = (1.0-dexp(-ntaut_nuebar))*(1.0-Vocthermal_nuebar)


         ntaut_numutau=0.5*(ntautmu+ntauttau)
         VoCthermal_numutau = 1.0/(3.0*ntaut_numutau)
         VoCthermal_numutau = dmin1(VoCthermal_numutau,1.0d0)
         ntaulocal_numutau = (1.0-dexp(-ntaut_numutau))*(1.0-Vocthermal_numutau)

c     old way
c     qtaulocal_nue = ddim((1.0-dexp(-(qtaut_nue))),0.0d0)
c     qtaulocal_nuebar = ddim((1.0-dexp(-(qtaut_nuebar))),0.0d0)
c     qtaulocal_numutau = ddim((1.0-dexp(-(0.5*(qtautmu+qtauttau)))),0.0d0)

c     ntaulocal_nue = ddim((1.0-dexp(-(ntaut_nue))),0.0d0)
c     ntaulocal_nuebar = ddim((1.0-dexp(-(ntaut_nuebar))),0.0d0)
c     ntaulocal_numutau = ddim((1.0-dexp(-(0.5*(ntautmu+ntauttau)))),0.0d0)


c     energy processes
c     <E_\nu^2>/kbtk**2
         qEsq_nue = qtaulocal_nue*qEsqth_nue + (1.0-qtaulocal_nue)*qEsqnt_nue
         qEsq_nuebar = qtaulocal_nuebar*qEsqth_nuebar + (1.0-qtaulocal_nuebar)*qEsqnt_nuebar

c     DEBUG:
c     write(*,*) 'qEsq_nue',qtaulocal_nue,qEsqth_nue,qEsqnt_nue

c     <E_\nu>/kbtk
         qE_nue = qtaulocal_nue*qEth_nue + (1.0-qtaulocal_nue)*qEnt_nue
         qE_nuebar = qtaulocal_nuebar*qEth_nuebar + (1.0-qtaulocal_nuebar)*qEnt_nuebar

c     <E_\nu^2>/kbtk**2
         qEsq_numutau = qtaulocal_numutau*qEsqth_numutau + (1.0-qtaulocal_numutau)*qEsqnt_numutau
c     <E_\nu>/kbtk
         qE_numutau = qtaulocal_numutau*qEth_numutau + (1.0-qtaulocal_numutau)*qEnt_numutau

c     number processes
c     <E_\nu^2>/kbtk**2
         nEsq_nue = ntaulocal_nue*nEsqth_nue + (1.0-ntaulocal_nue)*nEsqnt_nue
         nEsq_nuebar = ntaulocal_nuebar*nEsqth_nuebar + (1.0-ntaulocal_nuebar)*nEsqnt_nuebar

c     <E_\nu>/kbtk
         nE_nue = ntaulocal_nue*nEth_nue + (1.0-ntaulocal_nue)*nEnt_nue
         nE_nuebar = ntaulocal_nuebar*nEth_nuebar + (1.0-ntaulocal_nuebar)*nEnt_nuebar

c     <E_\nu^2>/kbtk**2
         nEsq_numutau = ntaulocal_numutau*nEsqth_numutau + (1.0-ntaulocal_numutau)*nEsqnt_numutau
c     <E_\nu>/kbtk
         nE_numutau = ntaulocal_numutau*nEth_numutau + (1.0-ntaulocal_numutau)*nEnt_numutau


      else if(ynumethod.eq.1 .OR. ynumethod.eq.2 .OR. ynumethod.eq.3) then

c     Then assume FD distribution is accurate
c     ynumethod.eq.3 is here because need to be as close to final HARM table as possible that uses ynumethod.eq.1.  Elocal is very different for _nue that changes lambdatot that changes everything else if I don't keep this consistent

         qEsq_nue = qEsqth_nue
         qEsq_nuebar = qEsqth_nuebar

         qE_nue = qEth_nue
         qE_nuebar = qEth_nuebar

         qEsq_numutau = qEsqth_numutau
         qE_numutau = qEth_numutau

         nEsq_nue = nEsqth_nue
         nEsq_nuebar = nEsqth_nuebar

         nE_nue = nEth_nue
         nE_nuebar = nEth_nuebar

         nEsq_numutau = nEsqth_numutau
         nE_numutau = nEth_numutau
      end if




ccccccccccccccccccccccccccccccc
c
c     Absorption Optical Depths
c
ccccccccccccccccccccccccccccc

c      if (0.AND. (tmev.lt.memev)) then
c         tauael=tauabs(qminusel,T11,H)*Exp(-memev/tmev)
c         tauamu=tauabs(qminusmu,T11,H)*Exp(-memev/tmev)
c         tauatau=tauabs(qminustau,T11,H)*Exp(-memev/tmev)
c      else
c      end if


c     energy terms
      qtaua_nue=dmax1(tauabs(qminus_nue,u_nue0,H),SMALL)
      qtaua_nuebar=dmax1(tauabs(qminus_nuebar,u_nuebar0,H),SMALL)
      qtauael = dmax1(tauabs(qminusel,u_nuel0,H),SMALL) ! diag
      qtauamu=dmax1(tauabs(qminusmu,u_numu0,H),SMALL)
      qtauatau=dmax1(tauabs(qminustau,u_nutau0,H),SMALL)



c     number terms
      ntaua_nue=dmax1(tauabs(nminus_nue,n_nue0,H),SMALL)
      ntaua_nuebar=dmax1(tauabs(nminus_nuebar,n_nuebar0,H),SMALL)
      ntauael = dmax1(tauabs(nminusel,n_nuel0,H),SMALL) ! diag
      ntauamu=dmax1(tauabs(nminusmu,n_numu0,H),SMALL)
      ntauatau=dmax1(tauabs(nminustau,n_nutau0,H),SMALL)


c     DEBUG:
c      write(*,*) 'preqtaua1',qminus_nue,u_nue0,H
c      write(*,*) 'preqtaua2',qminus_nuebar,u_nuebar0,H
c      write(*,*) 'qtaua',qtaua_nue,qtaua_nuebar,qtauael,qtauamu,qtauatau



ccccccccccccccccccccccccccccccc
c
c     Scattering
c
ccccccccccccccccccccccccccccc


ccccccccccccccccccccccccc
c     energy terms

c     tauscatt only uses rho10, not rho10xnuc since inside takes into account use of xnuc
c     electron normal-neutrinos
      Elocal = qEsq_nue
      qtausNP_nue=tauscattNP(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yefree)
      qtausNH_nue=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yeheav,kazaheav,nHnondeg)
      qtausNALPHA_nue=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yealpha,abaralpha,nALPHAnondeg)
      qtausN_nue = qtausNP_nue + qtausNH_nue + qtausNALPHA_nue

c     DEBUG:
c      write(*,*) qtausNP_nue,qtausNH_nue,qtausNALPHA_nue,qtausN_nue


c     electron anti-neutrinos
      Elocal = qEsq_nuebar
      qtausNP_nuebar=tauscattNP(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yefree)
      qtausNH_nuebar=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yeheav,kazaheav,nHnondeg)
      qtausNALPHA_nuebar=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yealpha,abaralpha,nALPHAnondeg)
      qtausN_nuebar = qtausNP_nuebar + qtausNH_nuebar + qtausNALPHA_nuebar

c     below scattering for *each* mu and tau for BOTH normal- and anti-particles
      Elocal = qEsq_numutau
      qtausNPnumutau=tauscattNP(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yefree)
      qtausNHnumutau=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yeheav,kazaheav,nHnondeg)
      qtausNALPHAnumutau=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yealpha,abaralpha,nALPHAnondeg)
c     Below is averaged over mu and tau so includes BOTH normal and anti-particles
      qtausNnumutau = qtausNPnumutau + qtausNHnumutau + qtausNALPHAnumutau


c     Electron scattering
      qtaus_enue=tauscattenue(qE_nue,rho10,T11,etae,xnuc,H,yetot)
      qtaus_enuebar=tauscattenuebar(qE_nuebar,rho10,T11,etae,xnuc,H,yetot)
c     energy-rate weighted average scattering optical depth to electron type neutrinos
      qtauseel = (qminus_nue*qtaus_enue + qminus_nuebar*qtaus_enuebar)/(qminus_nue+qminus_nuebar+SMALL) ! diag
c     Below is averaged over mu and tau so includes BOTH normal and anti-particles
      qtausemu=tauscattenu_mutau(qE_numutau,rho10,T11,etae,xnuc,H,yetot)
      qtausetau=qtausemu

c     DEBUG:
c      write(*,*) qtaus_enue,qtaus_enuebar,qtauseel



c     Total scattering term per species
      qtaustot_nue = qtausN_nue + qtaus_enue
      qtaustot_nuebar = qtausN_nuebar + qtaus_enuebar
c     energy-rate weighted average scattering optical depth to electron type neutrinos
      qtaustotel = (qminus_nue*qtaustot_nue + qminus_nuebar*qtaustot_nuebar)/(qminus_nue+qminus_nuebar+SMALL) ! diag
      qtaustotmu = qtausNnumutau + qtausemu
      qtaustottau = qtausNnumutau + qtausetau

c     Total scattering optical depth number-rate weighted averaged over all species (just for diagnostics)
      qtaus = (qminusel*qtaustotel + qminusmu*qtaustotmu + qminustau*qtaustottau)/(qminusel+qminusmu+qminustau+SMALL)


c      write(*,*) 'qtausel',qtausN_nue,qtaus_enue,qtausN_nuebar,qtaus_enuebar
c      write(*,*) 'Enu',qE_nue,qE_nuebar,qE_numutau
c      write(*,*) 'qtausmu',qtausNnumutau,qtausemu,qtausetau

ccccccccccccccccccccccccc
c     number terms (just replace qtau->ntau and qE -> nE from above)

c     tauscatt only uses rho10, not rho10xnuc since inside takes into account use of xnuc
c     electron normal-neutrinos
      Elocal = nEsq_nue
      ntausNP_nue=tauscattNP(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yefree)
      ntausNH_nue=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yeheav,kazaheav,nHnondeg)
      ntausNALPHA_nue=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yealpha,abaralpha,nALPHAnondeg)
      ntausN_nue = ntausNP_nue + ntausNH_nue + ntausNALPHA_nue

c      write(*,*) 'ntaus_nue',Elocal,ntausNP_nue,ntausNH_nue,ntausNALPHA_nue,ntausN_nue

cntaus_nue  6.000000000016335E+017   26574.2689412983
c   47842532.8190546        155539.585383746        48024646.6733796
c ntaut  3.150974618575872E-048   48024646.7354004



c     electron anti-neutrinos
      Elocal = nEsq_nuebar
      ntausNP_nuebar=tauscattNP(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yefree)
      ntausNH_nuebar=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yeheav,kazaheav,nHnondeg)
      ntausNALPHA_nuebar=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yealpha,abaralpha,nALPHAnondeg)
      ntausN_nuebar = ntausNP_nuebar + ntausNH_nuebar + ntausNALPHA_nuebar

c     below scattering for *each* mu and tau for BOTH normal- and anti-particles
      Elocal = nEsq_numutau
      ntausNPnumutau=tauscattNP(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yefree)
      ntausNHnumutau=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yeheav,kazaheav,nHnondeg)
      ntausNALPHAnumutau=tauscattA(Elocal,rho10,T11,etae,etap,etan,xnuc,H,yealpha,abaralpha,nALPHAnondeg)
c     Below is really averaged optical depth for including BOTH normal and anti-particles
      ntausNnumutau = ntausNPnumutau + ntausNHnumutau + ntausNALPHAnumutau

c     Electron scattering
      ntaus_enue=tauscattenue(nE_nue,rho10,T11,etae,xnuc,H,yetot)
      ntaus_enuebar=tauscattenuebar(nE_nuebar,rho10,T11,etae,xnuc,H,yetot)
c     number-rate weighted average scattering optical depth to electron type neutrinos
      ntauseel = (nminus_nue*ntaus_enue + nminus_nuebar*ntaus_enuebar)/(nminus_nue+nminus_nuebar+SMALL) ! diag
c     Below is really averaged optical depth for including BOTH normal and anti-particles
      ntausemu=tauscattenu_mutau(nE_numutau,rho10,T11,etae,xnuc,H,yetot)
      ntausetau=ntausemu

c     Total scattering term per species
      ntaustot_nue = ntausN_nue + ntaus_enue
      ntaustot_nuebar = ntausN_nuebar + ntaus_enuebar
c     number-rate weighted average scattering optical depth to electron type neutrinos
      ntaustotel = (nminus_nue*ntaustot_nue + nminus_nuebar*ntaustot_nuebar)/(nminus_nue+nminus_nuebar+SMALL) ! diag
      ntaustotmu = ntausNnumutau + ntausemu
      ntaustottau = ntausNnumutau + ntausetau

c     Total scattering optical depth number-rate weighted averaged over all species (just for diagnostics)
      ntaus = (nminusel*ntaustotel + nminusmu*ntaustotmu + nminustau*ntaustottau)/(nminusel+nminusmu+nminustau+SMALL)




ccccccccccccccccccccccccccccccc
c     
c     Total Absorption + Scattering optical depths for energy and number terms
c
ccccccccccccccccccccccccccccc

     
c     energy terms
      qtaut_nue=qtaua_nue   + qtaustot_nue
      qtaut_nuebar=qtaua_nuebar   + qtaustot_nuebar
      qtautel = (qminus_nue*qtaut_nue + qminus_nuebar*qtaut_nuebar)/(qminus_nue+qminus_nuebar+SMALL) ! for diagnostics only
      qtautmu=qtauamu   + qtaustotmu
      qtauttau=qtauatau + qtaustottau


c     number terms (elastic + inelastic : used preassigned scattering total)
      ntaut_nue=ntaua_nue   + ntaustot_nue
      ntaut_nuebar=ntaua_nuebar   + ntaustot_nuebar
      ntautel = (nminus_nue*ntaut_nue + nminus_nuebar*ntaut_nuebar)/(nminus_nue+nminus_nuebar+SMALL) ! for diagnostics only
      ntautmu=ntauamu   + ntaustotmu
      ntauttau=ntauatau + ntaustottau

c      write(*,*) 'ntaut',ntaua_nue,ntaustot_nue




ccccccccccccccccccccccccccccccc
c     
c     Total Absorption + Scattering optical depths for energy and number terms for INELASTIC terms only (used by thermalization timescale)
c
ccccccccccccccccccccccccccccc

c     number terms (inelastic only : Note separately assign scattering total)
c     Nucleon scattering is only inelastic for E_\nu>>\mu_b[with rest mass] \sim m_b c^2
c     Electron scattering is only inelastic for E_\nu>\mu_e[with rest mass] \sim m_e c^2
c     etap and etan with rest-mass are just with with extra $mb c^2$, not individual m_n and m_p since that's how etap and etan are defined.
      etapwithm=etap+mb*clight**2
      etanwithm=etan+mb*clight**2
      etabwithm=0.5*(etapwithm+etanwithm)
c     etae already has rest-mass
      etaewithm=etae
      
c     Note that nE_nue, etc. are per kbtk, so argument of exponential is dimensionless
      Ebaryonfactor_nue = exp(-etabwithm/nE_nue)
      Eelectronfactor_nue = exp(-etaewithm/nE_nue)

      Ebaryonfactor_nuebar = exp(-etabwithm/nE_nuebar)
      Eelectronfactor_nuebar = exp(-etaewithm/nE_nuebar)

      Ebaryonfactor_mutau = exp(-etabwithm/nE_numutau)
      Eelectronfactor_mutau = exp(-etaewithm/nE_numutau)

      ntautin_nue=ntaua_nue         + Eelectronfactor_nue*ntaus_enue       + Ebaryonfactor_nue*ntausN_nue
      ntautin_nuebar=ntaua_nuebar   + Eelectronfactor_nuebar*ntaus_enuebar + Ebaryonfactor_nuebar*ntausN_nuebar
      ntautinel = (nminus_nue*ntautin_nue + nminus_nuebar*ntautin_nuebar)/(nminus_nue+nminus_nuebar+SMALL) ! for diagnostics only
      ntautinmu=ntauamu             + Eelectronfactor_mutau*ntausemu       + Ebaryonfactor_mutau*ntausNnumutau
      ntautintau=ntauatau           + Eelectronfactor_mutau*ntausetau      + Ebaryonfactor_mutau*ntausNnumutau



cccccccccccccccccccccccccccccc
c
c     Compute final neutrino emission rates
c     Volume rates : note mistake in Kaz (2005) paper equation 58
c     Should be divided by H for volume rate or keep as surface rate
c
cccccccccccccccccccccccccccccc

c      Qmel=(1./hcm)*6.62d39*T11**4/(0.5d0*tautel+0.5774+1.d0/(3.d0*tauael))
c      Qmmu=(1./hcm)*6.62d39*T11**4/(0.5d0*tautmu+0.5774+1.d0/(3.d0*tauamu))
c      Qmtau=(1./hcm)*6.62d39*T11**4/(0.5d0*tauttau+0.5774+1.d0/(3.d0*tauatau))

c      Qm=Qmel+Qmmu+Qmtau


      Qm_nue   =rate_2stream(clight,hcm,u_nue0,qtaut_nue,qtaua_nue)
      Qm_nuebar=rate_2stream(clight,hcm,u_nuebar0,qtaut_nuebar,qtaua_nuebar)
      Qmel = Qm_nue + Qm_nuebar

c      write(*,*) 'qminus',qminus_nue,Qm_nue,qtaut_nue,qtaua_nue
c      write(*,*) 'qminus2',u_nue0,etanu,hcm

      Qmmu=rate_2stream(clight,hcm,u_numu0,qtautmu,qtauamu)
      Qmtau=rate_2stream(clight,hcm,u_nutau0,qtauttau,qtauatau)

c     does NO LONGER include photons!  But never use in any other way except as total cooling
c     e.g. might want to compute Qm/Nm for energy in HARM
      Qm=Qmel+Qmmu+Qmtau


c     number rates (in 1 direction)
      Nm_nue   =rate_2stream(clight,hcm,n_nue0,ntaut_nue,ntaua_nue)
      Nm_nuebar=rate_2stream(clight,hcm,n_nuebar0,ntaut_nuebar,ntaua_nuebar)
      Nmel = Nm_nue + Nm_nuebar

      Nmmu=rate_2stream(clight,hcm,n_numu0,ntautmu,ntauamu)
      Nmtau=rate_2stream(clight,hcm,n_nutau0,ntauttau,ntauatau)

      Nm=Nmel+Nmmu+Nmtau


c      if((Qm)/Nm/mev2erg .gt.10.0) then
c         write(*,*) 'Qm',Qm,Qphoton,Nm
c         write(*,*) 'eta',etanu,etae,qE_nue*kbtk/mev2erg
c      end if

cccccccccccccccccccccccccccccc
c
c     Compute dY_e/dt
c
c     dyedt = gamman2p - (gamman2p + gammap2n)*yetot
c     Nm = total number/sec/cc
c     Nm*H = total number/sec per unit area (for each face)
c     Nm/n = number/sec/particle
cccccccccccccccccccccccccccccc


c     Below for diagnostics
      gammapeglobal = nminuspe/nptotal
      gammaAeglobal = nminusAe/nptotal
      gammapnuglobal = npluspnu/nptotal
      gammapenuglobal = npluspenu/nptotal
      
      gammaneglobal=nminusne/nntotal
      gammanglobal=nminusn/nntotal
      gammannuglobal=nplusnnu/nntotal

c      write(*,*) 'WTF',nminuspe,nptotal,nntotal

c     Old way to get dY_e/dt
      if(1.eq.1) then
         taulocal_nue = ddim((1.0-dexp(-(ntaut_nue))),0.0d0)
         taulocal_nuebar = ddim((1.0-dexp(-(ntaut_nuebar))),0.0d0)
      else
c     Assume only used for optically thin calculation
         taulocal_nue = 0.0
         taulocal_nuebar = 0.0
      end if
      gammap2nthin = (gammapeglobal+gammaAeglobal
     1     + taulocal_nuebar*(gammapnuglobal+gammapenuglobal))
      gamman2pthin = (gammaneglobal+gammanglobal + taulocal_nue*(gammannuglobal))

c     Per volume change in Y_e for optically thin case
      if(tauiter.eq.0 .OR. tauiter.eq.1) then
         dyedtthin = gamman2pthin - (gamman2pthin + gammap2nthin)*yetot
      end if

c     Compute true dyedt with optical depth and neutrino chemical potential effects
c     Use generalized treatment of optical depth
      gamman2pglobal = Nm_nuebar/nntotal
      gammap2nglobal = Nm_nue/nptotal

c      DEBUG:
c      write(*,*) 'gamma',gamman2pglobal,gammap2nglobal,Nm_nuebar,nntotal,Nm_nue,nptotal
c     Per volume change in Y_e for optically thick case
      dyedt = gamman2pglobal - (gamman2pglobal + gammap2nglobal)*yetot

c     Per volume (of size H^3) change:
c     Per area is obtained as mb*(Nm_nuebar-Nm_nue)*hcm
c     and then this is flux per face
      graddotrhouye = mb*(Nm_nuebar - Nm_nue)

      if(tauiter.eq.0 .OR. tauiter.eq.1) then
c     optical depth effects included, but \eta_\nu = 0 so neutrinos are non-thermalized
         graddotrhouyenonthermal = graddotrhouye
      end if

c      Below does not assume neutrinos are perfectly thermalized since
c     only case if \tau_{\nu}>>1
      thermalnpratio = gammap2nglobal/(dabs(gamman2pglobal)+SMALL)
c     This is some kind of global thermal value, not an absolute thermal value
c     that would be obtained by computing \eta_e assuming the thermal Y_e
      thermalye = 1.0/(1.0+thermalnpratio)





cccccccccccccccccccccccccccccc
c
c     Compute final neutrino energy and number densities
c
cccccccccccccccccccccccccccccc

c     energy densities
      u_nue =density_2stream(u_nue0,qtaut_nue,qtaua_nue)
      u_nuebar =density_2stream(u_nuebar0,qtaut_nuebar,qtaua_nuebar)
      u_nuel = u_nue+u_nuebar ! diag only
      u_numu =density_2stream(u_numu0,qtautmu,qtauamu)
      u_nutau =density_2stream(u_nutau0,qtauttau,qtauatau)

c     number densities
      n_nue =density_2stream(n_nue0,ntaut_nue,ntaua_nue)
      n_nuebar =density_2stream(n_nuebar0,ntaut_nuebar,ntaua_nuebar)
      n_nuel = n_nue+n_nuebar ! diag only
      n_numu =density_2stream(n_numu0,ntautmu,ntauamu)
      n_nutau =density_2stream(n_nutau0,ntauttau,ntauatau)

c     total internal energy density
      rho_nu = u_nuel+u_numu+u_nutau
c     total pressure
      p_nu = rho_nu/3.d0
c     entropy density
c     Often plotted is entropy per unit k_b
      s_nu=4.d0/3.d0*rho_nu/mev4toecc/tmev*mev3ergKcc/kerg !in 1/cc
c     Electron neutrino fraction (other species cancel since their chemical potential is 0)
c     NOTE THAT THE BELOW DEPENDS UPON H, so Y_\nu depends on H.
c      Ynulocal = dmax1((n_nue-n_nuebar)/nb,SMALL)
      Ynulocal = (n_nue-n_nuebar)/nb ! can be negative or positive due to optical depth supression

c     DEBUG:
c      write(*,*) 'ns',eta_nue,eta_nuebar,n_nue,n_nuebar,nb,n_nue0,n_nuebar0


c     rho_nu,p_nu,s_nu assume thermalized neutrinos so is consistent with full thermalization choice for \eta_\nu.  So don't have separate "non-thermal" values.
c     GODMARK: Not sure.    For now assume only Y_e evolution affected by thermalization degree.

ccccccccccccccccccccccccccccccccccccccccc
c
c     Compute scattering+absorption diffusion time of dominant neutrino species
c
ccccccccccccccccccccccccccccccccccccccccc

c     Kaz KM02 equation 58
c     Nm: number/s/cc, where Nm*H = number/s/area is the true quantity.
c     Tdiff = (H/c)**2*Nm

c     Diffusion time-scale for dominant species
c     Kotake et al. (2003) equation 18, where the essence of equation 17 is accounted for
c     See also Tubbs & Schramm (1975)
c     See also Shapiro & Teukolsky page 528 equations 18.5.12 -18.5.14
c     See Ruffert et al. (2006) equation A20-A21 where Tdiff = 3H/c * \tau = 3H^2/(c \lambda) so comparable
c     Limit diffusion time to no smaller than based upon speed of light across H (Kotake equation 19)

c     Kaz KM02 equation 58
c     Nm: number/s/cc, where Nm*H = number/s/area is the true quantity.
c     Tdiff = (H/c)**2*Nm


c     Number terms

c     Mean free path
      nlambda_nue = H/ntaut_nue
      nlambda_nuebar = H/ntaut_nuebar
      nlambda_numu = H/ntautmu
      nlambda_nutau = H/ntauttau

c     Diffusion time-scale
      nTdiff_nue = tdifffromlambda(clight,H,nlambda_nue)
      nTdiff_nuebar = tdifffromlambda(clight,H,nlambda_nuebar)
      nTdiff_numu = tdifffromlambda(clight,H,nlambda_numu)
      nTdiff_nutau = tdifffromlambda(clight,H,nlambda_nutau)

c     Number-rate weighted mean free path (i.e. most important emission dominates estimate of \lambda)
c     Compute scattering+absorption diffusion time of dominant neutrino species
      inlambdatot = (
     1     +nminus_nue/nlambda_nue
     1     +nminus_nuebar/nlambda_nuebar
     1     +nminusmu/nlambda_numu
     1     +nminustau/nlambda_nutau
     1     )/(nminus_nue + nminus_nuebar + nminusmu + nminustau+SMALL)

c     DEBUG:
c      write(*,*) 'NH',Elocal,rho10,T11,etae,etap,etan,xnuc,H,yeheav,kazaheav,nHnondeg
c      write(*,*) 'ntausN',ntausNP_nue,ntausNH_nue,ntausNALPHA_nue
c      write(*,*) 'taus',ntausN_nue,ntaus_enue
c      write(*,*) 'ntau',ntaua_nue,ntaustot_nue
c      write(*,*) 'lambdatot',inlambdatot,nminus_nue,nlambda_nue
c     1     ,nminus_nuebar,nlambda_nuebar
c     1     ,nminusmu,nlambda_numu
c     1     ,nminustau,nlambda_nutau

      



      nlambdatot = 1.0/inlambdatot
      nTdifftot = dmax1(3.0*H**2/(1.0*clight*nlambdatot),H/clight)

c     Below used in HARM
      lambdatot = nlambdatot
      Tdifftot = nTdifftot

c     DEBUG:
c      write(*,*) 'ntaut',ntaua_nue,ntaustot_nue,ntaut_nue,lambdatot



c     DEBUG:
c      write(*,*) 'lambdatotcheck',lambdatot,nminus_nue,nlambda_nue
c     1     ,nminus_nuebar,nlambda_nuebar
c     1     ,nminusmu,nlambda_numu
c     1     ,nminustau,nlambda_nutau


c     Energy terms

c     Mean free path
      qlambda_nue = H/qtaut_nue
      qlambda_nuebar = H/qtaut_nuebar
      qlambda_numu = H/qtautmu
      qlambda_nutau = H/qtauttau

c     Diffusion time-scale
      qTdiff_nue = tdifffromlambda(clight,H,qlambda_nue)
      qTdiff_nuebar = tdifffromlambda(clight,H,qlambda_nuebar)
      qTdiff_numu = tdifffromlambda(clight,H,qlambda_numu)
      qTdiff_nutau = tdifffromlambda(clight,H,qlambda_nutau)

      iqlambdatot = (
     1     +qminus_nue/qlambda_nue
     1     +qminus_nuebar/qlambda_nuebar
     1     +qminusmu/qlambda_numu
     1     +qminustau/qlambda_nutau
     1     )/(qminus_nue + qminus_nuebar + qminusmu + qminustau+SMALL)

      qlambdatot = 1.0/iqlambdatot
      qTdifftot = dmax1(3.0*H**2/(1.0*clight*qlambdatot),H/clight)

c     DEBUG:
c      write(*,*) 'lambda',lambdatot,Tdifftot/H*clight

c     Diffusion velocity
c      Vdifftot = H/Tdifftot


ccccccccccccccccccccccccccccccccccccccccc
c
c     Compute neutrino loss time
c
c     See Ruffert et al. (1996) equation B20-B21 and B22-B23
c     Weighted version so dominant process controls loss time
c     Ruffert use this instead of optical depth suppression method
c
ccccccccccccccccccccccccccccccccccccccccc

c     number terms
      nTloss_nue = n_nue0/nminus_nue
      nTloss_nuebar = n_nuebar0/nminus_nuebar
      nTloss_numu = n_numu0/nminusmu
      nTloss_nutau = n_nutau0/nminustau

      nTlosstot = (
     1     +n_nue0
     1     +n_nuebar0
     1     +n_numu0
     1     +n_nutau0
     1     )/(nminus_nue + nminus_nuebar + nminusmu + nminustau+SMALL)


c     energy terms
      qTloss_nue = u_nue0/qminus_nue
      qTloss_nuebar = u_nuebar0/qminus_nuebar
      qTloss_numu = u_numu0/qminusmu
      qTloss_nutau = u_nutau0/qminustau

      qTlosstot = (
     1     +u_nue0
     1     +u_nuebar0
     1     +u_numu0
     1     +u_nutau0
     1     )/(qminus_nue + qminus_nuebar + qminusmu + qminustau+SMALL)




ccccccccccccccccccccccccccccccccccccccccc
c
c     Compute energy redistribution-type diffusion time of dominant neutrino species
c
ccccccccccccccccccccccccccccccccccccccccc

c     Mean free path between inelastic collisions
      lambdain_nue = H/ntautin_nue
      lambdain_nuebar = H/ntautin_nuebar
      lambdainmu = H/ntautinmu
      lambdaintau = H/ntautintau

c     Number-rate weighted mean free path (i.e. most important emission dominates estimate of \lambda)
c      ilambdaintot = (
c     1     +nminus_nue/lambdain_nue
c     1     +nminus_nuebar/lambdain_nuebar
c     1     +nminusmu/lambdainmu
c     1     +nminustau/lambdaintau
c     1     )/(nminus_nue + nminus_nuebar + nminusmu + nminustau+SMALL)

c     JCM: This is use to evolve Y_\nu, that only involves thermalizing electron types of neutrinos
c     JCM: So for now just ignore \mu and \tau types
c     Number-rate weighted mean free path (i.e. most important emission dominates estimate of \lambda)
      ilambdaintot = 1D-150 + (
     1     +nminus_nue/lambdain_nue
     1     +nminus_nuebar/lambdain_nuebar
     1     )/(nminus_nue + nminus_nuebar + SMALL)

      lambdaintot = 1.0/ilambdaintot

      Tthermaltot = 3.0*H**2/(1.0*clight*lambdaintot)
c     Limit redistribution timescale since can't stream faster than light so trapped if motion is relativistic
      Tthermaltot = dmax1(Tthermaltot,H/clight)

c      write(*,*) 'Thermal',Tdifftot,Tthermaltot




ccccccccccccccccccccccccccccccccccccccccc
c
c     Compute Ruffert version of effective emission rates
c
ccccccccccccccccccccccccccccccccccccccccc


c     number versions
      RufNm_nue = nminus_nue/(1.0 + nTdiff_nue/nTloss_nue)
      RufNm_nuebar = nminus_nuebar/(1.0 + nTdiff_nuebar/nTloss_nuebar)
      RufNm_numu = nminusmu/(1.0 + nTdiff_numu/nTloss_numu)
      RufNm_nutau = nminustau/(1.0 + nTdiff_nutau/nTloss_nutau)

c     energy versions
      RufQm_nue = qminus_nue/(1.0 + qTdiff_nue/qTloss_nue)
      RufQm_nuebar = qminus_nuebar/(1.0 + qTdiff_nuebar/qTloss_nuebar)
      RufQm_numu = qminusmu/(1.0 + qTdiff_numu/qTloss_numu)
      RufQm_nutau = qminustau/(1.0 + qTdiff_nutau/qTloss_nutau)


c     According to Ruffert et al. (1996), final rates are:
      RufNm = RufNm_nue+RufNm_nuebar+RufNm_numu+RufNm_nutau
c     does NOT include photons!
      RufQm = RufQm_nue+RufQm_nuebar+RufQm_numu+RufQm_nutau
      Rufgraddotrhouye = mb*(RufNm_nuebar - RufNm_nue)



ccccccccccccccccccccccccccccccccccccccccc
c
c     Store energy of neutrinos
c
c
ccccccccccccccccccccccccccccccccccccccccc

      if(1.eq.0) then
c     These energies are as generated within fluid
c     however, this isn't necessarily the energy of escaping neutrinos
c     Number-rate weighted energy
         Enuglobal = kbtk*(
     1        +nminus_nue*qE_nue
     1        +nminus_nuebar*qE_nuebar
     1        +nminusmu*qE_numutau
     1        +nminustau*qE_numutau
     1        )/(nminus_nue + nminus_nuebar + nminusmu + nminustau+SMALL)

         Enueglobal=kbtk*qE_nue
         Enuebarglobal=kbtk*qE_nuebar

c     if(Enuglobal/mev2erg.gt.10.0d0) then
c     write(*,*) 'Echeck',etanu,Enueglobal/mev2erg,Enuebarglobal/mev2erg
c     end if

      else

c     These are energies of escaping neutrinos
c     These energies are used for neutrino annihilation heating rates
      Enuglobal = dmax1(Qm/Nm,0.0d0)
      Enueglobal = dmax1(Qm_nue/Nm_nue,0.0d0)
      Enuebarglobal = dmax1(Qm_nuebar/Nm_nuebar,0.0d0)
      end if



cccccccccccccccccccccccccccccccccccccc
c
c     For any ynumethod, set Ynu but changes for whichhcmmethod
c
cccccccccccccccccccccccccccccccccccccc
      
c     Below Ynu0 assignment is always true
      Ynu0 = Ynu0local

c     Below assignment of Ynu depends upon whether Ynu is to be computed or is really Ynu0 that is independent variable for table
      if(whichhcmmethod.eq.0) then
c     Then scale-free with H so using Ynu0 instead of Ynu and will have to iteratively determine Ynu0/Ynu/Yl/Ye
         Ynuiter = Ynu0local
c         write(*,*) 'Ynu0set',Ynuiter,Ynu0local
      else
         Ynuiter = Ynulocal
c         write(*,*) 'Ynuset',Ynuiter,Ynulocal
      end if






c     quasi-thermalization sets Ynu directly and keeps it, iterating on \eta_\nu until consistent solution

cccccccccccccccccccccccccccc
c
c     Compute the neutrino chemical potentials
c
c     See KM07 and km07_opticaldepths.nb
c
c     Assume using these in iterative process
c
c
cccccccccccccccccccccccccccc

c     Equilibrium values when optically thick
c     KM07 do not have -Q/T term, but Ruffert et al. (1996) do.
c     LSEOS and Shen eos provide \eta_n and \eta_p offset from a single mass m_b c^2, so no extra term for us
      eta_nue_eq = etae + etap - etan
      eta_nuebar_eq = -eta_nue_eq

c     Final chemical potentials for electrons
c     JCM: Note, using ntautin since this eta_nue/nuebar are only thermalized by inelastic processes
c      eta_nue=eta_nue_eq*(1.0-dexp(-ntautin_nue))
c      eta_nuebar=eta_nuebar_eq*(1.0-dexp(-ntautin_nuebar))


      VoCthermal_nue = 1.0/(3.0*ntautin_nue)
      VoCthermal_nue = dmin1(VoCthermal_nue,1.0d0)
      eta_nue_factor = (1.0-dexp(-ntautin_nue))*(1.0-Vocthermal_nue)

      VoCthermal_nuebar = 1.0/(3.0*ntautin_nuebar)
      VoCthermal_nuebar = dmin1(VoCthermal_nuebar,1.0d0)
      eta_nuebar_factor = (1.0-dexp(-ntautin_nuebar))*(1.0-Vocthermal_nuebar)


      eta_nue=eta_nue_eq*eta_nue_factor
      eta_nuebar=eta_nuebar_eq*eta_nuebar_factor

c     mu and tau are assumed to have 0 chemical potentials

c     Below used to pass to outside function whether convergent
      if(ynumethod.eq.0) then
         etanu=eta_nue
      end if
      if(ynumethod.eq.3) then
         etanu=eta_nue
c     Enforce behavior to be like ynumetod.eq.1
         eta_nuebar = - eta_nue
      end if

c     DEBUG:
c      write(*,*) 'ntaut=',ntaut_nue,ntaut_nuebar

c     DEBUG:
c      write(*,*) 'eta_nue,eta_nuebar=',tauiter,eta_nue,eta_nuebar




ccccccccccccccccccccccccccccccccccccccccc
c
c     Assign H-independent ratios (H-independent strictly only when fixing Y_e for each rho,T rather than solving for Y_e, and setting chemical potential of neutrinos explicitly rather than using optical depths)
c
ccccccccccccccccccccccccccccccccccccccccc

c     need for output for Matlab conversion for HARM to compute things
c
c     u_photon0, tauphotonohcm, tauphotonabsohcm
c     all the below qtaut_i qtauta_i ntaut_i ntaua_i per unit hcm
c     u_nue0, u_nuebar0, u_numu0, u_nutau0
c     lambdaintot (to obtain Tthermaltot)
c     Because mu and tau treated same, only need to store one version
c     Note that u_photon0 = u_numu0=u_nutau0, so only need to store one of them
c
c
c     In the end want to store:
c
c     u_nue0,u_nuebar0,u_numu0
c     qtaut_nueohcm, qtaut_nuebarohcm, qtautmuohcm
c     qtaua_nueohcm, qtaua_nuebarohcm, qtauamuohcm
c     n_nue0,n_nuebar0,n_numu0
c     ntaut_nueohcm, ntaut_nuebarohcm, ntautmuohcm
c     ntaua_nueohcm, ntaua_nuebarohcm, ntauamuohcm
c     lambdatot,lambdaintot
c     tauphotonohcm, tauphotonabsohcm
c
c     In matlab these all have derivatives taken vs. [rho0,utotdiff,ptotdiff,chidiff,etc.] for HARM

      unumu0=u_numu0
      unue0=u_nue0
      unuebar0=u_nuebar0

      nnumu0=n_numu0
      nnue0=n_nue0
      nnuebar0=n_nuebar0

      qtautnueohcm   =   qtaut_nue/hcm
      qtautnuebarohcm   =   qtaut_nuebar/hcm
      qtautelohcm = qtautel/hcm ! diag only

      qtauanueohcm   =   qtaua_nue/hcm
      qtauanuebarohcm   =   qtaua_nuebar/hcm
      qtauaelohcm = qtauael/hcm ! diag only

      qtautmuohcm   =   qtautmu/hcm
      qtauamuohcm   =   qtauamu/hcm
      qtauttauohcm  =   qtauttau/hcm
      qtauatauohcm  =   qtauatau/hcm


      ntautnueohcm   =   ntaut_nue/hcm
      ntautnuebarohcm   =   ntaut_nuebar/hcm
      ntautelohcm = ntautel/hcm ! diag only

      ntauanueohcm   =   ntaua_nue/hcm
      ntauanuebarohcm   =   ntaua_nuebar/hcm
      ntauaelohcm = ntauael/hcm ! diag only

      ntautmuohcm   =   ntautmu/hcm
      ntauamuohcm   =   ntauamu/hcm
      ntauttauohcm  =   ntauttau/hcm
      ntauatauohcm  =   ntauatau/hcm


cccccccccccccccc
c     just for diagnostics
      qtausohcm     =   qtaus/hcm
      ntausohcm     =   ntaus/hcm



c     just always call to enforce check that code is correct
c     Outputs not used (at all!) if whichdatatype==4 (i.e. whichhcmmethod==0)
      call computefinal_fromhcm(
     1     clight,mb
     1     ,rhob,kbtk,hcm
     1     ,unue0,unuebar0,unumu0
     1     ,qtautnueohcm, qtautnuebarohcm, qtautmuohcm
     1     ,qtauanueohcm, qtauanuebarohcm, qtauamuohcm
     1     ,nnue0,nnuebar0,nnumu0
     1     ,ntautnueohcm, ntautnuebarohcm, ntautmuohcm
     1     ,ntauanueohcm, ntauanuebarohcm, ntauamuohcm
     1     ,lambdatot,lambdaintot
     1     ,tauphotonohcm, tauphotonabsohcm
     1     ,nnueth0,nnuebarth0
     1     ,Qphoton,Qm,graddotrhouye,Tthermaltot,Tdifftot,rho_nu,p_nu,s_nu,Ynulocal,Ynuthermal ! outputs
     1     ,Enuglobal,Enueglobal,Enuebarglobal ! more outputs
     1     )


c     Assign outputted Ynu from 2-stream approximation.  Should be same result as from prior 2-stream calls.  Not sure why do both -- why not just use function if want to test it?
      Ynu = Ynulocal


c     DEBUG:
c      write(*,*) 'lambdatotcheck2',lambdatot,nminus_nue,nlambda_nue
c     1     ,nminus_nuebar,nlambda_nuebar
c     1     ,nminusmu,nlambda_numu
c     1     ,nminustau,nlambda_nutau




c     DEBUG: (compare to harm extras/processed)
c     Both should be in the same order as in HARM table/output
c      write(*,*) 'DEBUGFORHARM',tauiter
c      write(*,*) 'extras'
c     1	,qtautnueohcm
c     1	,qtauanueohcm
c     1	,qtautnuebarohcm
c     1	,qtauanuebarohcm
c     1	,qtautmuohcm
c     1	,qtauamuohcm
c     1	,ntautnueohcm
c     1	,ntauanueohcm
c     1	,ntautnuebarohcm
c     1	,ntauanuebarohcm
c     1	,ntautmuohcm
c     1	,ntauamuohcm
c     1	,unue0
c     1	,unuebar0
c     1	,unumu0
c     1	,nnue0
c     1	,nnuebar0
c     1	,nnumu0
c     1	,lambdatot
c     1	,lambdaintot
c     1	,tauphotonohcm
c     1	,tauphotonabsohcm
c     1	,nnueth0
c     1	,nnuebarth0
c
c      write(*,*) 'processed'
c     1     ,Qphoton,Qm,graddotrhouye,Tthermaltot,Tdifftot,rho_nu,p_nu,s_nu,Ynulocal,Ynuthermal
c     1     ,Enuglobal,Enueglobal,Enuebarglobal


      return
      end



C=====================================================================
      subroutine compute_old_Ynp_Ypn(yp,etan,etap,Ynp,Ypn)

c     See Bruenn (1985) equation C14 and Ruffert et al. (1996) A13-A14
c     Assume Ye=Yp and Yn=1-Yp
C=====================================================================
      implicit none
      
C=====================================================================
      real*8 yp,etan,etap
      real*8 Ynp,Ypn
C=====================================================================

      include 'const.dek'

      Ynp = (2.0d0*yp-1.0d0)/(dexp(etap-etan)-1.0d0)
      Ypn = dexp(etap-etan)*Ynp

c      write(*,*) 'raw',yp,etap,etan
c      write(*,*) 'rawYnp',Ynp,Ypn

c     When etap~etan, can have problems.

c     Also, when nucleons not dominate, then using above with yp=Y_{p,free} and using full etap and etap (not just free!) can lead to issues
c     For example, Y_{p,f}>0.5 can occur, and then issues appear with both Ynp and Ypn being negative
c     Assume if negative, then free nucleons just unimportant

c     Now fix-up any issues with these calculations
      if(Ynp.gt.1.0d0) then
         Ynp = 1.0d0
      end if
      if(Ynp.lt.0.0d0) then
         Ynp = 0.0d0
      end if

      if(Ypn.gt.1.0d0) then
         Ypn = 1.0d0
      end if
      if(Ypn.lt.0.0d0) then
         Ypn = 0.0d0
      end if
c     Also enforce total number of baryons by diminishing Ypn
      if(Ypn+Ynp.gt.1.0d0) then
         Ypn = 1.0d0-Ynp
      end if

      return
      end

C=====================================================================
      subroutine compute_new_Ynp_Ypn(yp,etan,etap,Ynp,Ypn)

c     See Bruenn (1985) equation C14 and Ruffert et al. (1996) A13-A14
c     Assume Ye=Yp and Yn=1-Yp
c     Ynp is per unit n_n and Ypn is per unit n_p
C=====================================================================
      implicit none
      
C=====================================================================
      real*8 yp,etan,etap
      real*8 Ynp,Ypn
c     locals
      real*8 myyp
C=====================================================================

      include 'const.dek'

c     restrict to avoid NaN/Inf
      myyp=dmax1(yp,1.0D-3)

      Ynp = (2.0d0*myyp-1.0d0)/(1.0d0-myyp)/(dexp(etap-etan)-1.0d0)
      Ypn = (1.0d0-2.0d0*myyp)/(myyp)/(dexp(etan-etap)-1.0d0)

c      write(*,*) 'raw',yp,etap,etan
c      write(*,*) 'rawYnp',Ynp,Ypn

c     When etap~etan, can have accuracy problems since nuclear EOS can become inaccurate

c     Also, when nucleons not dominate, then using above with yp=Y_{p,free} and using full etap and etap (not just free!) can lead to issues
c     For example, Y_{p,f}>0.5 can occur, and then issues appear with both Ynp and Ypn being negative
c     Assume if negative, then free nucleons just unimportant

c     Below fix-ups make much more sense now using new Ynp and Ypn
c     Now fix-up any issues with these calculations 
      if(Ynp.gt.1.0d0) then
         Ynp = 1.0d0
      end if
      if(Ynp.lt.0.0d0) then
         Ynp = 0.0d0
      end if

      if(Ypn.gt.1.0d0) then
         Ypn = 1.0d0
      end if
      if(Ypn.lt.0.0d0) then
         Ypn = 0.0d0
      end if

      return
      end

      

      include 'tau_neededbyharm.f'







C=====================================================================
c
c Here begins the neutrino opacities
c
C=====================================================================








C=====================================================================
      real*8 function qdotne(T11,etae,etap,etan,eta_nuebar,extrablock,nnfree)


c     Cooling by neutronization, or URCA, reactions.  Calculate the
c     asymptotic nondegenerate and degenerate rates and interpolate for
c     the current value of deg.  This cooling applies only to electron
c     neutrinos.

c      implicit double precision (a-h,o-z)
      implicit none


C=====================================================================
      real*8 T11,etae,etap,etan,eta_nuebar,extrablock,nnfree
C=====================================================================
      real*8 qne2
C=====================================================================
c      include 'kazeos1.dek'
      include 'const.dek'

cc      qdotNenondeg=9.2d33*T11**6*rho10
c      qdotNenondeg=9.2d33*T11**6*rho10 * xnuc  !check1
c      qdotNedeg=1.1d31*etae**9*T11**9
c      qdotNe=qdotNenondeg/(1.+deg*deg)+deg*deg*qdotNedeg/(1.d0+deg*deg)


      qdotne=qne2(T11,etae,eta_nuebar,extrablock)*nnfree


      return
      end




C=====================================================================
      real*8 function qdotn(T11,etae,etap,etan,eta_nuebar,extrablock,nnfree)


c     Cooling by deneutronization, or URCA, reactions.  Calculate the
c     asymptotic nondegenerate and degenerate rates and interpolate for
c     the current value of deg.  This cooling applies only to electron
c     neutrinos.

c      implicit double precision (a-h,o-z)
      implicit none


C=====================================================================
      real*8 T11,etae,etap,etan,eta_nuebar,extrablock,nnfree
C=====================================================================
      real*8 qn2
C=====================================================================
c      include 'kazeos1.dek'
      include 'const.dek'

 
      qdotn=qn2(T11,etae,eta_nuebar,extrablock)*nnfree


      return
      end


C=====================================================================
      real*8 function qdotpe(T11,etae,etap,etan,eta_nue,extrablock,npfree)


c     Cooling by neutronization, or URCA, reactions.  Calculate the
c     asymptotic nondegenerate and degenerate rates and interpolate for
c     the current value of deg.  This cooling applies only to electron
c     neutrinos.

c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,etae,etap,etan,eta_nue,extrablock,npfree
C=====================================================================
      real*8 qpe2
C=====================================================================
c      include 'kazeos1.dek'
      include 'const.dek'

      qdotpe= qpe2(T11,etae,eta_nue,extrablock)*npfree


      return
      end





C=====================================================================
      real*8 function qdotAe(T11,etae,etap,etan,eta_nue
     1     ,abar,zbar,extrablock,nA)


c      implicit double precision (a-h,o-z)
      implicit none

c      common /check2/xnuc,xnucprev,dr

C=====================================================================
      real*8 T11,etae,etap,etan,eta_nue
     1     ,abar,zbar,extrablock,nA
C=====================================================================
      real*8 qAe2
C=====================================================================
c      include 'kazeos1.dek'
      include 'const.dek'

c     JCM: below takes into account fact that qAe2() returns rate per nuclei with atomic mass A
c     And Qfinal = qAe2 n_A = qAe2 (# density of bound baryons)/ABAR
      qdotAe = qAe2(T11,etae,etap,etan,eta_nue,abar,zbar,extrablock)*nA

      return
      end



C=====================================================================
      real*8 function ndotpair_ele(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11
      real*8 fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
C=====================================================================
      include 'const.dek'

c     See Ruffert et al. (1996) equation B8
c     See ruffert_emissionrates.nb
c     Including BOTH anti/normal particles
      
      ndotpair_ele = 1.3913D36*T11**8*fd3_epos*fd3_eneg


      return
      end


C=====================================================================
      real*8 function qdotpair_ele(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11
      real*8 fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
C=====================================================================
      real*8 ndotpair_ele ! function
      real*8 Efactor

      include 'const.dek'

c     Including BOTH anti/normal particles

      Efactor = 6.9035D-6*T11*(fd4_eneg/fd3_eneg + fd4_epos/fd3_epos)
     
      qdotpair_ele = Efactor*ndotpair_ele(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )


      return
      end


C=====================================================================
      real*8 function ndotpair_mutau(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11
      real*8 fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
C=====================================================================
      real*8 ndotpair_ele ! function

      include 'const.dek'

c     For EACH \mu, \tau (including BOTH anti/normal particles)


      ndotpair_mutau = 0.214749*ndotpair_ele(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )


      return
      end


C=====================================================================
      real*8 function qdotpair_mutau(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11
      real*8 fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
C=====================================================================
      real*8 qdotpair_ele ! function

      include 'const.dek'

c     For EACH \mu, \tau (including BOTH anti/normal particles)


      qdotpair_mutau = 0.214749*qdotpair_ele(T11
     1     ,fd2_eneg,fd2_epos,fd3_eneg,fd3_epos
     1     ,fd4_eneg,fd4_epos,fd5_eneg,fd5_epos
     1     )


      return
      end









C=====================================================================
      real*8 function qdotpair_ele_old(T11,deg)

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,deg
C=====================================================================

c     Cooling by pair annihilation processes.  This cooling applies to
c     all three neutrino species.


      qdotpair_ele_old=3.4d33*T11**9/(1.d0+deg*deg)    !total

      return
      end


C=====================================================================
      real*8 function qdotpair_mutau_old(T11,deg)

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,deg
C=====================================================================

c     Cooling by pair annihilation processes.  This cooling applies to
c     all three neutrino species.


      qdotpair_mutau_old=0.7d33*T11**9/(1.d0+deg*deg)    !total

      return
      end

C=====================================================================
      real*8 function ndotpair_ele_old(T11,deg)

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,deg
C=====================================================================
      real*8 qdotpair_ele_old
      real*8 qdotpairlocal

      include 'const.dek'

c     Cooling by pair annihilation processes.  This cooling applies to
c     all three neutrino species.

      qdotpairlocal = qdotpair_ele_old(T11,deg)

c     JCM: Fake version
      ndotpair_ele_old=qdotpairlocal/(3.15*kerg*T11*1D11)

      return
      end


C=====================================================================
      real*8 function ndotpair_mutau_old(T11,deg)

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,deg
C=====================================================================
      real*8 qdotpair_mutau_old
      real*8 qdotpairlocal

      include 'const.dek'

c     Cooling by pair annihilation processes.  This cooling applies to
c     all three neutrino species.

      qdotpairlocal = qdotpair_mutau_old(T11,deg)

c     JCM: assumes energy per neutrino of 3.15kbT (true if \eta_\nu=0)
      ndotpair_mutau_old=qdotpairlocal/(3.15*kerg*T11*1D11)

      return
      end






C=====================================================================
      subroutine qdotbrem_both (rho10,rho10nondeg,T11,qdotbremnondegen,qdotbremdegen,Xpnondeg,Xnnondeg)
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
c     Passed/Returned
      real*8 rho10,T11
      real*8 etan,etap
c     Locals:
      real*8 qdotbremnondegen,qdotbremdegen
      real*8 rho10nondeg,rho10deg
      real*8 kbtk
      real*8 Xpnondeg,Xnnondeg
      
      include 'const.dek'


      kbtk = kerg*T11*1.0D11

c     Cooling by neucleon-nucleon bremsstrahlung.  Applies to EACH
c     neutrino species and accounts for both neutrino and antineutrinos


c     nucleons that are non-degenerate still
c     divisor causes numerator to vary from 0 to rho10
c      rho10nondeg=rho10/dmax1(dexp(etanbruenn-etapbruenn),dexp(etapbruenn-etanbruenn),1.0d0)
c     nucleons that are degenerate
      rho10deg=dmax1(rho10-rho10nondeg,0.0d0)
      
c     JCM: my interpolation of KM02 formulae
c      qdotbremnondegen=1.5d27*rho10nondeg**2.0*T11**5.5d0
c     From BRT06:
      qdotbremnondegen=2.18d27*rho10nondeg**2.0*T11**5.5d0
      qdotbremnondegen=qdotbremnondegen*(Xnnondeg**2+Xpnondeg**2+28.0/3.0*Xnnondeg*Xpnondeg)
     
      qdotbremdegen = 3.4D32*rho10deg**(1.0/3.0)*T11**8.0



      return
      end



C=====================================================================
      real*8 function qdotbrem (rho10,rho10nondeg,T11,Xpnondeg,Xnnondeg)
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
c     Passed/Returned
      real*8 rho10,rho10nondeg,T11
c     Locals:
      real*8 qdotbremnondegen,qdotbremdegen
      real*8 Xpnondeg,Xnnondeg


c      write(*,*) 'qdotbrem',rho10,rho10nondeg,T11,Xpnondeg,Xnnondeg


      call qdotbrem_both(rho10,rho10nondeg,T11,qdotbremnondegen,qdotbremdegen,Xpnondeg,Xnnondeg)

c     1/3 is because Kaz says that his brem formulae apply to a sum over all species
      qdotbrem = (1.0/3.0)*(qdotbremnondegen + qdotbremdegen)




      return
      end


C=====================================================================
      real*8 function ndotbrem (rho10,rho10nondeg,T11,etap,etan,Elocalnue,Elocalnuebar,Elocalnumutau,Xpnondeg,Xnnondeg)
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
c     Passed/Returned
      real*8 rho10,rho10nondeg,T11
      real*8 etan,etap
c     Locals:
      real*8 qdotbremnondegen,qdotbremdegen
      real*8 ndotbremnondegen,ndotbremdegen
      real*8 Enondeg,Edeg,kbtk
      real*8 etanbruenn,etapbruenn
      real*8 Xpnondeg,Xnnondeg
      real*8 Elocalnue,Elocalnuebar,Elocalnumutau

      include 'const.dek'


      call qdotbrem_both(rho10,rho10nondeg,T11,qdotbremnondegen,qdotbremdegen,Xpnondeg,Xnnondeg)


      kbtk=kerg*T11*1.0D11
c      Enondeg = (3.15*kbtk)
      Enondeg = 0.333333*(Elocalnuebar+Elocalnuebar+Elocalnumutau)
c      Edeg = kbtk*dmax1((etan-etap),(etap-etan))
      Edeg = dmax1(etan*kbtk,etap*kbtk,Enondeg)

c     Approximate number rate, where as stated above qdotbrem_both returns both \nu and \bar{\nu} and all species e,\mu,\tau
      ndotbremnondegen=(qdotbremnondegen/6.0)/Elocalnue + (qdotbremnondegen/6.0)/Elocalnuebar + 
     ! (2.0*qdotbremnondegen/3.0)/Elocalnumutau
      ndotbremdegen = qdotbremdegen/Edeg

c     JCM: my interpolation of KM02 formulae
      ndotbrem = ndotbremnondegen + ndotbremdegen


      return
      end









C=====================================================================
      real*8 function qdotplasmonele (T11,etae)
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,etae,eta_nue,eta_nuebar
      real*8 gammap

C=====================================================================
c     Cooling by plasmon decay.  Applies only to electron neutrinos.
      include 'const.dek'
      real*8 ndotplasmonele ! function
      real*8 kbtk,Efactor

c     See Ruffert et al. (1996)

      gammap=3.21295d-2*sqrt(pi**2 + 3.d0*etae*etae)
      kbtk = (kerg*T11*1.0D11)
      Efactor = 0.5*kbtk*(2.0+gammap**2/(1.0+gammap))

      qdotplasmonele = Efactor*ndotplasmonele(T11,etae)

c     Kaz code:
c      qdotplasmonele=1.52428d32*T11**9*gammap**6*exp(-gammap)
c     &     *(2.d0+2.d0*gammap+gammap*gammap)


      return
      end


C=====================================================================
      real*8 function ndotplasmonele (T11,etae)
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,etae,eta_nue,eta_nuebar
      real*8 gammap

C=====================================================================
c     Cooling by plasmon decay.  Applies only to electron neutrinos.
      include 'const.dek'

c     See Ruffert et al. (1996).  Includes both normal and anti-neutrinos of electron type


      gammap=3.21295d-2*sqrt(pi**2 + 3.d0*etae*etae)
      ndotplasmonele=4.41713D37*T11**8*gammap**6*exp(-gammap)
     &     *(1.d0 + gammap)

      return
      end



C=====================================================================
      real*8 function qdotplasmonmutau (T11,etae)
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,etae
      real*8 qdotplasmonele !function

C=====================================================================
c     Cooling by plasmon decay.  Mu and tau each have this
      include 'const.dek'

c     See Ruffert et al. (1996)

      qdotplasmonmutau = 0.00173611 * qdotplasmonele(T11,etae)


      return
      end

C=====================================================================
      real*8 function ndotplasmonmutau (T11,etae)
c      implicit double precision (a-h,o-z)
      implicit none

C=====================================================================
      real*8 T11,etae
      real*8 ndotplasmonele !function

C=====================================================================
c     Cooling by plasmon decay.  Mu and tau each have this
      include 'const.dek'

c     See Ruffert et al. (1996)

      ndotplasmonmutau = 0.00173611 * ndotplasmonele(T11,etae)


      return
      end








C=====================================================================
      real*8 function tauabs (rate,density,H)

C=====================================================================
c     Given the cooling rate for a particular neutrino species,
c     calculated the absorptive opacity by: tau_a = q- H/(7/2) sigma T^4
c     Old method combined q-=q-(\nu_e)+q-(\bar{\nu_e})

c     Now each species is treated separately as in KM07

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none
      
C=====================================================================
      real*8 rate,density,H

C=====================================================================

c      tauabs=qminus*H/(1.98d40*T11**4)
      include 'const.dek'


      tauabs=rate/(clight*density)*H


      return
      end





C=====================================================================
      real*8 function tauscattenue(Enu,rho10,T11,etae,xnuc,H,yetotlocal)

C=====================================================================
c     From KM07, electron species (normal neutrino)
C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none
      
C=====================================================================
c     Passed:
      real*8 rho10,xnuc,T11,H
      real*8 etae,yetotlocal
      real*8 Enu


C=====================================================================

      tauscattenue = 2.53503D-8 * (4.0+etae)*rho10*T11**2*yetotlocal * H * (Enu/4.10596)

c      write(*,*) 'tauscattenue=',etae,rho10,T11,yetotlocal,H,Enu,tauscattenue

      return
      end



C=====================================================================
      real*8 function tauscattenuebar(Enu,rho10,T11,etae,xnuc,H,yetotlocal)

C=====================================================================
c     From KM07, electron species (anti-neutrino)
C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none
      
C=====================================================================
c     Passed:
      real*8 rho10,xnuc,T11,H
      real*8 etae,yetotlocal
      real*8 Enu


C=====================================================================

      tauscattenuebar = 1.06153D-6 * (4.0+etae)*rho10*T11**2*yetotlocal * H * (Enu/4.10596)

      return
      end

C=====================================================================
      real*8 function tauscattenu_mutau(Enu,rho10,T11,etae,xnuc,H,yetotlocal)

C=====================================================================
c     From KM07, for EACH \mu and \tau species
C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none
      
C=====================================================================
c     Passed:
      real*8 rho10,xnuc,T11,H
      real*8 etae, yetotlocal
      real*8 Enu


C=====================================================================

      tauscattenu_mutau = 7.72359D-9 * (4.0+etae)*rho10*T11**2*yetotlocal * H * (Enu/4.10596)

      return
      end







c     See km07_scatteringopacities.nb
c     Notice how degeneracy is taken into account differently than for pair capture
c     Shapiro & Tuekolsky
C=====================================================================
      real*8 function tauscattNP(Esqnu,rho10,T11,etae,etap,etan
     1     ,xnuc,H,yefree
     1     )

C=====================================================================
c     Scattering opacity for each neutrino species.  Taken from di
c     Matteo et al. (2002).  According to Tiziana, the rate given in the
c     paper is for each neutrino species.

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none
      
C=====================================================================
c     Passed:
      real*8 rho10,xnuc,T11,H
      real*8 etae,etap,etan,yefree,yebound,abarbound,nAnondeg
      real*8 Esqnu
c     locals:
      real*8 Ynn,Ypp,tauscattneut,tauscattprot
      real*8 rhoA,rhoA10
      real*8 kbtk,etanbruenn,etapbruenn
         
      include 'const.dek'

C=====================================================================

c     Original Kaz expression
c      tauscatt=2.7d-7*rho10*T11*T11*H

c     JCM: I disagree with this calculation.  Seems they assumed
c     Yn=1 AND Yp=1.  Yp=1/2 gives me 1.43D-7 H rho10 T11^2
c      tauscattNP=1.4d-7*rho10*T11*T11*H

c     Go ahead and use full Y_e dependence
c     Below is ok, but doesn't account for nucleon degeneracy
c      tauscattNP=1.51729D-7 * (xnuc*rho10) * T11**2 * (1.0-0.112749*yetotlocal) * H

c     Assume for consistency that this degeneracy parameter is like Ynp and Ypn
      kbtk=kerg*T11*1.0D11
c      etanruffert = etan + (amu-mn)*clight**2/kbtk
c      etapruffert = etap + (amu-mp)*clight**2/kbtk
      etanbruenn = etan
      etapbruenn = etap


c     use Bruenn1985-KM07-Ruffert expression
      Ynn = (1.0d0-yefree)/(1.0d0+(2.0d0/3.0d0)*dmax1(0.0d0,etanbruenn))
      tauscattneut = 2.28552D-7 * (xnuc*rho10*Ynn)*T11**2 * (Esqnu/20.8135D0) * H

      Ypp = yefree/(1.0d0+(2.0d0/3.0d0)*dmax1(0.0d0,etapbruenn))
      tauscattprot = 2.02783D-7 * (xnuc*rho10*Ypp)*T11**2 * (Esqnu/20.8135D0) * H

c     Considered total opacity to neutrinos of one species
      tauscattNP = tauscattneut + tauscattprot




c     DEBUG:
c      write(*,*) xnuc,rho10,abar,zbar,H
c      write(*,*) 'taus=',tauscatt,tauscattNP,tauscattA

      return
      end






c     See km07_scatteringopacities.nb
c     Shapiro & Tuekolsky
c
c
c     used for *each* species of heavy nuclei (A=4, heavy)

C=====================================================================
      real*8 function tauscattA(Esqnu,rho10,T11,etae,etap,etan
     1     ,xnuc,H,yeA,abarA
     1     ,nA
     1     )

C=====================================================================
c     Scattering opacity for each neutrino species.  Taken from di
c     Matteo et al. (2002).  According to Tiziana, the rate given in the
c     paper is for each neutrino species.

C=====================================================================
c      implicit double precision (a-h,o-z)
      implicit none
      
C=====================================================================
c     Passed:
      real*8 rho10,xnuc,T11,H
      real*8 etae,etap,etan,yeA,abarA,nA
      real*8 Esqnu
c     locals:
      real*8 rhoA,rhoA10
         
      include 'const.dek'

C=====================================================================



c     Also, from Burrows, A. (1990) review, equation 2b, I get using the Di Matteo averaing of E_\nu
c     tau = 1.03D-7 H rho10 T11^2, which is much closer to seemingly correct estimate

c      tauscattA = tauscatt*3.0*(A/56.0)
c     Use directly Burrows, A. (1990) equation (2a)
c      tauscattA = 3.07841D-7*rho10*T11*T11*(abarA/56.0)*H


c     Get baryon-mass density of particles
      rhoA = nA*(mb*abarA)
      rhoA10 = rhoA/1.0D10


      if(1.eq.1) then
c     Use Shapiro page 526 equation 18.5.6 from Tubbs & Schramm (1975)
c     Note that the cross section is given per nuclei and so when multiplying by \rho need to consider per baryon cross section
c     Other parameters needed are used as in Di Matteo et al. (2002)
c     Below is consistent with Burrows estimate if Z/A\sim 0.5
c     Below is true only for E_nu << 300 A^{-1/3} MeV, so even for Iron this works up to \sim 100MeV
         tauscattA = 2.1959D-6 * (rhoA10) * T11**2 * (abarA/56.0*(1.0-1.08*yeA)**2) * Esqnu/(20.8135D0) * H

c     Protect against large Y_e
         if(tauscattA.lt.0.0) then
            tauscattA=0.0
         end if


      else

c     Note that KM07 expression has neutrino energy dependence but no Y_e dependence.
c     KM07 has typo, should premultiply by A^2
c     In the end tauscattA\propto abarA if correcting KM07
c         tauscattA = 3.89025D-8 * (rhoA10/abarA)*T11**2*(Esqnu/20.8135D0) * H
c         tauscattA = 3.89025D-8 * (nA*mb/1.0D10)*T11**2*(Esqnu/20.8135D0) * H
         tauscattA = 2.17854D-6 * (abarA/56.0) * (rhoA10)*T11**2*(Esqnu/20.8135D0) * H

      end if

c     DEBUG:
c      write(*,*) 'tauscattA',abarA,nA,rhoA10,T11,Esqnu,H,tauscattA
c      write(*,*) xnuc,rho10,abar,zbar,H
c      write(*,*) 'taus=',tauscatt,tauscattNP,tauscattA

      return
      end


