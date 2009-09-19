c..
c.. compute some auxillary EOS quantities to make other EOSs more like Kaz EOS as far as variable names and some extra calculations
c
c
c     Note that rhob,tk,hcm,tdyn are set globally into kazeos.dek
c
c
      subroutine compute_kazlike_eos(jj)
      implicit none
      save
c..declare passed variables
      integer jj
      real*8 nb,neleposi,myye,yefit2
      double precision computeye ! function

      real*8 R,Re,a,K,K2
      real*8 mue,Ebin, m_N, Ebinrat,P_Nkaz

c..delcare local/global variables
      include 'eosparms.f'
      include 'vector_eos.dek'
      include 'kazeos.dek'
      include 'kaz_state.dek'
c      include 'vector_eos.extra4kaz.dek'
c..
      include 'const.dek'


    
      
      
c     write(6,10) jj


ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Translate HELM/TIMMES/LSEOS(KAZ converted) variables into Kaz variables in kazeos.dek
c     
ccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccc
c
c     Translate into Kaz-like variables the independent variables
c
cccccccccccccc
         rhob=den_row(jj)
         tk=temp_row(jj)
         hcm=hcm_row(jj)
c     never is TDYN used, but Ye is used to set zbar and abar, which eventually is used to set this
c         tdynorye=zbar_row(jj)/abarnum_row(jj)
c     Notice small error since recompute zbar,abarnum
         tdynorye=tdynorye_row(jj)
         tdynorynu=tdynorynu_row(jj)





cccccccccccccccc
c     Translate labels from HELM/TIMMES/(converted Kaz) format to Kaz format


ccccccccccc
c     Presure:
      p_photon=prad_row(jj)
      p_eleposi=pele_row(jj)+ppos_row(jj)
      p_N=pion_row(jj)+pcou_row(jj)
c     internal energy density for nucleons (always non-relativistic and non-degenerate)
c     u_N=p_N/(5.0/3.0-1.0)
      
      
ccccccccccc
c     Internal energy:
      rho_photon = erad_row(jj)*rhob
      rho_eleposi = (eele_row(jj) + epos_row(jj))*rhob
c     rho = rest-mass + internal energy (this is for free + bound nucleons)
      u_N = (eion_row(jj)+ecou_row(jj))*rhob
      rho_N = rhob*clight**2 + u_N
c     rest-mass negligible
      
ccccccccccc
c     Entropy densities in erg/K/kb/cc=1/cc:
      s_photon = (srad_row(jj)/kerg)*rhob
      s_eleposi = ((sele_row(jj)+spos_row(jj))/kerg)*rhob
      s_N = ((sion_row(jj) + scou_row(jj))/kerg)*rhob
c     s_N = 2.50*p_N/tk/kb         !in erg/K/kb/cc=1/cc





c     Other Kaz things that are outputted but stored in different variable
c     Store array into singles
      call storeback_row(jj)






      if( (whicheleeos.eq.0).OR.(whicheleeos.eq.1).OR.(whicheleeos.eq.2)
     1     ) then
         


         if(kazlikeneutrinos.eq.0) then

cccccccccccccccccccccccccccccccccccccccc
c     Compute neutrino emission rates (not using Kaz)
            call compute_kazlikenurates(jj, tk, rhob, Qm)
c     call compute_kazlikenurates(jj,Qm)
c     Presume the Kaz-like variables have already been set to 0 in init_row() in jon_lsbox.f

         else if(kazlikeneutrinos.eq.1) then


c     Call all of Kaz stuff assuming only 
c
c     etae=etaele_row(jj)
c     call kaz_physics_etae(etaele_row(jj))


c     Call ACTUAL Kaz code for neutrino-related things
c
c     This overrides any _row assignment from above except xnuc,npratiofree,yetot,yefree,yebound that (if used) are computed
c     by the nuclear EOS
c     In this way the results from this call are actually more correct than using full Kaz code since xnuc set more accurately


c     Ensure that neutrino terms aren't there by getting rid of them if they are
            ptot_row(jj) = ptot_row(jj) - p_nu_row(jj)
            etot_row(jj) = etot_row(jj) - rho_nu_row(jj)
            stot_row(jj) = stot_row(jj) - s_nu_row(jj)

c     Convert HELM EOS format of _row -> kaz like format of variables
            call preparecall2kazeos(jj)

c     Compute neutrino terms, perhaps iterating to get solution for a given optical depth or dynamical timescale
            call kaz_physics_neutrinos_etae(etae,etap,etan,etanu)

c     Convert neutrino stuff to HELM form from KAZ form
c     erg/K/g
            s_nu_row(jj) = kerg*s_nu/den_row(jj)
c     erg/g
            rho_nu_row(jj) = rho_nu/den_row(jj)
c     same
            p_nu_row(jj) = p_nu


c     Now assign totals
            ptot_row(jj) = ptot_row(jj) + p_nu_row(jj)
            etot_row(jj) = etot_row(jj) + rho_nu_row(jj)
c     entropy here is still in "HELM" form of erg/g/K
            stot_row(jj) = stot_row(jj) + s_nu_row(jj)
     
c     NOTE: All non-HELM quantities (pure Kaz quantities) are stored in singles, so no need to have a global _row for each of these!
c     That is, we just output singles just after this function is called

c            write(*,*) 'postkazneutrino',lambdatot
            

         end if



c     GODMARK: Below calculation isn't general
         if(0.eq.1) then
cccccccccccccccccccccccccccccccccccccccc
c     Compute npratiofree
c     
c     override npratiofree from kaz estimate
            nb = (rhob/mb)
c     xne actually includes both electrons and positrons
c     xnem is electrons associated with ions such that xnem/xni=Ye
c     neleposi = (xne_row(jj)-xnp_row(jj)+xnem_row(jj))
c     assume additional electrons and positrons always cancel
            neleposi = (xnem_row(jj))
c     Y_e = n_e/(n_n+n_p) = n_e/n_b
c     myye = neleposi/nb
            myye = zbar_row(jj)/abarnum_row(jj)
c     myye=0.5
c     This npratio is total npratio
c     npratio = (1.0d0 - myye)/myye
            yetot=myye

c     Using HELM version of npratiofree
c     call computeyetot(rhob,tk,npratiofree,xnuc,yetot) ! Y_e total

c     Was only applied if xnuc=1
            yefit2=computeye(rhob)
c     Now only free npratiofree
            npratiofree = (1.0d0 - yefit2)/yefit2
         end if

      end if

      

ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Kaz-like EOS w/ Kaz-like photodisintegration correction
c     
ccccccccccccccccccccccccccccccccccccccccccccc

      if (OUTPUTTYPE.eq.1 .OR. OUTPUTTYPE.eq.2) then
c     Some constants from Kaz code
         R=kerg/mp
         Re=kerg/me
         a=5.6704E-5 * 4 / clight
         K=(2*pi*hcgs*clight/3)*(3/(8*pi*mn))**(4/3)
         K2=9.9E12
c     Assumes 1 electron per (proton + neutron)
c     GODMARK: inconsistent with abar and zbar in HELM code
         mue=2.0
         Ebin = 28.3
         m_N = 938.919
         Ebinrat = Ebin/(4.d0*m_N)
c     
      end if

      

      if (OUTPUTTYPE.eq.1) then
c     
c     p_N=Pgas(rho10,T11)*(1.d0+3.d0*xnuc)/4.d0   !in erg/cc
c     Modified by JCM using Kaz's new formula
c     Below 2 in MeV
         P_Nkaz=8.26d28*(rhob/1E10)*(tk/1E11) !in erg/cc
         dissfactor=(1.d0+(3.d0-Ebinrat)*xnuc)/(1.d0-Ebinrat)
         pkaz_N=P_Nkaz*dissfactor/4.d0 !in erg/cc
         ukaz_N=pkaz_N/(5.0/3.0-1.0)
c     p_N=Pgas(rho10,T11)*(1.d0+3.d0*xnuc)/4.d0   !in erg/cc

c     s_N= 2.50*p_N/mev4toecc/tmev*mev3ergKcc/kerg !in 1/cc
         skaz_N= 2.50*pkaz_N/tk/kerg !in 1/cc

      end if

ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     PWF99 EOS w/ Kaz-like photodisintegration correction
c     
ccccccccccccccccccccccccccccccccccccccccccccc

      if (OUTPUTTYPE.eq.2) then
c     Assumes ele-pos are degenerate
         ppwf_ele = K*(rhob/mue)**(4.0/3.0)
         ppwf_rad = 11.0/12.0*a*tk**4.0
c     Same as Kaz term above
         ppwf_N = pkaz_N

         upwf_ele = 3.0*K*(rhob/mue)**(1.0/3.0)*rhob
         upwf_rad = ppwf_rad/(4.0/3.0-1.0)
         upwf_N = ppwf_N/(5.0/3.0-1.0)

         spwf_N = 2.5*pkaz_N/tk
c     Use HELM for entropy of electrons and radiation
         spwf_ele = s_eleposi
         spwf_rad = s_photon

      end if











ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     FINAL PRESSURE, INTERNAL ENERGY,ENTROPY, AND NEUTRINO EMISSION RATES
c     
ccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     HELM/TIMMES (or converted true Kaz)
c     
ccccccccccccccccccccccccccccccccccccccccccccc
      if (OUTPUTTYPE.eq.0) then
c     
c     HELM total pressure
         p_tot = ptot_row(jj)
         u_tot = etot_row(jj)*den_row(jj)
c     entropy density in erg/cc/K/kb = 1/cc
         s_tot = (stot_row(jj)/kerg)*den_row(jj)

c     Convert neutrino stuff to Kaz form
         s_nu = s_nu_row(jj)/kerg*den_row(jj)
         rho_nu = rho_nu_row(jj)*den_row(jj)
         p_nu = p_nu_row(jj)

      end if
ccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     HELM with fake Kaz-like entropy
c     
ccccccccccccccccccccccccccccccccccccccccccccc
      if (OUTPUTTYPE.eq.3) then
c     
c     HELM total pressure
         p_tot = ptot_row(jj)
         u_tot = etot_row(jj)*den_row(jj)
c     entropy density in erg/cc/K/kb = 1/cc
         s_tot = (stot_row(jj)/kerg)*den_row(jj)

c     Convert neutrino stuff to Kaz form
         s_nu = s_nu_row(jj)/kerg*den_row(jj)
         rho_nu = rho_nu_row(jj)*den_row(jj)
         p_nu = p_nu_row(jj)
      end if
ccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     fake Kaz-like
c     
ccccccccccccccccccccccccccccccccccccccccccccc
      if (OUTPUTTYPE.eq.1) then
c     
c     Kaz-like pressure (replace nucleon EOS)
         p_tot = ptot_row(jj) - p_N + pkaz_N
         u_tot = etot_row(jj)*den_row(jj) - u_N + ukaz_N
c     1/cc
         s_tot = (stot_row(jj)/kerg)*den_row(jj) - s_N + skaz_N

c     Convert neutrino stuff to Kaz form
         s_nu = s_nu_row(jj)/kerg*den_row(jj)
         rho_nu = rho_nu_row(jj)*den_row(jj)
         p_nu = p_nu_row(jj)
      end if

ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     PWF99 for all
c     
ccccccccccccccccccccccccccccccccccccccccccccc
      if (OUTPUTTYPE.eq.2) then
c     PWF pressure (replace all)
         p_tot = ppwf_ele + ppwf_rad + ppwf_N
         
         u_tot = upwf_ele + upwf_rad + upwf_N
         
         s_tot = spwf_ele + spwf_rad + spwf_N

c     Convert neutrino stuff to Kaz form
         s_nu = s_nu_row(jj)/kerg*den_row(jj)
         rho_nu = rho_nu_row(jj)*den_row(jj)
         p_nu = p_nu_row(jj)
      end if


c     write(6,100) rhob, tk, p_tot, u_tot, s_tot, Qm












 100  format(40(1pe26.15,' '))
 10   format(40(I6,' '))

      return
      end












c..
c.. Compute neutrino emission rates using Kaz pair capture rates
c
      subroutine compute_kazlikenurates(jj, tkinput, rhobinput, Qm)
      implicit none
      save

c     Passed
      integer jj
      real*8 tkinput,rhobinput

c     Returned
      real*8 Qm

c     Local
      real*8 nb,neleposi,myye,yefit2
      real*8 etae,xnuc
      real*8 Qtot,QNetot
      real*8 xnuccalc

c..delcare local/global variables
c      include 'const.dek'
      include 'vector_eos.dek'
c      include 'kazeos.dek'
c      include 'kaz_state.dek'
c..
      include 'const.dek'



c      write(6,10) jj




ccccccccccccccccccccccccccccccccccccccccccccccc
c MAKE LIKE KAZ EOS OUTPUT
c
c
c  Compute neutrino emission rates
c

c..get the neutrino losses (does not include pair capture nor electron degeneracy effects)
      call sneut5(tkinput,rhobinput,abar_row(jj),zbar_row(jj),
     1            snu,snudt,snudd,snuda,snudz)

      xnuc=xnuccalc(rhobinput,tkinput)
  
c     \eta_e calculation should be at least as good as by solving non-linear equation in Kohri & Mineshige (2002)
c     Assume has rest-mass of electrons
       etae=etaele_row(jj)
c     Compute Yp,xnuc,npratiofree,etc. (npratiofree seems wrong in degenerate case since ~1E3 -- overestimation)
       call state_calc(rhobinput,tkinput,xnuc)

ccccccccccccccccccccccccccccccccccccccccccccc
c
c      Itoh+Kaz2002 for neutrino emission rates
c
ccccccccccccccccccccccccccccccccccccccccccccc
c      erg/s/cc
c      GODMARK: Is "nefree" just xne_row(jj)?
       call kazQm_calc(rhobinput,tkinput,xnuc,etae,Qtot,QNetot)
c      add-in floor to Qm (units should imply negligible!) : GODMARK
       Qm = snu*rhobinput + QNetot + 1E-30
       
       return
       end
