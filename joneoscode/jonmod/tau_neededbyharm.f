

c     f2c -f -P tau_neededbyharm.f ; cp tau_neededbyharm.c tau_neededbyharm.P ~/latestcode


c     Duplicates those operations from tau_calc() that involve things dependent upon hcm so that hcm isn't independent variable in table
C======================================================================
c     ENSURE this list is the same for how used in tau_calc.f
      subroutine computefinal_fromhcm(
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
     1     ,Enu,Enue,Enuebar ! more outputs
     1     )
      
C======================================================================
      implicit none

      real*8 rate_2stream,density_2stream,tdifffromlambda ! functions


c     Passed to here
      real*8 clight,mb
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
     1     ,Enu,Enue,Enuebar ! more outputs

c     Locals
      real*8 H
      real*8 u_photon0
      real*8 u_nue0,u_nuebar0,u_numu0,u_nutau0
      real*8 n_nue0,n_nuebar0,n_numu0,n_nutau0
      real*8 n_nueth0,n_nuebarth0
      real*8 qtauttauohcm,qtauatauohcm
      real*8 ntauttauohcm,ntauatauohcm

      real*8 tauphoton,tauphotonabs

      real*8 Qm_nue,Qm_nuebar,Qmel,Qmmu,Qmtau
      real*8 Nm_nue,Nm_nuebar,Nmel,Nmmu,Nmtau
      real*8 Nm

      real*8 qtaut_nue ,
     1 qtaua_nue ,
     1 qtaut_nuebar ,
     1 qtaua_nuebar ,
     1 qtautmu ,
     1 qtauamu ,
     1 qtauttau ,
     1 qtauatau ,
     1 ntaut_nue ,
     1 ntaua_nue ,
     1 ntaut_nuebar ,
     1 ntaua_nuebar ,
     1 ntautmu ,
     1 ntauamu ,
     1 ntauttau ,
     1 ntauatau

      real*8 u_nue ,
     1 u_nuebar ,
     1 u_nuel ,
     1 u_numu ,
     1 u_nutau ,
     1 n_nue ,
     1 n_nuebar ,
     1 n_nuel ,
     1 n_numu ,
     1 n_nutau ,
     1 n_nueth ,
     1 n_nuebarth ,
     1 n_nuelth

      real*8 nb
      real*8 SMALL



c     Set some things
      nb = rhob/mb
      SMALL = 1D-50


c     Set some equivalences
      H=hcm
      u_photon0=unumu0

      u_nue0=unue0
      u_nuebar0=unuebar0
      u_numu0=unumu0
      u_nutau0=unumu0

      n_nue0=nnue0
      n_nuebar0=nnuebar0
      n_numu0=nnumu0
      n_nutau0=nnumu0

      n_nueth0=nnueth0
      n_nuebarth0=nnuebarth0

      qtauttauohcm=qtautmuohcm
      qtauatauohcm=qtauamuohcm


c     Set 2-stream approximation solutions



cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Photon emission energy-volume-rate
      tauphoton = tauphotonohcm*hcm
      tauphotonabs = tauphotonabsohcm*hcm
      Qphoton=rate_2stream(clight,hcm,u_photon0,tauphoton,tauphotonabs)


cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Neutrino emission energy-volume-rate

      qtaut_nue = qtautnueohcm*hcm
      qtaua_nue = qtauanueohcm*hcm
      qtaut_nuebar = qtautnuebarohcm*hcm
      qtaua_nuebar = qtauanuebarohcm*hcm
      qtautmu = qtautmuohcm*hcm
      qtauamu = qtauamuohcm*hcm
      qtauttau = qtauttauohcm*hcm
      qtauatau = qtauatauohcm*hcm

c     energy rates
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

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Neutrino emission number-volume-rate

c     Define similar quantities
      ntauttauohcm = ntautmuohcm
      ntauatauohcm = ntauamuohcm

      ntaut_nue = ntautnueohcm*hcm
      ntaua_nue = ntauanueohcm*hcm
      ntaut_nuebar = ntautnuebarohcm*hcm
      ntaua_nuebar = ntauanuebarohcm*hcm
      ntautmu = ntautmuohcm*hcm
      ntauamu = ntauamuohcm*hcm
      ntauttau = ntauttauohcm*hcm
      ntauatau = ntauatauohcm*hcm

c     number rates (in 1 direction)
      Nm_nue   =rate_2stream(clight,hcm,n_nue0,ntaut_nue,ntaua_nue)
      Nm_nuebar=rate_2stream(clight,hcm,n_nuebar0,ntaut_nuebar,ntaua_nuebar)
      Nmel = Nm_nue + Nm_nuebar

      Nmmu=rate_2stream(clight,hcm,n_numu0,ntautmu,ntauamu)
      Nmtau=rate_2stream(clight,hcm,n_nutau0,ntauttau,ntauatau)

      Nm=Nmel+Nmmu+Nmtau

c     Per volume (of size H^3) change:
c     Per area is obtained as mb*(Nm_nuebar-Nm_nue)*hcm
c     and then this is flux per face
c     Can be positive or negative
      graddotrhouye = mb*(Nm_nuebar - Nm_nue)


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

c     thermal number densities (using standard optical depth rather than thermalizing optical depth) GODMARK -- approximation to avoid tabulating thermalizing optical depths separately -- somewhat inconsistent with using lambdaintot for Tthermaltot
      n_nueth =density_2stream(n_nueth0,ntaut_nue,ntaua_nue)
      n_nuebarth =density_2stream(n_nuebarth0,ntaut_nuebar,ntaua_nuebar)
      n_nuelth = n_nueth+n_nuebarth ! diag only


c     total internal energy density
      rho_nu = u_nuel+u_numu+u_nutau
c     total pressure
      p_nu = rho_nu/3.d0
c     entropy density in 1/cc
      s_nu=4.d0/3.d0*rho_nu/(kbtk+SMALL) ! now in 1/cc
c     [k_b] form is s_nu/(rho/m_b) is "per baryon" entropy for rho and S given in cgs.

c     Electron neutrino fraction (other species cancel since their chemical potential is 0)
c     Should always be larger than 0.  For convergence process, let smallest be SMALL
c     can be negative due to optical depth supression
      Ynulocal = (n_nue-n_nuebar)/nb ! here this is just a diagnostic check
c     Could now check of Ynulocal is same as expected Ynu

c     Thermalized faction Y_\nu
c     can be negative due to optical depth supression
      Ynuthermal = (n_nueth-n_nuebarth)/nb


      Tdifftot =tdifffromlambda(clight,H,lambdatot)
      Tthermaltot =tdifffromlambda(clight,H,lambdaintot)


c     These are energies of *escaping* neutrinos
c     These energies are used for neutrino annihilation heating rates
      Enu = dmax1(Qm/(Nm+SMALL),0.0d0)
      Enue = dmax1(Qm_nue/(Nm_nue+SMALL),0.0d0)
      Enuebar = dmax1(Qm_nuebar/(Nm_nuebar+SMALL),0.0d0)


      return
      end
















c     Very similar to above but only computes energy densities (i.e. rho_nu, p_nu, s_nu)
C======================================================================
      subroutine computefinal_justdensities_fromhcm(
     1     clight,mb
     1     ,rhob,kbtk,hcm
     1     ,unue0,unuebar0,unumu0
     1     ,qtautnueohcm, qtautnuebarohcm, qtautmuohcm
     1     ,qtauanueohcm, qtauanuebarohcm, qtauamuohcm
     1     ,rho_nu,p_nu,s_nu ! outputs
     1     )
      
C======================================================================
      implicit none

      real*8 density_2stream ! functions


c     Passed to here
      real*8 clight,mb
     1     ,rhob,kbtk,hcm
     1     ,unue0,unuebar0,unumu0
     1     ,qtautnueohcm, qtautnuebarohcm, qtautmuohcm
     1     ,qtauanueohcm, qtauanuebarohcm, qtauamuohcm
     1     ,rho_nu,p_nu,s_nu

c     Locals
      real*8 H
      real*8 u_nue0,u_nuebar0,u_numu0,u_nutau0
      real*8 qtauttauohcm,qtauatauohcm

      real*8 qtaut_nue ,
     1 qtaua_nue ,
     1 qtaut_nuebar ,
     1 qtaua_nuebar ,
     1 qtautmu ,
     1 qtauamu ,
     1 qtauttau ,
     1 qtauatau

      real*8 u_nue ,
     1 u_nuebar ,
     1 u_nuel ,
     1 u_numu ,
     1 u_nutau

      real*8 nb
      real*8 SMALL



c     Set some things
      nb = rhob/mb
      SMALL = 1D-50


c     Set some equivalences
      H=hcm

      u_nue0=unue0
      u_nuebar0=unuebar0
      u_numu0=unumu0
      u_nutau0=unumu0

      qtauttauohcm=qtautmuohcm
      qtauatauohcm=qtauamuohcm


c     Set 2-stream approximation solutions



cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Neutrino emission energy-volume-rate

      qtaut_nue = qtautnueohcm*hcm
      qtaua_nue = qtauanueohcm*hcm
      qtaut_nuebar = qtautnuebarohcm*hcm
      qtaua_nuebar = qtauanuebarohcm*hcm
      qtautmu = qtautmuohcm*hcm
      qtauamu = qtauamuohcm*hcm
      qtauttau = qtauttauohcm*hcm
      qtauatau = qtauatauohcm*hcm

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

c     total internal energy density
      rho_nu = u_nuel+u_numu+u_nutau
c     total pressure
      p_nu = rho_nu/3.d0
c     entropy density in 1/cc
      s_nu=4.d0/3.d0*rho_nu/(kbtk+SMALL) ! now in 1/cc
c     [k_b] form is s_nu/(rho/m_b) is "per baryon" entropy for rho and S given in cgs.


      return
      end








C=====================================================================
      real*8 function rate_2stream(clight,H,density,tautot,tauabs)

c     2-stream approximation for volume rate
C=====================================================================
      implicit none
      
C=====================================================================
      real*8 clight,H,density,tautot,tauabs
      real*8 prefactor,osqrt3
      real*8 SMALL
C=====================================================================

      SMALL = 1D-50
c      SMALL = 0.0d0

c     /H means rates are per unit volume
      prefactor = (clight/4.0)/(3.0/4.0)/(H+SMALL)
      osqrt3=0.5773502691896

      rate_2stream=prefactor*density/(0.5d0*tautot + osqrt3 + 1.d0/(3.d0*tauabs+SMALL)+SMALL)

      return
      end


C=====================================================================
      real*8 function density_2stream(density,tautot,tauabs)

c     2-stream approximation for density
c     As tau->0, final density->0
c     As tau->\infty, final density -> density
C=====================================================================
      implicit none
      
C=====================================================================
      real*8 density,tautot,tauabs
      real*8 twosqrt3
      real*8 SMALL
C=====================================================================

      SMALL = 1D-50
c      SMALL = 0.0d0

c      twosqrt3 = 2.0*sqrt(3.0)
      twosqrt3 = 3.4641016151377545

      density_2stream = density*(3.0*tautot + twosqrt3)/(3.0*tautot + twosqrt3+2.0/(tauabs+SMALL)+SMALL)

      return
      end


C=====================================================================
      real*8 function tdifffromlambda(clight,H,lambda)

c     diffusion time limited by speed of light
C=====================================================================
      implicit none
      
C=====================================================================
      real*8 clight,H,lambda
      real*8 SMALL
C=====================================================================

      SMALL = 1D-50
c      SMALL = 0.0d0


      tdifffromlambda=dmax1(3.0*H*H/(SMALL+1.0*clight*lambda),H/clight)


      return
      end


