c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       force.inc
c    type:         include file
c    author:       f. douglas swesty
c    date:         2/29/92
c
c    purpose:      this include file contains the statement function
c                  definitions for the nuclear force expressions.  to
c                  change the nucleon-nucleon interaction only this file
c                  need be changed.  note:  we have assumed that the
c                  interaction is only density dependent, and that the
c                  interaction has no explicit temperature dependence!
c                  to introduce a temperature dependent interaction will
c                  require modification of the temperature derivatives
c                  in the routines nuceos & alfeos.
c
c
c    call line:    include 'force.inc/list'
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c                          nucleon-nucleon interaction parameters
      double precision aa, bb, cc, dd, dd3, scrdd3
c
c                          compression modulus, symmetry energy,
c                          binding energy, & saturation density
      double precision comp, symm, bind_e, nsubs
c
c                          surface symmetry energy & surface tension
      double precision sym_s, sig_s
c
c                          numerical coefficient (called alpha in
c                          the nucl. phys. a, vol. 535, pg. 331 paper)
      double precision skyrmc
c
c
c                   this common block contains the variables that
c                   specify the particular interaction used
      common /skyrme/ aa, bb, cc, dd, dd3, scrdd3,
     1                comp, symm, skyrmc, bind_e,
     2                nsubs, sym_s, sig_s
c
      double precision pv_e, dpvedn, dpvedx
      double precision pv_pr, dpvrdp, dpvrdn
      double precision pvp, dpvpdp, dpvpdn, dvp_dx, dvp_di
      double precision pvn, dpvndp, dpvndn, dvn_dx, dvn_di
      double precision denom, veterm, vterm1, vterm2
      double precision prot_d, nut_d, densit, nnn, xxx
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c                 these terms enter in many of the equations
c-----------------------------------------------------------------------
c
      denom(densit)=1.0+dd3*(densit**(dd-1.0))
c
      veterm(densit)= (cc*(1.0+dd)*(densit**dd)+
     1     2.0*cc*dd3*(densit**(2.0*dd-1.0)))/
     2     (denom(densit)**2)
c
c-----------------------------------------------------------------------
c                        internal energy stuff
c-----------------------------------------------------------------------
c
c                 internal energy (ie) due to interaction
      pv_e(prot_d,nut_d)=aa*((prot_d+nut_d)**2)+
     1     4.0*bb*prot_d*nut_d+( cc*((prot_d+nut_d)**(1.0+dd))/
     2     denom(prot_d+nut_d) )+prot_d*deltam
c
c                 derivative of ie w.r.t. density at fixed x
      dpvedn(nnn,xxx)=2.0*nnn*
     1    (aa+4.0*bb*xxx*(1.0-xxx))+veterm(nnn)+xxx*deltam
c
c                 derivative of ie w.r.t. x at fixed density
      dpvedx(nnn,xxx)=4.0*bb*(1.0-2.0*xxx)*(nnn**2)+
     1    nnn*deltam
c
c-----------------------------------------------------------------------
c                        interaction potential stuff
c-----------------------------------------------------------------------
c
c                 proton & neutron interaction potentials
      pvp(prot_d,nut_d)=2.0*aa*(prot_d+nut_d)+
     1     4.0*bb*nut_d+veterm(prot_d+nut_d)+deltam
c
      pvn(prot_d,nut_d)=2.0*aa*(prot_d+nut_d)+
     1     4.0*bb*prot_d+veterm(prot_d+nut_d)
c
c                 these terms enter in all of the equations
      vterm1(densit)=cc*(
     1    dd*(1.0+dd)*(densit**(dd-1.0))+
     2    2.0*dd3*(2.0*dd-1.0)*(densit**(2.0*dd-2.0)) 
     3    )/(denom(densit)**2)
      vterm2(densit)=-2.0*cc*dd3*(
     1    (dd**2-1.0)*(densit**(2.0*dd-2.0))+
     2    2.0*dd3*(dd-1.0)*(densit**(3.0*dd-3.0)) 
     3    )/(denom(densit)**3)
c
c                 derivative of proton potential w.r.t.
c                 proton density at fixed neutron density
      dpvpdp(prot_d,nut_d)=2.0*aa+
     1    vterm1(prot_d+nut_d)+vterm2(prot_d+nut_d)
c
c                 derivative of proton potential w.r.t.
c                 neutron density at fixed proton density
      dpvpdn(prot_d,nut_d)=2.0*aa+4.0*bb+
     1    vterm1(prot_d+nut_d)+vterm2(prot_d+nut_d)
c
c                 derivative of proton potential w.r.t.
c                 x at fixed density
cc      dvp_dx(nnn,xxx)=-4.0*bb*nnn
c
c                 derivative of proton potential w.r.t.
c                 density at fixed x
cc      dvp_di(nnn,xxx)=2.0*aa+4.0*bb*(1.0-xxx)+
cc     1    vterm1(nnn)+vterm2(nnn)
c
c                 derivative of neutron potential w.r.t.
c                 neutron density at fixed proton density
      dpvndn(prot_d,nut_d)=2.0*aa+
     1    vterm1(prot_d+nut_d)+vterm2(prot_d+nut_d)
c
c                 derivative of neutron potential w.r.t.
c                 proton density at fixed neutron density
      dpvndp(prot_d,nut_d)=2.0*aa+4.0*bb+
     1    vterm1(prot_d+nut_d)+vterm2(prot_d+nut_d)
c
c                 derivative of neutron potential w.r.t.
c                 x at fixed density
cc      dvn_dx(nnn,xxx)=4.0*bb*nnn
c
c                 derivative of neutron potential w.r.t.
c                 density at fixed x
cc      dvn_di(nnn,xxx)=2.0*aa+4.0*bb*xxx+
cc     1    vterm1(nnn)+vterm2(nnn)
c
c-----------------------------------------------------------------------
c                        pressure potential stuff
c-----------------------------------------------------------------------
c
c                 interaction contribution to the pressure (ip)
      pv_pr(prot_d,nut_d)=prot_d*pvp(prot_d,nut_d)+
     1    nut_d*pvn(prot_d,nut_d)-pv_e(prot_d,nut_d)
c
c                 derivative of ip w.r.t. proton density at
c                 fixed neutron density
      dpvrdp(prot_d,nut_d)=prot_d*
     1    dpvpdp(prot_d,nut_d)+nut_d*dpvndp(prot_d,nut_d)
c
c                 derivative of ip w.r.t. neutron density at
c                 fixed proton density
      dpvrdn(prot_d,nut_d)=prot_d*
     1    dpvpdn(prot_d,nut_d)+nut_d*dpvndn(prot_d,nut_d)
