

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c BELOW need never be touched by user
c Compared to original code, see JCM comments.
c JCM changed constants (used my const.dek)
c
c  Note that original LS code just called el_eos()
c  Their new wrapper, however, calls any_electron(), which by default calls the TIMMS EOS
c  
c
c
c
c      call any_electron(t,ye,brydns) 
c      call el_eos(t,ye,brydns)
c
c
c  any_electron has been moved to jon_lsbox.f so can make code more general and keep electron parts separate from nuclear parts
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc








c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c
c              lattimer & swesty eos (release # 2.7  9/1/95)
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         inveos
c    module:       inveos
c    type:         subroutine
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         5/23/90
c                  bug fixed on (5/24/90) (affected only performance
c                  of code not the results!)
c
c
c    call line:    call inveos(inpvar,t_old,ye,brydns,iflag,eosflg,xprev)
c
c    inputs:       inpvar = temp, internal eng, or entropy
c                  t_old = initial guess at the temperature
c                  ye = electron fraction
c                  brydns = baryon number density
c                  iflag = 1 --> inpvar is temperature
c                          2 --> inpvar is internal energy
c                          3 --> inpvar is entropy (not implem)
c
c    outputs       eosflg = 1 --> "no nuclei" eos
c                           2 --> general eos
c                           3 --> bulk eos for densities above nuclear
c                  xprev = unused
c                  p_prev = previous value of proton density
c
c
c
c 
c    include files:  eos_m4c.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine inveos(inpvar,t_old,ye,brydns,iflag,
     1                  eosflg,forflg,sf,xprev,p_prev)
c
c
c

c
      implicit none
c
      include 'eos_m4c.inc'
c
c                         local variables
c
      double precision inp_v, inp_vo, inp_vn, uftn, duftn, dt
      double precision t_old, t_new, t_temp, t_lb, t_ub, perdif
      integer loop, sf, new_f
c
      rsflag = 1
c                         input is the temperature; call the eos
c                         normally and then return
      if(iflag.eq.1) then
        call eos_m4c(inpvar,ye,brydns,1,eosflg,forflg,sf,
     1 xprev,p_prev)
        t_old = inpvar(1)
        return
      endif
c
c
c                         the input variable must be the internal
c                         energy so calc the internal energy for
c                         the initial guess at the temperature
        inp_v = inpvar(1)
c
        t_lb = 0.15
        t_ub = 50.0
c
        inpvar(1) = t_old
        call eos_m4c(inpvar,ye,brydns,1,eosflg,forflg,sf,
     1 xprev,p_prev)
ccc      call eos_m1d(t_old,ye,brydns,1,eosflg,xprev,p_prev)
c
c                         save the value of the internal energy
      if(iflag.eq.2) then
        inp_vo = utot
      elseif(iflag.eq.3) then
        inp_vo = stot
      endif
c
c
c                         tweak the initial guess slightly so as to
c                         get a new value of the internal energy
c
      t_new = 1.1*t_old
c
      new_f = 1
c
c JCM: Try more loops
c      do 20 loop=1,50,1
      do 20 loop=1,200,1

c JCM:
c         write(*,*) 'looping',loop
c
        inpvar(1) = t_new
        call eos_m4c(inpvar,ye,brydns,1,eosflg,forflg,sf,
     1 xprev,p_prev)
ccc        call eos_m1d(t_new,ye,brydns,1,eosflg,xprev,p_prev)
c
c
        if(sf.ne.1.and.new_f.eq.1) then
cc          write(*,*) 'inveos: eos fatally failed at try:'
cc          write(*,*) t_new,brydns,ye
          t_new = t_new-0.25
          t_lb = dmin1(t_lb,t_new-1.0d-1)
          goto 20
        elseif(sf.ne.1) then
           dt = 0.5*dt
           t_new = t_new+dt
        else
c
          new_f = 0
c
c                         save this value of the internal energy too
          if(iflag.eq.2) then
            inp_vn = utot
          elseif(iflag.eq.3) then
            inp_vn = stot
          endif
c
          if(inp_vn.lt.inp_v) then
            t_lb = t_new
c            write(*,*) 'l @ ',t_new,inp_vn,inp_v
          elseif(inp_vn.gt.inp_v) then
c            write(*,*) 'u @ ',t_new,inp_vn,inp_v
            t_ub = t_new
          endif
        endif
c
        uftn = inp_vn-inp_v
c
        if(loop.lt.20) then
c                         this is the function to be zeroed by the
c                         newton-raphson iteration
c
c                         numerical derivative of the above function
c                         w.r.t. the temperature
cc          duftn = ((inp_vn-inp_vo)/(t_new-t_old))+1.0d-15
c
c                         analytic derivatives
          if(iflag.eq.2) then
            duftn = dudt
          elseif(iflag.eq.3) then
            duftn = dsdt
          endif
c
c                         estimated correction to temperature
          dt = uftn/duftn
c
c
 10       continue
c                         temporarily store the new temperature
          t_temp = t_new-dt
c
c                         is the new temp within a valid range?
          if((t_temp.gt.t_lb).and.(t_temp.lt.t_ub)) then
            t_old = t_new
            inp_vo = inp_vn
            t_new = t_temp
          else
c                         if not cut the step size in half & try again
            dt = 0.5*dt
            if(t_temp.eq.t_new) then
              t_old = t_new
              t_new = 0.5*(t_lb+t_ub)
              dt = t_new-t_old
            else
              goto 10
            endif
          endif
c
        else
c
          t_old = t_new
          t_new = 0.5*(t_lb+t_ub)
          dt = t_new-t_old
c
        endif
c
c                         if relative change in t is less than 1.0e-5
c                         then quit out of loop
        if(abs(dt/t_new).lt.lstol1) goto 30
c
c                end of the do loop
 20   continue
c
c                didn't meet convergence criterion
      sf = 0
      write(*,*) ' inversion of eos failed to converge',dt,t_new
c
c                met the convergence criterion!!!
 30   continue
c
c                this stuff is commented out for speed reasons; it
c                virtually never gets tripped anyway
cc      inpvar(1) = t_old
cc      call eos_m4c(inpvar,ye,brydns,1,eosflg,forflg,sf,xprev,p_prev)
cc      if(iflag.eq.2) then
cc        perdif = inp_v-utot
cc      else
cc        perdif = inp_v-stot
cc      endif
cc      if(abs(perdif).gt.1.0d-4) then
cc        write(*,*) 'inveos: failure',inp_v,stot
cc        write(*,*) uftn,dt,loop
cc        write(*,*) 'try:',t_new,brydns,ye
cc        sf = 0
cc        return
cc      endif
c 
c                return this value for t
      inpvar(1) = inp_v
      t_old = t_new
c
c                time to call it quits!
 999  return
c
      end
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       eos_m4c
c    type:         subroutine
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         3/3/92  model 4c modifications completed
c                  12/15/90 modified from model 4a to include the
c                  phase boundary cutoffs and maxwell construction
c                  boundaries.
c                  7/13/90 modified from model 1-d to include maxwell
c                  construction
c                  5/25/90  model 1d
c
c                  please report any problems to me at:
c                  bitnet:  swesty@sunysbnp or
c                  internet: fswesty@astro.sunysb.edu or
c                            fswesty@sbast3.sunysb.edu
c
c
c    call line:    call eos_m4c(inpvar,ye,brydns,iflag,eosflg,fflag,
c                  xprev,p_prev)
c
c    inputs:       inpvar = temp, internal eng, or entropy
c                  ye = electron fraction
c                  brydns = baryon number density
c                  iflag = 1 --> inpvar is temperature
c                          2 --> inpvar is internal energy (not implem)
c                          3 --> inpvar is entropy (not implemented)
c                          (iflag=1 is now assumed at this level)
c                  fflag = "forcing flag"  0 --> no forcing
c                                          1 --> force a particular
c                                                scheme to be used
c
c
c    outputs:      eosflg = 1 --> not implemented in model 4b
c                           2 --> general eos
c                           3 --> bulk eos (includes alpha's)
c                  xprev = previous value of x (must be supplied on
c                          first call)
c                  p_prev = previous value of proton density (must be
c                          supplied on first call)
c
c
c
c 
c    include files:  eos_m4c.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine eos_m4c(inpvar,ye,brydns,iflag,eosflg,fflag,ssflag,
     1                   xprev,p_prev)
c
      implicit none
c
      double precision outvar(4)
c
c
c                       this include file contains all variable
c                       declarations.  note:: no implicit typing
c                       scheme is followed in this code; if you
c                       have any doubt as to a variables type check
c                       it!!!!.  also note that all variables are
c                       declared explicitly.
c
      include 'eos_m4c.inc'
c
c
c                         set the "switch" flag to zero
      swtflg = 0
c
c                         set t equal to the input variable (the entropy
c                         and internal energy options should go through
c                         inveos untill further notice)
      t = inpvar(1)
c
c
c                         if the "forcing" flag is set then skip
c                         the eos determination logic and go straight
c                         to the eos determined by eosflg
      if(fflag.eq.1) then
        goto 10
      else
c                         otherwise let the eos logic module determine
c                         the correct eos to use
        call eoslog(inpvar,ye,brydns,eosflg)
      endif
c
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                        try nuceos first and if not successfull
c                        then try bulk eos
 10   continue
      if(eosflg.eq.1) then
c
        call nuceos(inpvar,ye,brydns,xprev,p_prev,ssflag)
c
c                    if the nuclear eos failed and the reset flag is set
c                    then reset the initial guesses and try again
        if((ssflag.ne.1).and.(rsflag.eq.1)) then
          call reset(inpvar,ye,brydns,outvar)
          outvar(1) = inpvar(1)
          call nuceos(outvar,ye,brydns,xprev,p_prev,ssflag)
c
c
c                    make a last ditch effort at convergence
          if(ssflag.ne.1) then
            outvar(2) = 0.155
            outvar(3) = -15.0
            outvar(4) = -20.0
            call nuceos(outvar,ye,brydns,xprev,p_prev,ssflag)
          endif
c
        endif
c
c
c
        if((xh.gt.heavct).and.(ssflag.eq.1)) then
c                    set eos flag to full scheme
          eosflg = 2
c
c                    else if fraction of nuclei is less than the minimum
c                    or if nuceos was unsuccessful use the no nuclei eos
        else
          if(fflag.ne.1) then
c
            call alfeos(inpvar,ye,brydns,p_prev,ssflag)
c
            if((ssflag.ne.1).and.(fflag.eq.1)) then
              eosflg = 1
              write(*,*) 'a2 failed at try = ',t,brydns,ye
              goto 999
            endif
c
c                    set nuclei to bulk eos
            eosflg = 3
c                    save value of proton fraction
            p_prev = ye*brydns
c
            goto 999
c
          else
            if(nf_flg.eq.1) 
     1          write(*,*) 'nuc failed at t,rho = ',t,brydns
            goto 999
          endif
        endif
c
      endif
c
c
c          end of nuceos--bulk eos calculations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c
c
c
c
c
c
c
c
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                            calculate full eos (including nuclei)
      if(eosflg.eq.2) then
c
c                    call the nuclear eos
        call nuceos(inpvar,ye,brydns,xprev,p_prev,ssflag)
c
c
c                    if the nuclear eos failed and the reset flag is set
c                    then reset the initial guesses and try again
        if((ssflag.ne.1).and.(rsflag.eq.1)) then
cccc          write(*,*) ' eos_m4c:: r.i.gs.'
          call reset(inpvar,ye,brydns,outvar)
          outvar(1) = inpvar(1)
          call nuceos(outvar,ye,brydns,xprev,p_prev,ssflag)
c
c
c                    make a last ditch effort at convergence
          if(ssflag.ne.1) then
            outvar(2) = 0.155
            outvar(3) = -15.0
            outvar(4) = -20.0
            call nuceos(outvar,ye,brydns,xprev,p_prev,ssflag)
          endif
c
c
c
          if(ssflag.ne.1) then
cccc            write(*,*) '     r.i.gs. failure @ try: ',inpvar
            goto 999
          else
            inpvar(2) = outvar(2)
            inpvar(3) = outvar(3)
            inpvar(4) = outvar(4)
          endif
c                    otherwise quit and return
        elseif((ssflag.ne.1).and.(fflag.eq.1)) then
          goto 999
        endif
c
c
c
c
c                    if fraction of heavies is greater than the minimum
c                    parameter, then this eos is ok
        if((xh.gt.heavct).and.(ssflag.eq.1)) then
c                    set eos flag to full scheme
          eosflg = 2
c
c                    else if fraction of nuclei is less than the minimum
c                    or if nuceos was unsuccessful use the no nuclei eos
        else
c                    if the forcing flag is not set
          if(fflag.ne.1) then
c                    set nuclei to no nuclei eos
            eosflg = 3
c                    set flag to indicate switch is being made
            swtflg = 1
c
            write(*,*) ' nuceos failed at try =',t,brydns,ye
            write(*,*) ' where it shouldnt have; bulk eos was used'
            write(*,*) ' iv = ',inpvar
            write(*,*) ' '
c
c                    branch to bulk eos
            goto 50
c
c                    otherwise since forcing flag is set then declare
c                    a failure and return
          else
c                      if the failure message flag is set then announce
c                      the failure
            if(nf_flg.eq.1) 
     1          write(*,*) 'nuc failed at t,r = ',t,brydns
            goto 999
          endif
        endif
c
      endif
c                              end of full eos calulations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c
c
c
c
c
c
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                              calculate bulk eos
 50   continue
      if(eosflg.eq.3) then
c
        call alfeos(inpvar,ye,brydns,p_prev,ssflag)
c
        if((ssflag.eq.0).and.(fflag.eq.1).and.(nf_flg.eq.1)) then
          write(*,*) 'a1 failed at t,rho = ',t,brydns
          goto 999
        endif
c                           if this eos was used as a result of the
c                           nuclear eos failing then set the
c                           success flag to indicate a warning
        if(swtflg.eq.1) then
          ssflag = 2
        endif
c
c                           save the value of the proton fraction
        p_prev = ye*brydns
c
        goto 999
c
      endif
c                end of bulk eos calculations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c
c
c
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                              calculate via maxwell construction
      if(eosflg.eq.4) then
c
        call maxwel(inpvar,ye,brydns,xprev,p_prev,ssflag)
c
c                 save the value of the proton fraction
        p_prev = ye*brydns
c
c                 if maxwell eos failed then announce the failure
        if(ssflag.ne.1) then
          write(*,*) ' maxwel failed at try = '
          write(*,*) t,brydns,ye
        endif
c
          goto 999
c
      endif
c                end of maxwell construction calculations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c
 999  return
c
      end
c
cnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnu
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         nuceos.for
c
c***********************************************************************
c
c    module:       nuceos
c    type:         subroutine
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         7/13/90 modified from model 1-d
c
c                  bitnet:  swesty@sunysbnp or
c                  internet: fswesty@astro.sunysb.edu or
c                            fswesty@sbast3.sunysb.edu
c
c    call line:    call nuceos(inpvar,ye,brydns,x_prev,ssflag)
c
c    inputs:       inpvar = temp, internal eng, or entropy
c                  ye = electron fraction
c                  brydns = baryon number density
c
c    outputs:      xprev = previous value of x (must be supplied on
c                  first call)
c                  ssflag = success flag 0 --> failure
c                                        1 --> success
c
c
c 
c    include files:  eos_m4c.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine nuceos(inpvar,ye,brydns,xprev,p_prev,ssflag)
c
      implicit none
c
c
      include 'eos_m4c.inc'
      include 'el_eos.inc'
c
c
c                       function type declarations
c
      double precision f_1_2, f_3_2, finv12, fhalfi, fhalfo
      double precision fhalf
c
      double precision zng, zpg
      integer tcflag, ftflag
c
      integer kki,lli
      double precision result(5), r_check(5)
      double precision a_tmp(5,5)
      double precision ni_min
c
      integer cflag, schflg
      double precision dtst1, dtst2
      double precision break, dnsi, dtmp8
      double precision dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6,dtmp7
cc      double precision tbsph, tbph, tbnh, tbspl, tbpl, tbnl
cc      double precision dbspdx, dbpdx, dbndx, dbspdu, dbpdu, dbndu
cc      double precision tsgl, tsgh, thl, thh, dsgdx, dhfdx, ds2dx,dzdx
cc      double precision dpt1dx, dpt2dx
c
      include 'force.inc'
c     JCM: another version has this included
      include 'const.dek'

c
c                         set the scheme flag to zero
      schflg = 0
c
c
 5    continue
c
c
c
c
c                         set t equal to the input variable (the entropy
c                         and internal energy options are not implemented
c                         in this version)
      t = inpvar(1)
      nsubi = inpvar(2)
      eta_po = inpvar(3)
      eta_no = inpvar(4)
c
c
c                         calc the quantum concentration of nucleons
      nq = 2.36d-4*t**1.5 
c
c                         calc the fermi integral coefficent
      uq = 20.721
c
      mq = (t/uq)**1.5
c
      kq = ((t/uq)**2.5)/(2.0*pi**2)
c
      lq = uq*(mq**ovr53)/(3.0*(pi**2))
c
      etamax = 0.95*finv12(2.0*(pi**2)*brydns/mq)
c
      if(eta_po.ge.etamax) eta_po = etamax-0.1
      if(eta_no.ge.etamax) eta_no = etamax-0.1
      ni_min = dmax1(4.5d-2,brydns)
      if(nsubi.lt.ni_min) nsubi = ni_min+1.0d-3
c
      tcflag = 0
c
      cflag = 0
c
      newflg = 1
c
c                    start newton-raphson iteration here
c
c
      do 30 i=1,maxit,1
 
c JCM:
c        write(*,*) 'looping-i',i,maxit
c
        it_num = i
c                       set the "negative" flag
        ngflag = 0
c
c
c
        nnout = mq*f_1_2(eta_no)/(2.0*pi**2)
        npout = mq*f_1_2(eta_po)/(2.0*pi**2)
c
        nout = npout+nnout
c
c20        vnout = eiflag*(2.0*aa*nout+4.0*bb*npout+cc*(1.0+dd)*nout**dd)
        vnout = eiflag*pvn(npout,nnout)
c
c20        vpout = eiflag*(2.0*aa*nout+4.0*bb*nnout+
c20     1    cc*(1.0+dd)*nout**dd+deltam)
        vpout = eiflag*pvp(npout,nnout)
c
        f32_no = f_3_2(eta_no)
c
        f32_po = f_3_2(eta_po)
c
c20        bprout = lq*(f32_po+f32_no)+eiflag*(
c20     1    aa*(nout**2)+4.0*bb*npout*nnout+dd*cc*(nout**(1.0+dd)) )
        bprout = lq*(f32_po+f32_no)+eiflag*pv_pr(npout,nnout)
c
        mun_o = t*eta_no+vnout
c
        mup_o = t*eta_po+vpout
c
        mualfa = 2.0*mun_o+2.0*mup_o+balpha-bprout*v_alfa
c
        if(abs(mualfa/t).lt.30.0) then
          alfdns = 8.0*nq*dexp(mualfa/t)
        elseif((mualfa/t).lt.-30.0) then
          alfdns = 0.0
        else
          alfdns = 8.0*nq*dexp(3.0d1)
        endif
c
c
c                   these statements take out the alfas if the
c                   alpha particle enable flag is not set
        if(alflag.ne.1) then
          alfdns = 0.0
          mualfa = -300.0
        endif
c
c
c
        exalfa = 1.0-alfdns*v_alfa
c
c
        bpralf = alfdns*t
c
c---------------------------------------------------
c
c
c             calculate fraction of space occupied by nuclei
        u_nuc = (brydns-exalfa*nout-4.0*alfdns)/
     1        (nsubi-exalfa*nout-4.0*alfdns)
c
c
c            is volume occupied by nuclei within acceptable limits?
cc        if((u_nuc.lt.0.0).or.((u_nuc-1.0).gt.-1.0e-20)) then
cc        if((u_nuc.lt.0.0).or.(u_nuc.gt.0.996)) then
        if((u_nuc.lt.1.0d-17).or.(u_nuc.gt.0.996)) then
          ngflag = 1
          goto 29
        endif
c
c
c            volume exclusion factor due to nuclei
        exclu = 1.0-u_nuc
c
c
c            if calculated nucleon and alfa densities are larger
c            than the baryon density then reduce the eta's
        if((exclu*exalfa*nout+exclu*4.0*alfdns).gt.brydns) then
          ngflag = 1
          goto 29
        endif
c
c
c            calculate the internal (inside nuclei) proton fraction
c
        x = (brydns*ye-(1.0-u_nuc)*(exalfa*npout+2.0*alfdns))/
     1    (u_nuc*nsubi)
        compx = 1.0-x
c
c
c            is x within reasonable (but not necessarily correct)
c            limits? (ye may not be the lower bound on x !!!)
cccc        x_min = dmax1(1.0d-2,(ye-0.05))
        x_min = dmax1(1.0d-2,(0.8*ye))
cc        x_min = 0.95*ye
        if((x.lt.x_min).or.(x.gt.0.6)) then
          ngflag = 1
          goto 29
        endif
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c                     calculate critical temperature & its x derivative
        tsc_12 = 87.76*((comp/375.0)**0.5)*((0.155/nsubs)**ovr3)
c
cccdebug      tsc_12 = 1.0d8
c
        tsubc = tsc_12*x*compx
        dtcdx = tsc_12*(1.0-2.0*x)
        dtcdxx = -2.0*tsc_12
        h = 1.0-2.0*(t/tsubc)**2+(t/tsubc)**4
c
cc        tsubc = tsc_12*0.25
cc        dtcdx = 0.0
cc        dtcdxx = 0.0
c
c
cc        tsubc = 80.0*x*compx
cc        dtcdx = 80.0*(1.0-2.0*x)
cc        dtcdxx = -160.0
c
c                     if the x is such that t is greater than the
c                     critical temperature then fix nsubi so that
c                     it lies in the bounds of acceptable parameter
c                     space
        ftflag = 0
        if(((t.gt.tsubc).or.(h.le.0.0)).and.(schflg.eq.0)) then
c                       if this is an initial guess, then lower
c                       nsubi untill we get a good x
          if(newflg.eq.1) then
cc        write(*,*) ' nuc exceeded tc'
cc        write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc
            zng = 2.0*(pi**2)*brydns*(1.0-0.1*ye)/(1.5*mq)
            zpg = 2.0*(pi**2)*brydns*0.1*ye/(1.5*mq)
            if(tcflag.ne.1) then
              tcflag = 1
              eta_po = finv12(zpg)-0.0
              eta_no = finv12(zng)-0.0
            else
              eta_po = eta_po-2.0/t
              eta_no = eta_no-2.0/t
              nsubi = dmax1(0.9*nsubi,5.1d-2)
            endif
            if(dbflag.eq.1) then
              write(*,2000) '1',i,nsubi,eta_po,eta_no,dnsubi
            endif
            goto 30
          else
c                       otherwise go back and cut the stepsize in
c                       half since it was obviously too big
            ngflag = 1
            goto 29
          endif
        elseif(((t.gt.tsubc).or.(h.le.0.0)).and.(schflg.eq.1)) then
          ftflag = 1
          tsubc = 80.0*(0.25+0.5*ye)*(0.75-0.25*ye)
c
        endif
c
c
        r_0 = (0.75/(pi*nsubs))**ovr3
        q = (384.0*pi*(r_0**2)*sig_s/sym_s)-16.0
c
c                        calculate surface functions of the internal
c                        (nuclear) proton fraction, x
        sigma = 1.0/(q+1.0/(x**3)+1.0/(compx**3))
        ovrx4 = (1.0/x**4)-(1.0/compx**4)
        dsigdx = 3.0*(sigma**2)*ovrx4
        sigsgp = dsigdx/sigma
        sigsg2 = 18.0*(sigma**2)*ovrx4**2-12.0*sigma*((1.0/x**5)+
     1  (1.0/compx**5))
c
c                        if t is less than critical temp then
        if(t.lt.tsubc) then
c                        calculate the surface energy temperature factor
c                        and its x and t derivatives
          h = 1.0-2.0*(t/tsubc)**2+(t/tsubc)**4
          hprim = -4.0*t/(tsubc**2)+4.0*((t/tsubc)**3)/tsubc
          hpprim = -4.0/(tsubc**2)+12.0*(t**2)/(tsubc**4)
          dhdx = 4.0*(t**2/tsubc**3-t**4/tsubc**5)*dtcdx
          dhdxx = 4.0*(t**2/tsubc**3-t**4/tsubc**5)*dtcdxx+
     1    4.0*(-3.0*t**2/tsubc**4+5.0*t**4/tsubc**6)*(dtcdx**2)
          hx = dhdx/h
          dhdtdx = 8.0*(t/tsubc**3-2.0*(t**3)/tsubc**5)*dtcdx
c
c
c                        x independent version of tzero
c          tzero = 0.25*tsc_12
c          dtzdx = 0.0
c          dtzdxx = 0.0
c                        x dependent version of tzero
c          tzero = tsubc
c          dtzdx = dtcdx
c          dtzdxx = dtcdxx
c
c
c
c                        coulomb liquid correction factors and their
c                        derivatives
c          w = 1-(t/tzero)**2
c          dwdx = 2.0*(t**2)*dtzdx/(tzero**3)
c          dwdt = -2.0*t/(tzero**2)
c          dwdtdx = 4.0*t*dtzdx/(tzero**3)
c          dwdxdx = 2.0*(t**2)*
c     1    (dtzdxx/(tzero**3)-3.0*(dtzdx**2)/(tzero**4))
c          dwdtdt = -2.0/(tzero**2)
c
          w = 1.0
          dwdt = 0.0
          dwdx = 0.0
          dwdtdx = 0.0
          dwdxdx = 0.0
          dwdtdt = 0.0
c
c
c
c                        calc lattice factor & derivatives & products
c
          exclu = 1.0-u_nuc
          compu = 1.0-u_nuc
c
          du = dmax1(1.0d-15, (1.0-1.5*w*u_nuc**ovr3+0.5*u_nuc))
          dmu = dmax1(1.0d-15,(1.0-1.5*w*(1.0-u_nuc+1.0e-20)**ovr3+
     1    0.5*(1.0-u_nuc)))
c
          dup = -0.5*w*u_nuc**m2ovr3+0.5
          dmup =-0.5*w*(1.0-u_nuc+1.0e-20)**m2ovr3+0.5
          dupp = ovr3*w*((u_nuc+1.0d-20)**m5ovr3)
          dmupp = ovr3*w*((1.0-u_nuc)+1.0e-20)**m5ovr3
c
c                derivatives w.r.t. t
c
          dut = -1.5*dwdt*u_nuc**ovr3
          dmut = -1.5*dwdt*(1.0-u_nuc+1.0e-20)**ovr3
          dupt = -0.5*dwdt*u_nuc**m2ovr3
          dmupt = -0.5*dwdt*(1.0-u_nuc+1.0e-20)**m2ovr3
c
c                derivatives w.r.t. x
c
          dux = -1.5*dwdx*u_nuc**ovr3
          dmux = -1.5*dwdx*(1.0-u_nuc+1.0e-20)**ovr3
          dupx = -0.5*dwdx*u_nuc**m2ovr3
          dmupx = -0.5*dwdx*(1.0-u_nuc+1.0e-20)**m2ovr3
c
c                second derivatives w.r.t. x
c
          duxx = -1.5*dwdxdx*u_nuc**ovr3
          dmuxx = -1.5*dwdxdx*(1.0-u_nuc+1.0e-20)**ovr3
c
c                second derivatives w.r.t. t
c
          dutt = -1.5*dwdtdt*u_nuc**ovr3
          dmutt = -1.5*dwdtdt*(1.0-u_nuc+1.0e-20)**ovr3
c
c                second derivatives w.r.t. x & t
c
          duxt = -1.5*dwdtdx*u_nuc**ovr3
          dmuxt = -1.5*dwdtdx*(1.0-u_nuc+1.0e-20)**ovr3
c
c
          tmp1 = (u_nuc**2)+(compu**2)+0.6*(u_nuc*compu)**2
          tmp1p = 4.0*u_nuc-2.0+
     1    2.0*0.6*(u_nuc*compu**2-compu*u_nuc**2)
          tmp1pp = 4.0+2.0*0.6*(compu**2-4.0*u_nuc*compu+u_nuc**2)
c
          tmp2 = compu*(du**ovr3)
          tmp2p = -1.0*du**ovr3+ovr3*compu*(du**m2ovr3)*dup
          tmp2pp = -ovr23*(du**m2ovr3)*dup-ovr29*compu*
     1    (du**m5ovr3)*dup**2+ovr3*compu*(du**m2ovr3)*dupp
c
          tmp2t = ovr3*compu*(du**m2ovr3)*dut
          tmp2x = ovr3*compu*(du**m2ovr3)*dux
          tmp2xx = ovr3*compu*(du**m2ovr3)*duxx+
     1        m2ovr3*ovr3*compu*(du**m5ovr3)*(dux**2)
          tmp2tt = ovr3*compu*(du**m2ovr3)*dutt+
     1        m2ovr3*ovr3*compu*(du**m5ovr3)*(dut**2)
          tmp2xt = ovr3*compu*(du**m2ovr3)*duxt+
     1        m2ovr3*ovr3*compu*(du**m5ovr3)*dux*dut
          tmp2pt = -ovr3*(du**m2ovr3)*dut+
     1        m2ovr3*ovr3*compu*(du**m5ovr3)*dup*dut+
     2        ovr3*compu*(du**m2ovr3)*dupt
          tmp2px = -ovr3*(du**m2ovr3)*dux+
     1        m2ovr3*ovr3*compu*(du**m5ovr3)*dup*dux+
     2        ovr3*compu*(du**m2ovr3)*dupx
c
c
c
          tmp3 = u_nuc*(dmu**ovr3)
          tmp3p = (dmu**ovr3)-ovr3*u_nuc*(dmu**m2ovr3)*dmup
          tmp3pp = -ovr23*(dmu**m2ovr3)*dmup-ovr29*u_nuc*
     1    (dmu**m5ovr3)*(dmup**2)+ovr3*u_nuc*(dmu**m2ovr3)*dmupp
c
          tmp3t = ovr3*u_nuc*(dmu**m2ovr3)*dmut
          tmp3x = ovr3*u_nuc*(dmu**m2ovr3)*dmux
          tmp3xx = ovr3*u_nuc*(dmu**m2ovr3)*dmuxx+
     1        m2ovr3*ovr3*u_nuc*(dmu**m5ovr3)*(dmux**2)
          tmp3tt = ovr3*u_nuc*(dmu**m2ovr3)*dmutt+
     1        m2ovr3*ovr3*u_nuc*(dmu**m5ovr3)*(dmut**2)
          tmp3xt = ovr3*u_nuc*(dmu**m2ovr3)*dmuxt+
     1        m2ovr3*ovr3*u_nuc*(dmu**m5ovr3)*dmux*dmut
          tmp3pt = ovr3*(dmu**m2ovr3)*dmut-ovr3*m2ovr3*u_nuc*
     1      (dmu**m5ovr3)*dmup*dmut-ovr3*u_nuc*(dmu**m2ovr3)*dmupt

          tmp3px = ovr3*(dmu**m2ovr3)*dmux-ovr3*m2ovr3*u_nuc*
     1      (dmu**m5ovr3)*dmup*dmux-ovr3*u_nuc*(dmu**m2ovr3)*dmupx
c
c
c                 combination d function
c
          scrdu = u_nuc*compu*(tmp2+tmp3)/tmp1
          scrdut = u_nuc*compu*(tmp2t+tmp3t)/tmp1
          scrdux = u_nuc*compu*(tmp2x+tmp3x)/tmp1
          scrdxx = u_nuc*compu*(tmp2xx+tmp3xx)/tmp1
          scrdtt = u_nuc*compu*(tmp2tt+tmp3tt)/tmp1
          scrdxt = u_nuc*compu*(tmp2xt+tmp3xt)/tmp1
c
          scrd = scrdu/u_nuc
          scrdt = scrdut/u_nuc
          scrdx = scrdux/u_nuc
c
          scrd2 = scrdu/compu
          scrd2t = scrdut/compu
          scrd2x = scrdux/compu
c
          scrdup = scrd-scrd2+u_nuc*compu*
     1    ((tmp2p+tmp3p)/tmp1-(tmp2+tmp3)*tmp1p/tmp1**2)
c
          scrdpt = scrdt-scrd2t+u_nuc*compu*
     1    ((tmp2pt+tmp3pt)/tmp1-(tmp2t+tmp3t)*tmp1p/tmp1**2)
c
          scrdpx = scrdx-scrd2x+u_nuc*compu*
     1    ((tmp2px+tmp3px)/tmp1-(tmp2x+tmp3x)*tmp1p/tmp1**2)
c
          scrdpp = (scrdup-scrd)/u_nuc-(scrd2+scrdup)/compu+
     1    (1.0-2.0*u_nuc)*
     2    ((tmp2p+tmp3p)/tmp1-(tmp2+tmp3)*tmp1p/tmp1**2)+u_nuc*compu*
     3    ((tmp2pp+tmp3pp)/tmp1-2.0*(tmp2p+tmp3p)*tmp1p/tmp1**2-
     4    (tmp2+tmp3)*tmp1pp/tmp1**2+
     5    2.0*(tmp2+tmp3)*(tmp1p**2)/tmp1**3)
c
c
c
c           bubble d function
cbub          scrdu = (1.0-u_nuc)*dmu**ovr3
cbub          scrd = scrdu/u_nuc
cbub          scrd2 = dmu**ovr3
cbub          scrdup = -1.0*dmu**ovr3-
cbub     1    ovr3*(1.0-u_nuc)*dmup*dmu**m2ovr3
cbub          scrdpp = ovr23*dmup*dmu**m2ovr3-ovr29*(1.0-u_nuc)*
cbub     1    dmu**m5ovr3*dmup**2+ovr3*(1.0-u_nuc)*dmu**m2ovr3*dmupp
c
c         
c           nuclei d function
cnuc          scrdu = u_nuc*du**ovr3
cnuc          scrd = du**ovr3
cnuc          scrd2 = scrdu/(1.0-u_nuc)
cnuc          scrdup = du**ovr3+ovr3*u_nuc*dup*du**m2ovr3
cnuc          scrdpp = ovr23*dup*du**m2ovr3-ovr29*u_nuc*
cnuc     1    (du**m5ovr3)*(dup**2)+ovr3*u_nuc*(du**m2ovr3)*dupp
c
c
c
          zeta_0 = csscal*6.035204*(sig_s*(16.0+q))**ovr23
c
c                        surface energy coefficent
          zeta = zeta_0*(h*sigma*x*nsubi)**ovr23
c
c                        derivative of zeta w.r.t. x
          dzdt = ovr23*zeta*hprim/h
c
c                        derivative of zeta w.r.t. x
          dzdx = ovr23*zeta*(dhdx/h+sigsgp+1.0/x)
c
c                        derivative of zeta w.r.t. nsubi
          dzdni = ovr23*zeta/nsubi
c
c
c
c                        nuclear radius
          rsubn = 9.0*h*sigma*sig_s*(16.0d0+q)*u_nuc*(1.0-u_nuc)/
     1    (2.0*zeta*scrdu)
c
c                        nuclear volume
          vsubn = 4.0*pi*(rsubn**3)/3.0
c
c                        atomic number
          a = nsubi*vsubn
c
c                        now calc surface, coulomb free energies
c
          fsubsc = zeta*scrdu/brydns
          fsubs = ovr23*zeta*scrdu/brydns
          fsubc = ovr3*zeta*scrdu/brydns
c
c
c
c                   translational chemical potential
          musubt = trscal*
     1        t*dlog((1.0-u_nuc)*(u_nuc*nsubi)/(nq*azero**2.5))
c
c                   derivative of trans. chem. potential w.r.t. t
          dmutdt = trscal*(musubt/t-1.5)
c
c                   translational free energy per baryon
          ftrans = trscal*h*(musubt-t)/azero
c
c                            if t is above the critical temperature
        else
          a = 0.0
          rsubn = 0.0
          vsubn = 0.0
          fsubs = 0.0
          fsubc = 0.0
          ftrans = 0.0
        endif
c                            calc ratio of nsubi to nsubs
        nratio = nsubi/nsubs
c
c
c20        vni = 2.0*aa*nsubi+4.0*bb*x*nsubi+cc*(1.0+dd)*nsubi**dd
        vni = pvn(x*nsubi,(1.0-x)*nsubi)
c
c20        vpi = 2.0*aa*nsubi+4.0*bb*(1.0-x)*nsubi+
c20     1    cc*(1.0+dd)*nsubi**dd+deltam
        vpi = pvp(x*nsubi,(1.0-x)*nsubi)
c
c---------------------------------------------------
c
        zni = 2.0*(pi**2)*nsubi*(1.0-x)/mq
c
        zpi = 2.0*(pi**2)*nsubi*x/mq
c
        eta_ni = finv12(zni)
c
        eta_pi = finv12(zpi)
c
        mun_i = t*eta_ni+vni
c
        mup_i = t*eta_pi+vpi
c
        f32_ni = f_3_2(eta_ni)
c
        f32_pi = f_3_2(eta_pi)
c
c20        psubi = lq*(f32_ni+f32_pi)+
c20     1    (nsubi**2)*(aa+4.0*bb*x*(1.0-x))+dd*cc*nsubi**(1.0+dd)
        psubi = lq*(f32_ni+f32_pi)+pv_pr(x*nsubi,(1.0-x)*nsubi)
c
c
        bn = ovr23*zeta*scrd*(sigsgp+hx+1.5*scrdux/scrdu)*x/nsubi-
     1  trscal*(1.0-u_nuc)*(musubt*(h-x*dhdx)/azero+x*dhdx*t/azero)
c
        bp = -ovr23*zeta*scrd*
     1 ((sigsgp+hx+1.5*scrdux/scrdu)*compx+1.0/x)/nsubi-
     1 trscal*(1.0-u_nuc)*
     2 (musubt*(h+dhdx*compx)/azero-dhdx*t*compx/azero)
c
        bsubp = zeta*scrdup-ovr23*zeta*scrd-
     1        trscal*u_nuc*nsubi*h*musubt/azero
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c
cc        gpi = 2.0*fhalfi(eta_pi)
cc        gpi = 2.0*fhalfi(2.0*(pi**2)*x*nsubi/mq)
cc        gni = 2.0*fhalfi(eta_ni)
cc        gni = 2.0*fhalfi(2.0*(pi**2)*(1.0-x)*nsubi/mq)
c
cc        gpo = 2.0*fhalfo(eta_po)
cc        gno = 2.0*fhalfo(eta_no)
c
c
        gpo = 2.0*fhalf(eta_po)
        gno = 2.0*fhalf(eta_no)
        gpi = 2.0*fhalf(eta_pi)
        gni = 2.0*fhalf(eta_ni)
c
c                  derivatives of inside potentials
c
c20        dvpidp = 2.0*aa+dd*(1.0+dd)*cc*(nsubi**(dd-1.0))
        dvpidp = dpvpdp(x*nsubi,(1.0-x)*nsubi)
c20        dvpidn = 2.0*aa+4.0*bb+dd*(1.0+dd)*cc*(nsubi**(dd-1.0))
        dvpidn = dpvpdn(x*nsubi,(1.0-x)*nsubi)
c20        dvnidp = dvpidn
        dvnidp = dpvndp(x*nsubi,(1.0-x)*nsubi)
c20        dvnidn = dvpidp
        dvnidn = dpvndn(x*nsubi,(1.0-x)*nsubi)
c
c                  derivatives of outside potentials
c
c20        dvpodp = eiflag*(2.0*aa+dd*(1.0+dd)*cc*(nout**(dd-1.0)) )
        dvpodp = eiflag*dpvpdp(npout,nnout)
c20        dvpodn = eiflag*(2.0*aa+4.0*bb+dd*(1.0+dd)*cc*(nout**(dd-1.0)))
        dvpodn = eiflag*dpvpdn(npout,nnout)
c20        dvnodp = dvpodn
        dvnodp = eiflag*dpvndp(npout,nnout)
c20        dvnodn = dvpodp
        dvnodn = eiflag*dpvndn(npout,nnout)
c
c                  derivatives of inside k.e. densities
c
        msscon = 3.0*massn/((hbar*c)**2)
        dtpidp = msscon*t*gpi
        dtpidn = 0.0
        dtnidp = 0.0
        dtnidn = msscon*t*gni
c
c                  derivatives of outside k.e. densities
c
        dtpodp = msscon*t*gpo
        dtpodn = 0.0
        dtnodp = 0.0
        dtnodn = msscon*t*gno
c
c
c                  derivatives of inside chem. potentials
c
        dmpidp = t*gpi/(x*nsubi)+dvpidp
        dmpidn = dvpidn
        dmnidp = dvnidp
        dmnidn = t*gni/((1.0-x)*nsubi)+dvnidn
c
c                  derivatives of outside chem. potentials
c
        dmpodp = t+dvpodp*npout/gpo
        dmpodn = dvpodn*nnout/gno
        dmnodp = dvnodp*npout/gpo
        dmnodn = t+dvnodn*nnout/gno
c
c                  derivatives of inside pressure
c
        dpidp = x*nsubi*dmpidp+(1.0-x)*nsubi*dmnidp
        dpidn = x*nsubi*dmpidn+(1.0-x)*nsubi*dmnidn
c
c                  derivatives of outside pressure
c
        dpodp = npout*dmpodp+nnout*dmnodp
        dpodn = npout*dmpodn+nnout*dmnodn
c
c                  derivatives of alpha pressure
c
        dpadp = alfdns*
     1  ( (2.0-npout*v_alfa)*dmpodp+(2.0-nnout*v_alfa)*dmnodp )
        dpadn = alfdns*
     1  ( (2.0-npout*v_alfa)*dmpodn+(2.0-nnout*v_alfa)*dmnodn )
c
c
        n1 = nsubi-exalfa*(nnout+npout)-4.0*alfdns
        n2 = nsubi*x-exalfa*npout-2.0*alfdns
c
c                  derivatives of u
c
        dudpo = -exclu*(exalfa*npout/gpo+
     1           (4.0-nout*v_alfa)*dpadp/t)/n1
        dudno = -exclu*(exalfa*nnout/gno+
     1           (4.0-nout*v_alfa)*dpadn/t)/n1
        dudni = -u_nuc/n1
c
c                  derivatives of x
c
        dxdpo = -(n2*dudpo+exclu*(exalfa*npout/gpo+
     1           (2.0-npout*v_alfa)*dpadp/t))/(u_nuc*nsubi)
        dxdno = -(n2*dudno+exclu*(2.0-npout*v_alfa)*dpadn/t)/
     1           (u_nuc*nsubi)
        dxdni = (n2-x*n1)/(nsubi*n1)
c
c                  derivatives of b's w.r.t. nsubi
c
        db1dni = trscal*( -u_nuc*h*(musubt+t)/azero )+
     1      ovr23*zeta*(scrdup-ovr23*scrd)/nsubi
c
        tmp4 = (sigsgp+dhdx/h+1.5*scrdux/scrdu)*(x-1.0)-1.0/x
c
        db2dni = -2.0*zeta*scrd*tmp4/(9.0*nsubi**2)-
     1  trscal*( (compu*t/(azero*nsubi))*(h+compx*dhdx) )
c
        tmp4 = (sigsgp+dhdx/h+1.5*scrdux/scrdu)
        db3dni = -2.0*zeta*scrd*x*tmp4/(9.0*nsubi**2)-
     1          trscal*( ((compu*t)/(azero*nsubi))*(h-x*dhdx) )
c
c
c
c                  derivatives of b's w.r.t. x
c
        db1dx = ovr23*zeta*(scrdup-ovr23*scrd)*(sigsgp+dhdx/h+1.0/x)+
     1  ovr23*zeta*(scrdpx-ovr23*scrdx)-
     2  trscal*( u_nuc*nsubi*dhdx*musubt/azero )
c
c
c
c
c
        tmp4 = (sigsgp+dhdx/h+1.5*scrdux/scrdu)*(x-1.0)-1.0/x
c
        tmp5 = sigsgp+dhdx/h+1.5*scrdux/scrdu+(x**(-2))+(x-1.0)*
     1  (sigsg2-(sigsgp**2)-(dhdx/h)**2+dhdxx/h+1.5*scrdxx/scrdu-
     2  1.5*(scrdux/scrdu)**2)
c
c
        db2dx = ovr23*(zeta*scrdux+scrdu*dzdx)*tmp4/(u_nuc*nsubi)+
     1      ovr23*zeta*scrd*tmp5/nsubi-trscal*( 
     2      compu*(dhdx*musubt+(dhdxx*(1.0-x)-dhdx)*(musubt-t))/azero)
c
c
c
c
        tmp4 = (sigsgp+dhdx/h+1.5*scrdux/scrdu)*x
c
        tmp5 = sigsgp+dhdx/h+1.5*scrdux/scrdu+x*
     1         (sigsg2-(sigsgp**2)-(dhdx/h)**2+dhdxx/h+
     2       1.5*scrdxx/scrdu-1.5*(scrdux/scrdu)**2)
c
        db3dx = ovr23*(zeta*scrdux+scrdu*dzdx)*tmp4/(u_nuc*nsubi)+
     1      ovr23*zeta*scrd*tmp5/nsubi-
     2      trscal*( compu*(dhdx*t-x*dhdxx*(musubt-t))/azero )
c
c
c
c                  derivatives of b's w.r.t. u_nuc
c
        db1du = zeta*(scrdpp-ovr23*scrdup/u_nuc+ovr23*scrd/u_nuc)-
     1  trscal*( nsubi*h*(musubt+t*(1.0-2.0*u_nuc)/(1.0-u_nuc))/azero )
c
        tmp4 = (sigsgp+dhdx/h+1.5*scrdux/scrdu)*(x-1.0)-1.0/x
        tmp5 = (x-1.0)*1.5*(scrdpx/scrdu-scrdux*scrdup/scrdu**2)
        db2du = (ovr23*zeta*scrd/nsubi)*tmp4*(scrdup/scrdu-1.0/u_nuc)+
     1    ovr23*zeta*scrdu*tmp5/(u_nuc*nsubi)+
     1    trscal*( (h*musubt+dhdx*compx*(musubt-t))/azero-
     2    (t*(1.0-2.0*u_nuc)/u_nuc)*(h+dhdx*compx)/azero )
c
        tmp4 = (sigsgp+dhdx/h+1.5*scrdux/scrdu)*x
        tmp5 = x*1.5*(scrdpx/scrdu-scrdup*scrdux/scrdu**2)
        db3du = ovr23*zeta*scrd*tmp4*(u_nuc*scrdup/scrdu-1.0)/
     1 (u_nuc*nsubi)+ovr23*zeta*scrdu*tmp5/(u_nuc*nsubi)+
     2  trscal*( (h*musubt-x*dhdx*(musubt-t))/azero-
     3 t*(1.0-2.0*u_nuc)*(h-x*dhdx)/(azero*u_nuc) )
c
c
c                      a1 derivatives
c
        da1id1 = x*dpidp+(1.0-x)*dpidn+nsubi*(dpidp-dpidn)*dxdni
        da1id2 = nsubi*(dpidp-dpidn)*dxdpo
        da1id3 = nsubi*(dpidp-dpidn)*dxdno
c
        da1od1 = 0.0
        da1od2 = dpodp+dpadp
        da1od3 = dpodn+dpadn
c
        db1d1 = db1dni+db1dx*dxdni+db1du*dudni
        db1d2 = db1dx*dxdpo+db1du*dudpo
        db1d3 = db1dx*dxdno+db1du*dudno
c
        da1d1 = da1id1-db1d1-da1od1
        da1d2 = da1id2-db1d2-da1od2
        da1d3 = da1id3-db1d3-da1od3
c
c                      a3 derivatives
c
        da3id1 = x*dmnidp+(1.0-x)*dmnidn+nsubi*(dmnidp-dmnidn)*dxdni
        da3id2 = nsubi*(dmnidp-dmnidn)*dxdpo
        da3id3 = nsubi*(dmnidp-dmnidn)*dxdno
c
        da3od1 = 0.0
        da3od2 = dmnodp
        da3od3 = dmnodn
c
        db3d1 = db3dni+db3dx*dxdni+db3du*dudni
        db3d2 = db3dx*dxdpo+db3du*dudpo
        db3d3 = db3dx*dxdno+db3du*dudno
c
        da3d1 = da3id1-db3d1-da3od1
        da3d2 = da3id2-db3d2-da3od2
        da3d3 = da3id3-db3d3-da3od3
c
c                      a2 derivatives
c
        da2id1 = x*dmpidp+(1.0-x)*dmpidn+nsubi*(dmpidp-dmpidn)*dxdni
        da2id2 = nsubi*(dmpidp-dmpidn)*dxdpo
        da2id3 = nsubi*(dmpidp-dmpidn)*dxdno
c
        da2od1 = 0.0
        da2od2 = dmpodp
        da2od3 = dmpodn
c
        db2d1 = db2dni+db2dx*dxdni+db2du*dudni
        db2d2 = db2dx*dxdpo+db2du*dudpo
        db2d3 = db2dx*dxdno+db2du*dudno
c
        da2d1 = da2id1-db2d1-da2od1
        da2d2 = da2id2-db2d2-da2od2
        da2d3 = da2id3-db2d3-da2od3
c
c
c                      eta derivatives
c
        dndetn = nnout/gno
        dpdetp = npout/gpo
c
        da1dn = da1d1
        da1etp = da1d2
        da1etn = da1d3
c
        da2dn = da2d1
        da2etp = da2d2
        da2etn = da2d3
c
        da3dn = da3d1
        da3etp = da3d2
        da3etn = da3d3
c
c
c
c
        a1 = psubi-bsubp-bprout-bpralf
        a2 = mup_i-bp-mup_o
        a3 = mun_i-bn-mun_o
c
c
c                          unset the "new" flag
        newflg = 0
c
        determ = da1dn*(da2etp*da3etn-da2etn*da3etp)-
     1           da1etp*(da2dn*da3etn-da2etn*da3dn)+
     2           da1etn*(da2dn*da3etp-da2etp*da3dn)
c
        dnsubi = -1.0*(a1*(da2etp*da3etn-da2etn*da3etp)+
     1           a2*(da3etp*da1etn-da1etp*da3etn)+
     2           a3*(da1etp*da2etn-da1etn*da2etp))/determ
c
c
        detap = -1.0*(a1*(da2etn*da3dn-da2dn*da3etn)+
     1          a2*(da1dn*da3etn-da1etn*da3dn)+
     2          a3*(da1etn*da2dn-da1dn*da2etn))/determ
c
c
        detan = -1.0*(a1*(da2dn*da3etp-da2etp*da3dn)+
     1          a2*(da1etp*da3dn-da1dn*da3etp)+
     2          a3*(da1dn*da2etp-da1etp*da2dn))/determ
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c JCM: lieb insert next 4 lines
c        if (DNSUBI.lt.1.d-99) then
c           SSFLAG = 0
c           goto 999
c        endif
c
c
c
c                        check the step size in nsubi
        if(abs(dnsubi/nsubi).gt.0.04) then
          dnsubi = 0.04*dnsubi*nsubi/abs(dnsubi)
        endif
 26     continue
        nsubin = nsubi+dnsubi
c JCM:
c        write(*,*) maxit,i
        if((nsubin.lt.dmax1(4.5d-2,brydns)).or.(nsubin.gt.0.25)) then
          dnsubi = 0.5*dnsubi
          goto 26
        endif
c
c                        check the step size in eta_po
        if(abs(detap).gt.4.0) then
          detap = 4.0*detap/abs(detap)
        endif
 27     continue
        netap = eta_po+detap
c JCM: next 4 lines lieb insert
c Causes no solution to be in larger domain
c        if (DETAP.lt.1.d-99) then
c           SSFLAG = 0.
c           goto 999
c        endif

        if((netap.lt.-5000.0).or.(netap.gt.etamax)) then
          detap = 0.5*detap
          goto 27
        endif
c
c                        check the step size in eta_no
        if(abs(detan).gt.4.0) then
          detan = 4.0*detan/abs(detan)
        endif
 28     continue
        netan = eta_no+detan
c JCM: next 4 lines, lieb insert
c Causes solution to not exist in much larger domain
c        if (DETAN.lt.1.d-99) then
c           SSFLAG = 0
c           goto 999
c        endif

        if((netan.lt.-5000.0).or.(netan.gt.etamax)) then
          detan = 0.5*detan
          goto 28
        endif
c
c
c                        update the variables
ccc        if(i.lt.30) write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc
 1205   format(i3,1p9e21.14)
c
        nsubi = nsubi+dnsubi
        eta_po = eta_po+detap
        eta_no = eta_no+detan
c
c
c
c                        if the required tolarences have been met
c                        break out of the loop
        if((abs(dnsubi).lt.nsiacc).and.(abs(detap).lt.prtacc)
     1    .and.(abs(detan).lt.nutacc) ) then

c          write(6,*) 'converged going  to 40'
c          write(6,*) abs(dnsubi),nsiacc
c          write(6,*) abs(detap),prtacc
c          write(6,*) abs(detan),nutacc

          goto 40
        else
      if(dbflag.eq.1) then
        write(*,2000) '2',i,nsubi,eta_po,eta_no,dnsubi
      endif
          goto 30
        endif
c
c
 29     continue
        if(newflg.ne.1) then
          cflag = cflag+1
          dnsubi = 0.5*dnsubi
          nsubi = nsubi-dnsubi
          detap = 0.5*detap
          eta_po = eta_po-detap
          detan = 0.5*detan
          eta_no = eta_no-detan
          if(dbflag.eq.1) then
            write(*,2000) '3',i,nsubi,eta_po,eta_no,dnsubi
          endif
          goto 30
        else
          nsubi = nsubs
cc          eta_po = eta_po-0.5/t
cc          eta_no = eta_no-0.5/t
          eta_po = eta_po-2.0/t
          eta_no = eta_no-2.0/t
        endif
c
c
      if(dbflag.eq.1) then
        write(*,2000) '4',i,nsubi,eta_po,eta_no,dnsubi
      endif
 2000   format(t2,a,1x,i3,1x,f8.5,3(1x,g13.5))
c
c
 30   continue
c
c            if scheme 1 has failed try scheme 2
      if(schflg.eq.0) then
        schflg = 1
        goto 5
      endif
c
c
      ssflag = 0
      goto 999
c
c                    branch label to break out of do 30 iteration
 40   continue
c
c
c                    the following logic determines whether this was
c                    the correct scheme to use, and if not then which
c                    one should be used
c
      if(ftflag.ne.0) then
        ssflag = 4
        goto 999
      endif
c
c                    if calculated critical temperature is less than t,
c                    then switch to the scheme with no nuclei
      if(t.ge.tsubc) then
c                    set flag to indicate failure
        ssflag = 0
        goto 999
      endif
c
c
c                    if fraction of nuclei present is zero and no switch
c                    has been made then switch to the no nuclei scheme
      if(u_nuc.le.0.0) then
c                    set flag to indicate failure
        ssflag = 0
        goto 999
      elseif(u_nuc.gt.1.0) then
c                    set flag to indicate failure
        ssflag = 0
        goto 999
      else
c                    set flag to indicate success
        ssflag = 1
      endif
c
c
c
c                    if eqns aren't really zeroed then fail
c
c
      if( (abs(a1).gt.1.0d-5).or.(abs(a2).gt.1.0d-5).or.
     1    (abs(a3).gt.1.0d-5) ) then
        ssflag = 0
        write(*,*) ' nuceos: false convg; a = ',a1,a2,a3
        goto 999
      endif
c
c
c
c
      if(nsubi.lt.0.05) then
        write(*,*) 'nuceos:: <<warning>> nsubi getting close to lb'
      endif


c..      write(6,*) ' success: maxit = ',i

c
c
c
      zni = 2.0*(pi**2)*nsubi*(1.0-x)/mq
c
      zpi = 2.0*(pi**2)*nsubi*x/mq
c
      eta_ni = finv12(zni)
c
      eta_pi = finv12(zpi)
c
      mun_i = t*eta_ni+vni
c
      mup_i = t*eta_pi+vpi
c
      f32_ni = f_3_2(eta_ni)
c
      f32_pi = f_3_2(eta_pi)
c
      exclu = 1.0-u_nuc
      exalfa = 1.0-alfdns*v_alfa
c
c
c
c        JCM: These are particle fractions per unit baryon mass (i.e. not by number, which is why alphas have 4.0)
c                    calculate particle fractions
c
      xalfa = 4.0*exclu*alfdns/brydns
      xnut = nnout*exclu*exalfa/brydns
      xprot = npout*exclu*exalfa/brydns
      xh = 1.0-xprot-xnut-xalfa
      xhchk = u_nuc*nsubi/brydns
c
      if((xh.lt.heavct).or.(xhchk.lt.heavct)) then
c                    set flag to indicate switch is being made
        ssflag = 0
cc        write(*,*) ' xh,xhchk = ',xh,xhchk
        goto 999
      endif
c
      if((xalfa.lt.0.0).or.(xh.lt.0.0).or.
     1   (xnut.lt.0.0).or.(xprot.lt.0.0)) then   
        ssflag = 0
        write(*,*) ' xs hnpa = ',xh,xnut,xprot,xalfa
        goto 999
      endif
c
c
c
c
c                    baryons
c
c
      muprot = mup_o
      mun = mun_o
      muhat = mun-muprot
c JCM: dimensionless chemical potentials for protons and neutrons
      etapls=muprot/t
      etanls=mun/t
c
c 
      if(abs((xh-xhchk)/xhchk).gt.1.0d-4) then
        ssflag = 0
        goto 999
ccc        write(*,*) ' inconsistencey in xh at',t,brydns,ye,xh,xhchk
      endif
c   
      nucdns = brydns*xh
c
      tau_po = kq*f32_po
      tau_pi = kq*f32_pi
c
      tau_no = kq*f32_no
      tau_ni = kq*f32_ni
c
      if(nout.gt.0.0) xout = npout/nout
c
c
c                    calculate internal energy of outside nucleons,
c                    alpha particles, and nuclei (per baryon)
c
c20      buout = (exclu*exalfa/brydns)*( uq*(tau_po+tau_no)+eiflag*
c20     1    ( (nout**2)*aa+4.0*bb*npout*nnout+
c20     2    cc*nout**(1.0+dd)+npout*deltam) )
      buout = (exclu*exalfa/brydns)*( 
     1    uq*(tau_po+tau_no)+eiflag*pv_e(npout,nnout) )
c
c20      bunuc = xh*( ( uq*(tau_pi+tau_ni)+(nsubi**2)*
c20     1 (aa+4.0*bb*x*(1.0-x))+cc*nsubi**(1.0+dd)+x*nsubi*deltam )/
c20     2 nsubi)+fsubsc*(1.0-t*(scrdut/scrdu+ovr23*hprim/h))+
c20     3 trscal*
c20     4 (1.0-u_nuc)*xh*(ftrans*(1.0-t*hprim/h)-h*(musubt-2.5*t)/azero)
      bunuc = xh*( (uq*(tau_pi+tau_ni)+pv_e(x*nsubi,(1.0-x)*nsubi))/
     2 nsubi)+fsubsc*(1.0-t*(scrdut/scrdu+ovr23*hprim/h))+
     3 trscal*
     4 (1.0-u_nuc)*xh*(ftrans*(1.0-t*hprim/h)-h*(musubt-2.5*t)/azero)
c

c
      bualfa = 0.25*xalfa*(1.5*t-balpha)
c
      bu = buout+bualfa+bunuc
c
c
      bsout = (exclu*exalfa/brydns)*( (5.0*uq/(3.0*t))*(tau_no+tau_po)-
     1 nnout*eta_no-npout*eta_po )
c
c
c                    calculate entropy of alpha particles (per baryon)
      bsalfa = -0.25*xalfa*(mualfa/t-2.5)
c
c
      bsnuc = xh*( (5.0*uq/(3.0*t))*(tau_ni+tau_pi)-
     1 nsubi*(1.0-x)*eta_ni-nsubi*x*eta_pi )/nsubi-
     2 fsubsc*(scrdut/scrdu+ovr23*hprim/h)-
     3 xh*trscal*(1.0-u_nuc)*
     4 ((ftrans*hprim/h)+h*(musubt/t-2.5)/azero)
c
c                    calculate total baryon entropy (per baryon)
      bs = bsout+bsnuc+bsalfa
c
c                    calculate free energy of outside nucleons (per baryon)
      bfout = buout-t*bsout
c
c                    calculate free energy of alpha particles (per baryon)
      bfalfa = bualfa-t*bsalfa
c
c                    calculate free energy of nuclei (per baryon)
      bfnuc = bunuc-t*bsnuc
c
c                    calculate total baryon free energy (per baryon)
      bftot = bfout+bfnuc+bfalfa
c
c                    calculate pressure due to nuclei
      bprnuc = -zeta*(scrdu-u_nuc*scrdup)+
     1 trscal*u_nuc*nsubi*h*((1.0-u_nuc)*t-u_nuc*musubt)/azero
c
c
c                    calculate total baryon pressure
      bpress = bprout+bpralf+bprnuc
c
c
c                    leptons & photons
c

c JCM:
c      write(*,*) 'before anyelectron'
      call any_electron(t,ye,brydns) 

c      call el_eos(t,ye,brydns)



c
c
c
c                    total pressure and eng/ent per baryon
c
      fbary = bftot + fsube
      pbary = bpress+epress
      mubary = ye*muprot+(1.0-ye)*mun
      mu_mat = ye*(muprot+musube)+(1.0-ye)*mun


c..filter

      bftot  = bftot*ionadd
      bu     = bu*ionadd
      bs     = bs*ionadd
      bpress = bpress*ionadd

      fsube  = fsube*eleadd
      eu     = eu*eleadd
      es     = es*eleadd
      epress = epress*eleadd

      pf     = pf*radadd
      pu     = pu*radadd
      ps     = ps*radadd
      ppress = ppress*radadd


c..totals
      ftot = bftot  + fsube  + pf
      utot = bu     + eu     + pu
      stot = bs     + es     + ps
      ptot = bpress + epress + ppress



c      write(6,801) ' ftot ',
c     1 ftot, bftot, 100.*(ftot-bftot)/ftot, 
c     1       fsube, 100.*(ftot-fsube)/ftot,
c     1       pf,    100.*(ftot-pf)/ftot

c      write(6,801) ' utot ',
c     1 utot, bu, 100.*(utot-bu)/utot, 
c     1       eu, 100.*(utot-eu)/utot,
c     1       pu, 100.*(utot-pu)/utot

c      write(6,801) ' stot ',
c     1 stot, bs, 100.*(stot-bs)/stot, 
c     1       es, 100.*(stot-es)/stot,
c     1       ps, 100.*(stot-ps)/stot

c      write(6,801) ' ptot ',
c     1 ptot, bpress, 100.*(ptot-bpress)/ptot, 
c     1       epress, 100.*(ptot-epress)/ptot,
c     1       ppress, 100.*(ptot-ppress)/ptot
c801   format(a6,1pe8.1,3(1pe8.1,0pf7.1))



c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c                derivatives of thermodynamic variables
c-----------------------------------------------------------------------
c
c                 ------------------------------------
c                 !      derivatives of exterior     !
c                 !      quantities                  !
c                 !      (w.r.t. temp. and eta's)    !
c                 !                                  !
c                 ------------------------------------
c
c
c                  derivatives of exterior potentials
c                  w.r.t. particle densities
c20      dvpodp = eiflag*(2.0*aa+dd*(1.0+dd)*cc*(nout**(dd-1.0)) )
      dvpodp = eiflag*dpvpdp(npout,nnout)
c20      dvpodn = eiflag*(2.0*aa+4.0*bb+dd*(1.0+dd)*cc*(nout**(dd-1.0)))
      dvpodn = eiflag*dpvpdn(npout,nnout)
c20      dvnodp = dvpodn
      dvnodp = eiflag*dpvndp(npout,nnout)
c20      dvnodn = dvpodp
      dvnodn = eiflag*dpvndn(npout,nnout)
c
c
c                  derviatives of exterior chem. pot. w.r.t. eta's
c                  (at fixed t)
      dmpdep = t+dvpodp*npout/gpo
      dmpden = dvpodn*nnout/gno
      dmndep = dvnodp*npout/gpo
      dmnden = t+dvnodn*nnout/gno
c
c                  derivatives of pressure potential w.r.t.
c                  particle densities
c20      dv_dpo = eiflag*
c20     1    (2.0*aa*nout+4.0*bb*nnout+cc*dd*(1.0+dd)*(nout**dd) )
      dv_dpo = eiflag*dpvrdp(npout,nnout)
c20      dv_dno = eiflag*
c20     1    (2.0*aa*nout+4.0*bb*npout+cc*dd*(1.0+dd)*(nout**dd) )
      dv_dno = eiflag*dpvrdn(npout,nnout)
c
c                  derivatives of pressure potential w.r.t. eta's
c                  (at fixed t)
      dv_dep = dv_dpo*npout/gpo
      dv_den = dv_dno*nnout/gno
c
c                  derivatives of outside pressure w.r.t. eta's
c                  (at fixed t)
      dpodep = npout*t+dv_dep
      dpoden = nnout*t+dv_den
c
c                  derivatives of alpha density w.r.t. eta's
c                  (at fixed t)
      dnadep = alfdns*(2.0*dmpdep+2.0*dmndep-v_alfa*dpodep)/t
      dnaden = alfdns*(2.0*dmpden+2.0*dmnden-v_alfa*dpoden)/t
c
c                  derivatives of alpha pressure w.r.t. eta's
c                  (at fixed t)
      dpadep = t*dnadep
      dpaden = t*dnaden
c
c                  derivatives of particle densities w.r.t. t
c                  (at fixed eta's)
      dnpodt = 1.5*npout/t
      dnnodt = 1.5*nnout/t
c
c                  derivatives of exterior chem. pot. w.r.t. t
c                  (at fixed eta's)
      dmpodt = eta_po+dvpodp*dnpodt+dvpodn*dnnodt
      dmnodt = eta_no+dvnodp*dnpodt+dvnodn*dnnodt
c
c                  derivative of pressure potential w.r.t. t
c                  (at fixed eta's)
      dv_dt = dv_dpo*dnpodt+dv_dno*dnnodt
c
c                  derivative of outside pressure w.r.t. t
c                  (at fixed eta's)
      dpodt = ovr23*uq*2.5*(tau_po+tau_no)/t+dv_dt
c
c                  derivative of alpha chem. pot. w.r.t. t
c                  (at fixed eta's)
      dmuadt = 2.0*dmpodt+2.0*dmnodt-v_alfa*dpodt
c
c                  derivative of alpha particle density w.r.t. t
c                  (at fixed eta's)
      dnadt = 1.5*alfdns/t-alfdns*mualfa/(t**2)+alfdns*dmuadt/t
c
c                  derivative of alpha particle pressure w.r.t. t
c                  (at fixed eta's)
      dpadt = alfdns+t*dnadt
c
c
c                 ------------------------------------
c                 !      derivatives of interior     !
c                 !      quantities                  !
c                 !      (w.r.t. temp. and density)  !
c                 !                                  !
c                 ------------------------------------
c
c
c                   derivatives of kinetic energy densities w.r.t. t
c                   (holding the number densities (x & nsubi) fixed)
      dtpidt =2.5*tau_pi/t-2.25*x*nsubi*gpi/uq
      dtnidt =2.5*tau_ni/t-2.25*(1.0-x)*nsubi*gni/uq
c
c                   derivatives of pressures w.r.t. t
c                   (holding the number densities (x & nsubi) fixed)
      dpidt = ovr23*uq*(dtpidt+dtnidt)
c
c                   derivatives of interior chem. pot. w.r.t. t
c                   (holding the number densities (x & nsubi) fixed)
      dmpidt = eta_pi-1.5*gpi
      dmnidt = eta_ni-1.5*gni
c
c
c                  derivatives of inside potentials w.r.t.
c                  interior proton and neutron densities
c                  (at fixed t)
c20      dvpidp = 2.0*aa+dd*(1.0+dd)*cc*(nsubi**(dd-1.0))
      dvpidp = dpvpdp(x*nsubi,(1.0-x)*nsubi)
c20      dvpidn = 2.0*aa+4.0*bb+dd*(1.0+dd)*cc*(nsubi**(dd-1.0))
      dvpidn = dpvpdn(x*nsubi,(1.0-x)*nsubi)
c20      dvnidp = dvpidn
      dvnidp = dpvndp(x*nsubi,(1.0-x)*nsubi)
c20      dvnidn = dvpidp
      dvnidn = dpvndn(x*nsubi,(1.0-x)*nsubi)
c
c
c                   derivatives of interior chemical potentials
c                   w.r.t. interior neutron and proton densities
c                  (at fixed t)
      dmpidp = t*gpi/(x*nsubi)+dvpidp
      dmpidn = dvpidn
      dmnidp = dvnidp
      dmnidn = t*gni/((1.0-x)*nsubi)+dvnidn
c
c                   derivatives of interior pressure
c                   w.r.t. interior neutron and proton densities
c                  (at fixed t)
      dpidp = x*nsubi*dmpidp+(1.0-x)*nsubi*dmnidp
      dpidn = x*nsubi*dmpidn+(1.0-x)*nsubi*dmnidn
c
c
c
c
c                 ------------------------------------
c                 !      derivatives of "b" terms    !
c                 !      from the chemical and       !
c                 !      pressure equilibrium        !
c                 !      equations                   !
c                 !                                  !
c                 !      (w.r.t. temperature )       !
c                 !                                  !
c                 ------------------------------------
c
c
c             derivative of term from pressure equilibrium eqn.
c
      db1dt = ovr23*zeta*(scrdup-ovr23*scrd)*hprim/h+
     1    zeta*(scrdpt-ovr23*scrdt)-
     2    trscal*u_nuc*nsubi*(hprim*musubt+h*dmutdt)/azero
c
c
c             derivative of term from proton equilibrium eqn.
c
      tmp4 = (sigsgp+dhdx/h+1.5*scrdux/scrdu)*(x-1.0)-1.0/x
      tmp5 = dhdtdx/h-dhdx*hprim/h**2+
     1 1.5*scrdxt/scrdu-1.5*scrdux*scrdut/scrdu**2
c
      db2dt = ovr49*(zeta*scrd*hprim/(h*nsubi))*tmp4+
     1    ovr23*zeta*scrdt*tmp4/nsubi+
     2    ovr23*(zeta*scrd/nsubi)*(x-1.0)*tmp5-
     3    trscal*exclu*(dmutdt*(h+dhdx*(1.0-x))+musubt*
     4    (hprim+dhdtdx*(1.0-x))-dhdx*(1.0-x)-t*dhdx*(1.0-x))/azero
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c             derivative of term from neutron equilibrium eqn.
c
      tmp4 = sigsgp+dhdx/h+1.5*scrdux/scrdu
      tmp5 = dhdtdx/h-dhdx*hprim/h**2+
     1 1.5*scrdxt/scrdu-1.5*scrdux*scrdut/scrdu**2    
      db3dt = ovr49*(zeta*scrd*hprim/(h*nsubi))*x*tmp4+
     1        ovr23*(zeta*scrdt/nsubi)*x*tmp4+
     2        ovr23*(zeta*scrd/nsubi)*x*tmp5-
     3        trscal*exclu*(hprim*musubt+h*dmutdt-x*dhdtdx*(musubt-t)-
     4        x*dhdx*(dmutdt-1.0))/azero
c
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c
c                 ------------------------------------
c                 !      derivatives of constraint   !
c                 !      and equilibrium equations   !
c                 !      with respect to the five    !
c                 !      compositional variables     !
c                 !      (u,x,n_i,eta_po,eta_no)     !
c                 !      and the three independent   !
c                 !      variables                   !
c                 !      (baryon density, t, and ye) !
c                 !                                  !
c                 ------------------------------------
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                equation 1 (baryon conservation)
c
      dfdom(1,1) = nout*exalfa+4.0*alfdns-nsubi
c
      dfdom(1,2) = 0.0
c
      dfdom(1,3) = -u_nuc
c
      dfdom(1,4) = -exclu*exalfa*npout/gpo+
     1             v_alfa*dnadep*exclu*nout-4.0*exclu*dnadep
c
      dfdom(1,5) = -exclu*exalfa*nnout/gno+
     1             v_alfa*dnaden*exclu*nout-4.0*exclu*dnaden
c
c
c
      dfdl_1(1) = -1.0
c
      dfdl_2(1) = exclu*exalfa*(dnpodt+dnnodt)-exclu*v_alfa*nout*dnadt+
     1     4.0*exclu*dnadt
c            
      dfdl_3(1) = 0.0     
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                equation 2 (charge conservation)
c
      dfdom(2,1) = exalfa*npout+2.0*alfdns-x*nsubi
c
      dfdom(2,2) = -u_nuc*nsubi
c
      dfdom(2,3) = -x*u_nuc
c
      dfdom(2,4) = -exclu*exalfa*npout/gpo+
     1     v_alfa*exclu*npout*dnadep-2.0*exclu*dnadep
c
      dfdom(2,5) = v_alfa*exclu*npout*dnaden-2.0*exclu*dnaden             
c
c
c
      dfdl_1(2) = -1.0*ye
c
      dfdl_2(2) = exclu*exalfa*dnpodt-v_alfa*exclu*npout*dnadt+
     1     2.0*exclu*dnadt
c            
      dfdl_3(2) = -1.0*brydns
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                equation 3 (proton chemical equilibrium)
c
      dfdom(3,1) = -db2du
c
      dfdom(3,2) = nsubi*(dmpidp-dmpidn)-db2dx
c
      dfdom(3,3) = (1.0-x)*dmpidn+x*dmpidp-db2dni
c
      dfdom(3,4) = -dmpdep
c
      dfdom(3,5) = -dmpden
c
      dfdl_1(3) = 0.0
      dfdl_2(3) = -1.0*(dmpidt-dmpodt-db2dt)
      dfdl_3(3) = 0.0
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                equation 4 (neutron chemical equilibrium)
c
      dfdom(4,1) = -db3du
c
      dfdom(4,2) = nsubi*(dmnidp-dmnidn)-db3dx
c
      dfdom(4,3) = (1.0-x)*dmnidn+x*dmnidp-db3dni
c
      dfdom(4,4) = -dmndep
c
      dfdom(4,5) = -dmnden
c
      dfdl_1(4) = 0.0
      dfdl_2(4) = -1.0*(dmnidt-dmnodt-db3dt)
      dfdl_3(4) = 0.0
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                equation 5 (pressure equilibrium)
c
      dfdom(5,1) = -db1du
c
      dfdom(5,2) = nsubi*(dpidp-dpidn)-db1dx
c
      dfdom(5,3) = (1.0-x)*dpidn+x*dpidp-db1dni
      ncomp = dfdom(5,3)
c
      dfdom(5,4) = -dpodep-dpadep
c
      dfdom(5,5) = -dpoden-dpaden
c
      dfdl_1(5) = 0.0
      dfdl_2(5) = -1.0*(dpidt-dpodt-dpadt-db1dt)
      dfdl_3(5) = 0.0
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                    lu decompose the dfdom matrix
      call matlud(dfdom,dfdmlu,5,ipvt)
c
c                    now solve the lu decomposed linear system
c                    to get the density derivatives
      call mluslv(dfdmlu,result,dfdl_1,ipvt,5)
c
      du_dn = result(1)
      dx_dn = result(2)
      dni_dn = result(3)
      dep_dn = result(4)
      den_dn = result(5)
c
c
c                    now solve the lu decomposed linear system
c                    to get the temperature derivatives
      call mluslv(dfdmlu,result,dfdl_2,ipvt,5)
c
      du_dt = result(1)
      dx_dt = result(2)
      dni_dt = result(3)
      dep_dt = result(4)
      den_dt = result(5)
c
c                    now solve the lu decomposed linear system
c                    to get the ye derivatives
      call mluslv(dfdmlu,result,dfdl_3,ipvt,5)
c
      du_dy = result(1)
      dx_dy = result(2)
      dni_dy = result(3)
      dep_dy = result(4)
      den_dy = result(5)
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                 ------------------------------------
c                 !      derivatives of finite size  !
c                 !      terms in the internal       !
c                 !      energy and entropy          !
c                 !      densities w.r.t. to u,x,n_i !
c                 !      and t.  these are used in   !
c                 !      calculating the derivatives !
c                 !      w.r.t. the independant vars !
c                 !      (baryon density, t, and ye) !
c                 !                                  !
c                 ------------------------------------
c
c                        free energy surface & coulomb terms
c                                  (densities)
c
      f_sc = zeta*scrdu
c
      dfscdu = zeta*scrdup
c
      dfscdx = zeta*scrdux+scrdu*dzdx
c
      dfscdn = scrdu*dzdni
c
      dfscdt = zeta*scrdut+scrdu*dzdt
c
c
c                        free energy translational terms
c                                  (densities)
      ftr = u_nuc*exclu*nsubi*ftrans
c
      dftrdt = ftr*(hprim/h+1.0/t)-
     1    1.5*trscal*u_nuc*exclu*nsubi*h/azero
c
      dftrdx = ftr*dhdx/h
c
      dftrdu = ftr/u_nuc-ftr/exclu+
     1    trscal*nsubi*h*(1.0-2.0*u_nuc)/azero
c
      dftrdn = ftr/nsubi+trscal*u_nuc*exclu*h*t/azero
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
c                        internal energy surface & coulomb terms
c                                  (densities)
c
      tmp4 = 1.0-t*scrdut/scrdu-ovr23*t*hprim/h
c
      e_sc = f_sc*tmp4
c
      descdu = dfscdu*tmp4+
     1    f_sc*(t*scrdut*scrdup/scrdu**2-t*scrdpt/scrdu)
c
      descdx = dfscdx*tmp4+
     1    f_sc*(t*scrdut*scrdux/scrdu**2-t*scrdxt/scrdu+
     2    ovr23*t*hprim*dhdx/h**2-ovr23*t*dhdtdx/h)
c
      descdn = dfscdn*tmp4
c
      descdt = dfscdt*tmp4+f_sc*
     1   (t*(scrdut**2)/scrdu**2-scrdut/scrdu-t*scrdtt/scrdu+
     2    ovr23*t*(hprim**2)/h**2-ovr23*hprim/h-ovr23*t*hpprim/h)
c
c                        internal energy translational terms
c                                  (densities)
c
      tmp4 = 1.5*h*t/azero-t*hprim*(musubt-t)/azero
c
      e_tr = trscal*exclu*brydns*xh*tmp4
c
      detrdu = trscal*(nsubi*(1.0-2.0*u_nuc)*tmp4-
     1    nsubi*(t**2)*hprim*(1.0-2.0*u_nuc)/azero)
c
      detrdx = trscal*brydns*xh*exclu*
     1    (1.5*t*dhdx/azero-t*(musubt-t)*dhdtdx/azero)
c
      detrdn = trscal*(u_nuc*exclu*tmp4-
     1    brydns*xh*exclu*(t**2)*hprim/(nsubi*azero))
c
      detrdt = trscal*brydns*xh*exclu*
     1    (1.5*(h+t*hprim)/azero-(hprim+t*hpprim)*(musubt-t)/azero-
     2    t*hprim*(musubt/t-2.5)/azero )
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
c                        entropy surface & coulomb terms
c                                  (densities)
c
      s_sc = (e_sc-f_sc)/t
c
      dsscdu = (descdu-dfscdu)/t
c
      dsscdx = (descdx-dfscdx)/t
c
      dsscdn = (descdn-dfscdn)/t
c
      dsscdt = (descdt-dfscdt)/t-(e_sc-f_sc)/t**2
c
c                        entropy translational terms
c                                  (densities)
c
      tmp4 = musubt*(hprim+h/t)/azero-(t*hprim+2.5*h)/azero
c
      s_tr = -trscal*brydns*xh*exclu*tmp4
c
      dstrdu = -trscal*(nsubi*(1.0-2.0*u_nuc)*tmp4+
     1    nsubi*t*(1.0-2.0*u_nuc)*(hprim+h/t)/azero)
c
      dstrdx = -trscal*brydns*xh*exclu*
     1    (musubt*(dhdtdx+dhdx/t)/azero-
     2    (t*dhdtdx+2.5*dhdx)/azero)
c
      dstrdn = -trscal*
     1    (u_nuc*exclu*tmp4+u_nuc*exclu*t*(hprim+h/t)/azero)
c
      dstrdt = -(brydns*xh*exclu*((musubt/t-1.5)*(hprim+h/t)/azero+
     1    musubt*(hpprim+hprim/t-h/t**2)/azero-
     2    (3.5*hprim+t*hpprim)/azero ))*trscal
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                 -------------------------------------
c                 !      derivatives of interior bulk !
c                 !      terms in the internal        !
c                 !      energy and entropy           !
c                 !      densities w.r.t. to u,x,n_i  !
c                 !      and t.  these are used in    !
c                 !      calculating the derivatives  !
c                 !      w.r.t. the independant vars  !
c                 !      (baryon density, t, and ye)  !
c                 !                                   !
c                 -------------------------------------
c
c
c
      s_nuc =(ovr53*uq/t)*(tau_ni+tau_pi)-
     1    nsubi*((1.0-x)*eta_ni+x*eta_pi)
c
c20      e_nuc = uq*(tau_pi+tau_ni)+(nsubi**2)*(aa+4.0*bb*x*(1.0-x))+
c20     1    cc*nsubi**(1.0+dd)+x*nsubi*deltam
      e_nuc = uq*(tau_pi+tau_ni)+pv_e(x*nsubi,(1.0-x)*nsubi)
c
c
c                    interior particle densties
      npi = x*nsubi
      nni = (1.0-x)*nsubi
c
      dtpidt = 2.5*tau_pi/t-2.25*npi*gpi/uq
      dtnidt = 2.5*tau_ni/t-2.25*nni*gni/uq
c
c               derivative of interior entropy density w.r.t. t
      dsidt = uq*(dtpidt+dtnidt)/t
c
c               derivative of interior internal energy density w.r.t. t
      deidt = t*dsidt
c
c
c
c
c                    derivatives of eta's w.r.t. x and nsubi
      detpdx = gpi/x
      detndx = -gni/(1.0-x)
      detpdn = gpi/nsubi
      detndn = gni/nsubi
c
c                    derivatives of tau's w.r.t. x and nsubi
      dtpidx = 1.5*t*npi*detpdx/uq
      dtnidx = 1.5*t*nni*detndx/uq
      dtpdni = 1.5*t*npi*detpdn/uq
      dtndni = 1.5*t*nni*detndn/uq
c
c
c
c           derivative of interior entropy density w.r.t. x
      dsidx = ovr53*uq*(dtpidx+dtnidx)/t-nsubi*(eta_pi-eta_ni)-
     1    nsubi*((1.0-x)*detndx+x*detpdx)
c
c           derivative of interior internal energy density w.r.t. x
c20      deidx = uq*(dtpidx+dtnidx)+
c20     1    (nsubi**2)*4.0*bb*(1.0-2.0*x)+nsubi*deltam
      deidx = uq*(dtpidx+dtnidx)+dpvedx(nsubi,x)
c
c
c           derivative of interior entropy density w.r.t. nsubi
      dsidn = ovr53*uq*(dtpdni+dtndni)/t-((1.0-x)*eta_ni+x*eta_pi)-
     1    nsubi*((1.0-x)*detndn+x*detpdn)
c
c
c           derivative of interior internal energy density w.r.t. nsubi
c20      deidn = uq*(dtpdni+dtndni)+2.0*nsubi*(aa+4.0*bb*x*(1.0-x))+
c20     1    cc*(1.0+dd)*(nsubi**dd)+x*deltam
      deidn = uq*(dtpdni+dtndni)+dpvedn(nsubi,x)
c
c
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                 -------------------------------------
c                 !      derivatives of exterior bulk !
c                 !      nucleon internal energy &    !
c                 !      entropy densities and the    !
c                 !      chem. pot.  w.r.t. to eta_p, !
c                 !      ate_n & t. these are used in !
c                 !      calculating the derivatives  !
c                 !      w.r.t. the independant vars  !
c                 !      (baryon density, t, and ye)  !
c                 !                                   !
c                 -------------------------------------
c
c
      s_out =(ovr53*uq/t)*(tau_no+tau_po)-nnout*eta_no-npout*eta_po
c
c20      e_out = uq*(tau_po+tau_no)+eiflag*
c20     1((nout**2)*aa+4.0*bb*npout*nnout+cc*nout**(1.0+dd)+npout*deltam)
      e_out = uq*(tau_po+tau_no)+eiflag*pv_e(npout,nnout)
c
c                   derivative of exterior entropy density w.r.t. t
      dsodt =  ovr53*uq*(1.5*(tau_po+tau_no)/(t**2))-
     1     1.5*(npout*eta_po+nnout*eta_no)/t
c
      deodt = t*dsodt
c
c                    derivatives of exterior particle densities w.r.t.
c                    temperature (eta's fixed)
      dnpodt = 1.5*npout/t
      dnnodt = 1.5*nnout/t
c
      dmpodt = eta_po+dvpodp*dnpodt+dvpodn*dnnodt
      dmnodt = eta_no+dvnodp*dnpodt+dvnodn*dnnodt
c
c
      dnpdep = npout/gpo
      dnnden = nnout/gno
c
      dtpdep = 1.5*t*npout/uq
      dtnden = 1.5*t*nnout/uq
c
      dsodep = (ovr53*uq/t)*dtpdep-npout-eta_po*dnpdep
      dsoden = (ovr53*uq/t)*dtnden-nnout-eta_no*dnnden
c
c
c                    exterior particle potentials
c20      vnout = eiflag*(2.0*aa*nout+4.0*bb*npout+cc*(1.0+dd)*nout**dd )
      vnout = eiflag*pvn(npout,nnout)
c20      vpout = eiflag*
c20     1    (2.0*aa*nout+4.0*bb*nnout+cc*(1.0+dd)*nout**dd+deltam)
      vpout = eiflag*pvp(npout,nnout)
c
c
      deodep = uq*dtpdep+vpout*dnpdep
      deoden = uq*dtnden+vnout*dnnden
c
      dmpdep = t+dvpodp*npout/gpo
      dmpden = dvpodn*nnout/gno
      dmndep = dvnodp*npout/gpo
      dmnden = t+dvnodn*nnout/gno
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                 -------------------------------------
c                 !      derivatives of alpha         !
c                 !      particle internal energy &   !
c                 !      entropy densities and the    !
c                 !      chem. pot.  w.r.t. to eta_p, !
c                 !      ate_n & t. these are used in !
c                 !      calculating the derivatives  !
c                 !      w.r.t. the independant vars  !
c                 !      (baryon density, t, and ye)  !
c                 !                                   !
c                 -------------------------------------
c
c
c
c
      s_alfa = alfdns*(2.5-mualfa/t)
c
c
      e_alfa = alfdns*(1.5*t-balpha)
c
c                  derivative of pressure potential w.r.t. t
      dv_dt = dv_dpo*dnpodt+dv_dno*dnnodt
c
c                  derivative of outside pressure w.r.t. t
      dpodt = ovr23*uq*2.5*(tau_po+tau_no)/t+dv_dt
c
c
      dmuadt = 2.0*dmpodt+2.0*dmnodt-v_alfa*dpodt
c
c                  derivative of alpha particle density w.r.t. t
      dnadt = 1.5*alfdns/t-alfdns*mualfa/(t**2)+alfdns*dmuadt/t
c
c
      dsadt = dnadt*(2.5-mualfa/t)-alfdns*dmuadt/t+alfdns*mualfa/t**2
c
      deadt = dnadt*(1.5*t-balpha)+1.5*alfdns
c
c
      dv_dep = dv_dpo*npout/gpo
      dv_den = dv_dno*nnout/gno
c
      dpodep = ovr23*uq*dtpdep+dv_dep
      dpoden = ovr23*uq*dtnden+dv_den
c
      dmadep = 2.0*dmpdep+2.0*dmndep-v_alfa*dpodep
      dmaden = 2.0*dmpden+2.0*dmnden-v_alfa*dpoden
c
      dnadep = alfdns*dmadep/t
      dnaden = alfdns*dmaden/t
c
      dsadep = dnadep*(2.5-mualfa/t)-alfdns*dmadep/t
      dsaden = dnaden*(2.5-mualfa/t)-alfdns*dmaden/t
c
      deadep = dnadep*(1.5*t-balpha)
      deaden = dnaden*(1.5*t-balpha)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c23456789012345678901234567890123456789012345678901234567890123456789012
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      s_dens = u_nuc*s_nuc+exclu*exalfa*s_out+exclu*s_alfa+s_sc+s_tr
c
      e_dens = u_nuc*e_nuc+exclu*exalfa*e_out+exclu*e_alfa+e_sc+e_tr
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
c                 ------------------------------------
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 !      temperature derivatives     !
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 ------------------------------------
c
      dna_dt = dnadt+dnadep*dep_dt+dnaden*den_dt
c
c
      dbsdt = (du_dt*s_nuc-
     1    du_dt*exalfa*s_out-exclu*v_alfa*dna_dt*s_out
     2    -du_dt*s_alfa+
     3    u_nuc*(dsidt+dsidx*dx_dt+dsidn*dni_dt)+
     4    exclu*exalfa*(dsodt+dsodep*dep_dt+dsoden*den_dt)+
     5    exclu*(dsadt+dsadep*dep_dt+dsaden*den_dt)+
     6    dsscdt+dsscdu*du_dt+dsscdx*dx_dt+dsscdn*dni_dt+
     7    dstrdt+dstrdu*du_dt+dstrdx*dx_dt+dstrdn*dni_dt)/brydns
c
c
c~~~~~~~~~~~~~~~~~~
c
      dbudt = t*dbsdt
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbfdt = dbudt-s_dens/brydns-t*dbsdt
c
c~~~~~~~~~~~~~~~~~~
c
      dbmudt = ye*(dmpodt+dmpdep*dep_dt+dmpden*den_dt)+
     1    (1.0-ye)*(dmnodt+dmndep*dep_dt+dmnden*den_dt)
c
c~~~~~~~~~~~~~~~~~~
c
      dbpdt = brydns*(dbmudt-dbfdt)
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c
c                 ------------------------------------
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 !       density derivatives        !
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 ------------------------------------
c
c
      dna_dn = dnadep*dep_dn+dnaden*den_dn
c
c
      dbsdn = (du_dn*s_nuc-
     1    du_dn*exalfa*s_out-exclu*v_alfa*dna_dn*s_out-du_dn*s_alfa+
     2    u_nuc*(dsidx*dx_dn+dsidn*dni_dn)+
     3    exclu*exalfa*(dsodep*dep_dn+dsoden*den_dn)+
     4    exclu*(dsadep*dep_dn+dsaden*den_dn)+
     5    dsscdu*du_dn+dsscdx*dx_dn+dsscdn*dni_dn+
     6    dstrdu*du_dn+dstrdx*dx_dn+dstrdn*dni_dn)/brydns-
     7    s_dens/brydns**2
c
c
c
c~~~~~~~~~~~~~~~~~~
c
      dbudn = (du_dn*e_nuc-
     1    du_dn*exalfa*e_out-exclu*v_alfa*dna_dn*e_out-du_dn*e_alfa+
     2    u_nuc*(deidx*dx_dn+deidn*dni_dn)+
     3    exclu*exalfa*(deodep*dep_dn+deoden*den_dn)+
     4    exclu*(deadep*dep_dn+deaden*den_dn)+
     5    descdu*du_dn+descdx*dx_dn+descdn*dni_dn+
     6    detrdu*du_dn+detrdx*dx_dn+detrdn*dni_dn)/brydns-
     7    e_dens/brydns**2
c
c
c
c
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbfdn = dbudn-t*dbsdn
c
c~~~~~~~~~~~~~~~~~~
c
      dbmudn = ye*(dmpdep*dep_dn+dmpden*den_dn)+
     1    (1.0-ye)*(dmndep*dep_dn+dmnden*den_dn)
c
c~~~~~~~~~~~~~~~~~~
c
      dbpdn = brydns*(dbmudn-dbfdn)+mubary-bftot
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c
c                 ------------------------------------
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 !         ye derivatives           !
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 ------------------------------------
c
c
c
c
      dna_dy = dnadep*dep_dy+dnaden*den_dy
c
c
      dbsdy = (du_dy*s_nuc-
     1    du_dy*exalfa*s_out-exclu*v_alfa*dna_dy*s_out-du_dy*s_alfa+
     2    u_nuc*(dsidx*dx_dy+dsidn*dni_dy)+
     3    exclu*exalfa*(dsodep*dep_dy+dsoden*den_dy)+
     4    exclu*(dsadep*dep_dy+dsaden*den_dy)+
     5    dsscdu*du_dy+dsscdx*dx_dy+dsscdn*dni_dy+
     6    dstrdu*du_dy+dstrdx*dx_dy+dstrdn*dni_dy)/brydns
c
c
c~~~~~~~~~~~~~~~~~~
c
      dbudy = (du_dy*e_nuc-
     1    du_dy*exalfa*e_out-exclu*v_alfa*dna_dy*e_out-du_dy*e_alfa+
     2    u_nuc*(deidx*dx_dy+deidn*dni_dy)+
     3    exclu*exalfa*(deodep*dep_dy+deoden*den_dy)+
     4    exclu*(deadep*dep_dy+deaden*den_dy)+
     5    descdu*du_dy+descdx*dx_dy+descdn*dni_dy+
     6    detrdu*du_dy+detrdx*dx_dy+detrdn*dni_dy)/brydns
c
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbfdy = dbudy-t*dbsdy
c
c~~~~~~~~~~~~~~~~~~
c
      dbmudy = ye*(dmpdep*dep_dy+dmpden*den_dy)+muprot+
     1    (1.0-ye)*(dmndep*dep_dy+dmnden*den_dy)-mun
c
c~~~~~~~~~~~~~~~~~~
c
      dbpdy = brydns*(dbmudy-dbfdy)
c
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c-----------------------------------------------------------------------
c                end of derivatives of thermodynamic variables
c-----------------------------------------------------------------------
c
c                  total derivatives
c                  (baryons+electrons+photons)
c
     
      dbudt = dbudt*ionadd
      deudt = deudt*eleadd
      dpudt = dpudt*radadd
      dudt  = dbudt + deudt + dpudt

      dbudn = dbudn*ionadd
      deudn = deudn*eleadd
      dpudn = dpudn*radadd
      dudn  = dbudn + deudn + dpudn

      dbudy = dbudy*ionadd
      deudy = deudy*eleadd
      dpudy = dpudy*radadd
      dudy  = dbudy + deudy + dpudy


      dbsdt = dbsdt*ionadd
      desdt = desdt*eleadd
      dpsdt = dpsdt*radadd
      dsdt  = dbsdt + desdt + dpsdt

      dbsdn = dbsdn*ionadd
      desdn = desdn*eleadd
      dpsdn = dpsdn*radadd
      dsdn  = dbsdn + desdn + dpsdn

      dbsdy = dbsdy*ionadd
      desdy = desdy*eleadd
      dpsdy = dpsdy*radadd
      dsdy  = dbsdy + desdy + dpsdy


      dbpdt = dbpdt*ionadd
      depdt = depdt*eleadd
      dppdt = dppdt*radadd
      dpdt  = dbpdt + depdt + dppdt

      dbpdn = dbpdn*ionadd
      depdn = depdn*eleadd
      dppdn = dppdn*radadd
      dpdn  = dbpdn + depdn + dppdn

      dbpdy = dbpdy*ionadd
      depdy = depdy*eleadd
      dppdy = dppdy*radadd
      dpdy  = dbpdy + depdy + dppdy

c
c
      dmudt = dbmudt + ye*demudt
      dmudn = dbmudn + ye*demudn
      dmudy = dbmudy + ye*demudy
c
c                calculate the adiabatic index
      gam_s = brydns*dpdn/ptot+t*(dpdt**2)/(brydns*ptot*dudt)
c
c
c                set the value of xprev to x for use the next 
c                time through
c
      xprev = x
c
c                save the value of the proton density to be used
c                by the "no nuclei" scheme on the next call
      p_prev = npout
c
c
c                return the three internal compositional variables
      inpvar(2) = nsubi
      inpvar(3) = eta_po
      inpvar(4) = eta_no
c
c
c
c                rejoice for this routine is finished!!!!!!!
 999  return
c
c
      end




c***********************************************************************
c
c    file:         alfeos.for
c
c***********************************************************************
c
c    module:       alfeos
c    type:         subroutine
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         8/30/90 modified from model 4-a
c
c                  please report any problems to me at:
c                  bitnet:  swesty@sunysbnp or
c                  internet: fswesty@astro.sunysb.edu
c                            fswesty@sbast3.sunysb.edu
c
c
c    call line:    call alfeos(inpvar,ye,brydns)
c
c    inputs:       inpvar = temp, internal eng, or entropy
c                  ye = electron fraction
c                  brydns = baryon number density
c
c    outputs:
c
c 
c    include files:  eos_m4c.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine alfeos(inpvar,ye,brydns,p_prev,ssflag)
c
      implicit none
c
c
      include 'eos_m4c.inc'
      include 'el_eos.inc'
c
c                       "zero" flag
      integer zflag
c
c                       "negative" flag
      integer nflag
c
c                       function type declarations
c
      double precision f_1_2, f_3_2, finv12, fhalfo, fhalf
c
c                       include the nucleon-nucleon interaction
c                       statement function definitions
      include 'force.inc'
      include 'const.dek'
c
c
c
c
c
c                         unset the zero flag
      zflag = 0
c
c                         ratio of baryon density to saturation density
      y = brydns/nsubs
c
c
c                         set t equal to the input variable (the entropy
c                         and internal energy options are not implemented
c                         in this version)
      t = inpvar(1)
c
c
c                         calc the quantum concentration of nucleons
      nq = 2.36d-4*t**1.5 
c
c                         calc the fermi integral coefficent
      uq = 20.721
c
      mq = (t/uq)**1.5
c
      kq = ((t/uq)**2.5)/(2.0*pi**2)
c
      lq = uq*(mq**ovr53)/(3.0*(pi**2))
c
      etamax = 0.95*finv12(2.0*(pi**2)*brydns/mq)
c
c
c
c                              set the proton density to its old value
      npout = p_prev
c
      if(brydns.gt.(0.98*2.0/(ye*v_alfa))) then
        npout = ye*brydns
        nnout = (1.0-ye)*brydns
        nout = brydns
c
c
c20        vnout = eiflag*(2.0*aa*nout+4.0*bb*npout+cc*(1.0+dd)*nout**dd)
        vnout = eiflag*pvn(npout,nnout)
c
c20        vpout = eiflag*(2.0*aa*nout+4.0*bb*nnout+
c20     1    cc*(1.0+dd)*nout**dd+deltam)
        vpout = eiflag*pvp(npout,nnout)
c
c
        zno = 2.0*(pi**2)*nnout/mq
c
        zpo = 2.0*(pi**2)*npout/mq
c
        eta_no = finv12(zno)
c
        eta_po = finv12(zpo)
c
        f32_no = f_3_2(eta_no)
        f32_po = f_3_2(eta_po)
c
        tau_no = kq*f32_no
        tau_po = kq*f32_po
c
c
c20        bprout = lq*(f32_po+f32_no)+eiflag*(
c20     1    aa*(nout**2)+4.0*bb*npout*nnout+dd*cc*(nout**(1.0+dd)) )
        bprout = lq*(f32_po+f32_no)+eiflag*pv_pr(npout,nnout)
c
c
        mun_o = t*eta_no+vnout
        mun = mun_o
c
        mup_o = t*eta_po+vpout
        muprot = mup_o
c
c                              calculate diff. of chem. potentials
        muhat = mun-muprot

c JCM: dimensionless chemical potentials for protons and neutrons
      etapls=muprot/t
      etanls=mun/t

c
c
c                              calculate the alpha particle
c                              chemical potential
        mualfa = 2.0*mun+2.0*muprot+balpha-bprout*v_alfa
c
        alfdns = 0.0
c
        exalfa = 1.0-alfdns*v_alfa
c
      else
c
c                              calculate the neutron density
        nnout = 2.0*brydns*(1.0-2.0*ye)/(2.0-brydns*ye*v_alfa)+
     1            npout*(2.0-(1.0-ye)*brydns*v_alfa)/
     2            (2.0-brydns*ye*v_alfa)
c
c                              calculate density of outside nucleons
        nout = npout+nnout
c
c20        vnout = eiflag*(2.0*aa*nout+4.0*bb*npout+cc*(1.0+dd)*nout**dd)
        vnout = eiflag*pvn(npout,nnout)
c
c20        vpout = eiflag*(2.0*aa*nout+4.0*bb*nnout+
c20     1    cc*(1.0+dd)*nout**dd+deltam)
        vpout = eiflag*pvp(npout,nnout)
c
c
        zno = 2.0*(pi**2)*nnout/mq
c
        zpo = 2.0*(pi**2)*npout/mq
c
        eta_no = finv12(zno)
c
        eta_po = finv12(zpo)
c
        f32_no = f_3_2(eta_no)
        f32_po = f_3_2(eta_po)
c
        tau_no = kq*f32_no
        tau_po = kq*f32_po
c
c
c20        bprout = lq*(f32_po+f32_no)+eiflag*(
c20     1    aa*(nout**2)+4.0*bb*npout*nnout+dd*cc*(nout**(1.0+dd)) )
        bprout = lq*(f32_po+f32_no)+eiflag*pv_pr(npout,nnout)
c
c
        mun_o = t*eta_no+vnout
        mun = mun_o
c
        mup_o = t*eta_po+vpout
        muprot = mup_o
c
c                              calculate diff. of chem. potentials
        muhat = mun-muprot

c JCM: dimensionless chemical potentials for protons and neutrons
      etapls=muprot/t
      etanls=mun/t

c
c
c                              calculate the alpha particle
c                              chemical potential
        mualfa = 2.0*mun+2.0*muprot+balpha-bprout*v_alfa
c
c                              calculate density of alpha particles
c
        if(abs(mualfa/t).lt.30.0) then
          alfdns = 8.0*nq*dexp(mualfa/t)
        elseif((mualfa/t).lt.-30.0) then
          alfdns = 0.0
        else
          alfdns = 8.0*nq*dexp(3.0d1)
        endif
c
c
        exalfa = 1.0-alfdns*v_alfa
c
c                              calculate "non-zeroness" of baryon
c                              conservation equation and save the
c                              value to be used in the finite
c                              difference approximation of dgdprt
        gold = brydns-exalfa*(nnout+npout)-4.0*alfdns
        prtold = npout
c
c                              take a small step to get derivative
        npout = npout+0.001*brydns
c
        do 11 i=1,30,1

c JCM:
c           write(*,*) 'looping2',i
c
c                              unset the negative flag
          nflag = 0
c
 14       continue
c                              calculate the neutron density
          nnout = 2.0*brydns*(1.0-2.0*ye)/(2.0-brydns*ye*v_alfa)+
     1            npout*(2.0-(1.0-ye)*brydns*v_alfa)/
     2            (2.0-brydns*ye*v_alfa)
c     

c JCM:
c           write(*,*) 'looping2b',i

          if((nnout.lt.0.0).and.(i.eq.1)) then
            npout = prtold-0.5*dprt
          elseif((nnout.lt.0.0).and.(i.eq.1).and.(nflag.ne.1)) then
            npout = 0.99*p_prev
            nflag = 1
          elseif((nnout.lt.0.0).and.(i.eq.1).and.(nflag.eq.1)) then
            ssflag = 0
            goto 999
          endif
c                              calculate density of outside nucleons
          nout = npout+nnout
c
c20          vnout = eiflag*
c20     1      (2.0*aa*nout+4.0*bb*npout+cc*(1.0+dd)*nout**dd)
          vnout = eiflag*pvn(npout,nnout)
c
c20          vpout = eiflag*(2.0*aa*nout+4.0*bb*nnout+
c20     1      cc*(1.0+dd)*nout**dd+deltam)
          vpout = eiflag*pvp(npout,nnout)
c
c
          zno = 2.0*(pi**2)*nnout/mq
c
          zpo = 2.0*(pi**2)*npout/mq
c
          eta_no = finv12(zno)
c
          eta_po = finv12(zpo)
c
          f32_no = f_3_2(eta_no)
c
          f32_po = f_3_2(eta_po)
c
          tau_no = kq*f32_no
          tau_po = kq*f32_po
c
c20          bprout = lq*(f32_po+f32_no)+eiflag*
c20     1      (aa*(nout**2)+4.0*bb*npout*nnout+dd*cc*(nout**(1.0+dd)) )
          bprout = lq*(f32_po+f32_no)+eiflag*pv_pr(npout,nnout)
c
          mun_o = t*eta_no+vnout
          mun = mun_o
c
          mup_o = t*eta_po+vpout
          muprot = mup_o
c
c                              calc difference of potentials
          muhat = mun-muprot

c JCM: dimensionless chemical potentials for protons and neutrons
      etapls=muprot/t
      etanls=mun/t

c
c                              calc alpha particle chemical potentials
          mualfa = 2.0*mun+2.0*muprot+balpha-bprout*v_alfa
c
c                              calc alpha particle density
c
          if(abs(mualfa/t).lt.30.0) then
            alfdns = 8.0*nq*dexp(mualfa/t)
          elseif((mualfa/t).lt.-30.0) then
            alfdns = 0.0
          else
            alfdns = 8.0*nq*dexp(3.0d1)
          endif
c
c
          exalfa = 1.0-alfdns*v_alfa
c
c                              calc "non-zeroness" of baryon cons. eqn.
          g = brydns-exalfa*(nnout+npout)-4.0*alfdns
c
c                              calculate derivative of baryon conservation
c                              equation w.r.t. proton density by finite
c                              diference approximation
          dgdprt = (g-gold)/(npout-prtold)
c
c                              if rate of change is near zero
c                              and zero flag is not set
          if((abs(dgdprt).lt.1.0d-25).and.(zflag.eq.0)) then
c                              tweak the step size
            npout = prtold-0.5*dprt
c                              and set the zero flag
            zflag = 1
c                              and go back and try again
            goto 14
c                              if failure occurs again
          elseif((abs(dgdprt).lt.1.0d-25).and.(zflag.eq.1)) then
c                              declare an eos failure
            ssflag = 0
c                              and return
            goto 999
          endif
c
c                              calculate new newton-raphson step
          dprt = g/dgdprt
c
c                              save old value of proton density & g
          prtold = npout
          gold = g
c
c
 13       continue
c
c     JCM: GODMARK
c           write(*,*) 'looping2c',i,dprt
c                              potential "new" value of proton density
          prtnew = npout-dprt
c
c JCM: lieb insert next 4 lines
c          if (DPRT.lt.1.d-199) then
c             SSFLAG = 0.
c             goto 999
c          endif
c                              if new proton density is less than the
c                              baryon density and greater than zero 
c                              then update the proton density
          if(prtnew*(brydns-prtnew).gt.0.0) then
            npout = npout-dprt
c                              else cut the step size in half and try again
          else
            dprt = dprt*0.5
c     JCM: infinite loop occurs here
            goto 13
          endif
c
c                              if step size is small enough break out of
c                              the do 11 loop, otherwise continue
          if(abs(dprt/npout).lt.10e-11) goto 12

c JCM:
c           write(*,*) 'looping2d',i

 11     continue
c
c      write(*,*) 'a failed to converge; switching to f' ! take out later
        ssflag = 0
        goto 999
c
c JCM:
c           write(*,*) 'looping2e',i
c
 12     continue
c
      endif
c                              set the success flag
      ssflag = 1
c
c
c                              calc outside nucleon density
      nout = nnout+npout
c
c                              calc outside nucleon fraction
      xout = npout/nout
c
c                              calculate particle fractions
      xalfa = 4.0*alfdns/brydns
      xprot = exalfa*npout/brydns
      xnut = exalfa*nnout/brydns
      xh = 0.0
c
c                              baryons
c
      f32_no = f_3_2(eta_no)
c
      f32_po = f_3_2(eta_po)
c
      tau_po = kq*f32_po
c
      tau_no = kq*f32_no
c
c
c
c
c
c
c                    calculate internal energy of outside nucleons
c20      buout = (xnut+xprot)*( uq*(tau_po+tau_no)+
c20     1    eiflag*((nout**2)*aa+
c20     2   4.0*bb*npout*nnout+cc*nout**(1.0+dd)+npout*deltam) )/nout
      buout = (xnut+xprot)*( uq*(tau_po+tau_no)+
     1    eiflag*pv_e(npout,nnout) )/nout
c
c
c                                calc alfa particle internal energy
      bualfa = 0.25*xalfa*(1.5*t-balpha)
c
c                                set nuclei internal energy to zero
      bunuc = 0.0

c                                calculate total baryon internal energy
c                                (per baryon)
      bu = buout+bualfa+bunuc
c
c
c                                calc entropy of outside nucleons
      bsout = (xnut+xprot)*( (5.0*uq/(3.0*t))*(tau_no+tau_po)-
     1   nnout*eta_no-npout*eta_po )/nout
c
c                                calc alpha particle entropy
      bsalfa = 0.25*xalfa*(2.5-mualfa/t)
c
c                                set nuclei entropy to zero
      bsnuc = 0.0
c
c                                calc total baryon entropy (per baryon)
      bs = bsout+bsalfa+bsnuc
c
c
c
c                                calc outside free energy
      bfout = buout-t*bsout
c                                calc alpha particle free energy
      bfalfa = bualfa-t*bsalfa
c                                set nuclei free energy to zero
      bfnuc = bunuc-t*bsnuc
c                                calc total baryon free energy (per nucleon)
      bftot = bfout+bfalfa+bfnuc
c
c
c
c
c
c                                calc outside pressure
c20      bprout = lq*(f32_po+f32_no)+eiflag*(
c20     1    aa*(nout**2)+4.0*bb*npout*nnout+dd*cc*(nout**(1.0+dd)))
      bprout = lq*(f32_po+f32_no)+eiflag*pv_pr(npout,nnout)
c
c                                calc alfa particle pressure
      bpralf = alfdns*t
c
c                                set nuclei pressure to zero
      bprnuc = 0.0
c
c                                calc total baryon pressure
      bpress = bprout+bpralf+bprnuc
c
c
c
c
c
c
c
c
c                           leptons & photons

c JCM:
c      write(*,*) 'before anyelectron2'
      call any_electron(t,ye,brydns) 

c      call el_eos(t,ye,brydns)


c
c
c
c                           total pressure and eng/ent per baryon
c
      fbary = bftot+fsube
      pbary = bpress+epress
      mubary = ye*muprot+(1.0-ye)*mun
      mu_mat = ye*(muprot+musube)+(1.0-ye)*mun
c

      bftot  = bftot*ionadd
      bu     = bu*ionadd
      bs     = bs*ionadd
      bpress = bpress*ionadd

      fsube  = fsube*eleadd
      eu     = eu*eleadd
      es     = es*eleadd
      epress = epress*eleadd

      pf     = pf*radadd
      pu     = pu*radadd
      ps     = ps*radadd
      ppress = ppress*radadd

c..totals
      ftot = bftot  + fsube  + pf
      utot = bu     + eu     + pu
      stot = bs     + es     + ps
      ptot = bpress + epress + ppress


c
c
c
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c                derivatives of thermodynamic variables
c-----------------------------------------------------------------------
c
c
c
cc      gpo = 2.0*fhalfo(eta_po)
cc      gno = 2.0*fhalfo(eta_no)
c
c
      gpo = 2.0*fhalf(eta_po)
      gno = 2.0*fhalf(eta_no)
c
c
c                 ------------------------------------
c                 !      derivatives of exterior     !
c                 !      quantities                  !
c                 !      (w.r.t. temp. and eta's)    !
c                 !                                  !
c                 ------------------------------------
c
c
c                  derivatives of exterior potentials
c                  w.r.t. particle densities
c20      dvpodp = eiflag*(2.0*aa+dd*(1.0+dd)*cc*(nout**(dd-1.0)) )
      dvpodp = eiflag*dpvpdp(npout,nnout)
c20      dvpodn = eiflag*(2.0*aa+4.0*bb+dd*(1.0+dd)*cc*(nout**(dd-1.0)) )
      dvpodn = eiflag*dpvpdn(npout,nnout)
c20      dvnodp = dvpodn
      dvnodp = eiflag*dpvndp(npout,nnout)
c20      dvnodn = dvpodp
      dvnodn = eiflag*dpvndn(npout,nnout)
c
c
c                  derviatives of exterior chem. pot. w.r.t. eta's
c                  (at fixed t)
      dmpdep = t+dvpodp*npout/gpo
      dmpden = dvpodn*nnout/gno
      dmndep = dvnodp*npout/gpo
      dmnden = t+dvnodn*nnout/gno
c
c                  derivatives of pressure potential w.r.t.
c                  particle densities
c20      dv_dpo = eiflag*
c20     1    (2.0*aa*nout+4.0*bb*nnout+cc*dd*(1.0+dd)*(nout**dd) )
      dv_dpo = eiflag*dpvrdp(npout,nnout)
c20      dv_dno = eiflag*
c20     1    (2.0*aa*nout+4.0*bb*npout+cc*dd*(1.0+dd)*(nout**dd) )
      dv_dno = eiflag*dpvrdn(npout,nnout)
c
c                  derivatives of pressure potential w.r.t. eta's
c                  (at fixed t)
      dv_dep = dv_dpo*npout/gpo
      dv_den = dv_dno*nnout/gno
c
c                  derivatives of outside pressure w.r.t. eta's
c                  (at fixed t)
      dpodep = npout*t+dv_dep
      dpoden = nnout*t+dv_den
c
c                  derivatives of alpha density w.r.t. eta's
c                  (at fixed t)
      dnadep = alfdns*(2.0*dmpdep+2.0*dmndep-v_alfa*dpodep)/t
      dnaden = alfdns*(2.0*dmpden+2.0*dmnden-v_alfa*dpoden)/t
c
c
c                  derivatives of particle densities w.r.t. t
c                  (at fixed eta's)
      dnpodt = 1.5*npout/t
      dnnodt = 1.5*nnout/t
c
c                  derivatives of exterior chem. pot. w.r.t. t
c                  (at fixed eta's)
      dmpodt = eta_po+dvpodp*dnpodt+dvpodn*dnnodt
      dmnodt = eta_no+dvnodp*dnpodt+dvnodn*dnnodt
c
c                  derivative of pressure potential w.r.t. t
c                  (at fixed eta's)
      dv_dt = dv_dpo*dnpodt+dv_dno*dnnodt
c
c                  derivative of outside pressure w.r.t. t
c                  (at fixed eta's)
      dpodt = ovr23*uq*2.5*(tau_po+tau_no)/t+dv_dt
c
c                  derivative of alpha chem. pot. w.r.t. t
c                  (at fixed eta's)
      dmuadt = 2.0*dmpodt+2.0*dmnodt-v_alfa*dpodt
c
c                  derivative of alpha particle density w.r.t. t
c                  (at fixed eta's)
      dnadt = 1.5*alfdns/t-alfdns*mualfa/(t**2)+alfdns*dmuadt/t
c
c                  derivative of alpha particle pressure w.r.t. t
c                  (at fixed eta's)
      dpadt = alfdns+t*dnadt
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c
c                 ------------------------------------
c                 !      derivatives of constraint   !
c                 !      and equilibrium equations   !
c                 !      with respect to the five    !
c                 !      compositional variables     !
c                 !      (u,x,n_i,eta_po,eta_no)     !
c                 !      and the three independent   !
c                 !      variables                   !
c                 !      (baryon density, t, and ye) !
c                 !                                  !
c                 ------------------------------------
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                equation 1 (baryon conservation)
c
c
c
      dg1do1 = -exalfa*npout/gpo+(v_alfa*nout-4.0)*dnadep
c
      dg1do2 = -exalfa*nnout/gno+(v_alfa*nout-4.0)*dnaden
c
c
      dg1dl1 = 1.0
c
      dg1dl2 = -exalfa*(dnnodt+dnpodt)+(v_alfa*nout-4.0)*dnadt
c
      dg1dl3 = 0.0
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                equation 2 (charge conservation)
c
c
      dg2do1 = -exalfa*npout/gpo+(v_alfa*npout-2.0)*dnadep
c
      dg2do2 = (v_alfa*npout-2.0)*dnaden
c
c
      dg2dl1 = ye
c
      dg2dl2 = -exalfa*dnpodt+(v_alfa*npout-2.0)*dnadt
c            
      dg2dl3 = brydns
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      det_gt = dg1do1*dg2do2-dg1do2*dg2do1
c
c
      dep_dn = (dg1do2*dg2dl1-dg2do2*dg1dl1)/det_gt
      den_dn = (dg2do1*dg1dl1-dg1do1*dg2dl1)/det_gt
c
c
      dep_dt = (dg1do2*dg2dl2-dg2do2*dg1dl2)/det_gt
      den_dt = (dg2do1*dg1dl2-dg1do1*dg2dl2)/det_gt
c
c
      dep_dy = (dg1do2*dg2dl3-dg2do2*dg1dl3)/det_gt
      den_dy = (dg2do1*dg1dl3-dg1do1*dg2dl3)/det_gt
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                 -------------------------------------
c                 !      derivatives of exterior bulk !
c                 !      nucleon internal energy &    !
c                 !      entropy densities and the    !
c                 !      chem. pot.  w.r.t. to eta_p, !
c                 !      ate_n & t. these are used in !
c                 !      calculating the derivatives  !
c                 !      w.r.t. the independant vars  !
c                 !      (baryon density, t, and ye)  !
c                 !                                   !
c                 -------------------------------------
c
c
      s_out =(ovr53*uq/t)*(tau_no+tau_po)-nnout*eta_no-npout*eta_po
c
c20      e_out = uq*(tau_po+tau_no)+eiflag*
c20     1((nout**2)*aa+4.0*bb*npout*nnout+cc*nout**(1.0+dd)+npout*deltam)
      e_out = uq*(tau_po+tau_no)+eiflag*pv_e(npout,nnout)
c
c                   derivative of exterior entropy density w.r.t. t
      dsodt =  ovr53*uq*(1.5*(tau_po+tau_no)/(t**2))-
     1     1.5*(npout*eta_po+nnout*eta_no)/t
c
      deodt = t*dsodt
c
c                    derivatives of exterior particle densities w.r.t.
c                    temperature (eta's fixed)
      dnpodt = 1.5*npout/t
      dnnodt = 1.5*nnout/t
c
      dmpodt = eta_po+dvpodp*dnpodt+dvpodn*dnnodt
      dmnodt = eta_no+dvnodp*dnpodt+dvnodn*dnnodt
c
c
      dnpdep = npout/gpo
      dnnden = nnout/gno
c
      dtpdep = 1.5*t*npout/uq
      dtnden = 1.5*t*nnout/uq
c
      dsodep = (ovr53*uq/t)*dtpdep-npout-eta_po*dnpdep
      dsoden = (ovr53*uq/t)*dtnden-nnout-eta_no*dnnden
c
c
c                    exterior particle potentials
c20      vnout = eiflag*(2.0*aa*nout+4.0*bb*npout+cc*(1.0+dd)*nout**dd )
      vnout = eiflag*pvn(npout,nnout)
c20      vpout = eiflag*
c20     1    (2.0*aa*nout+4.0*bb*nnout+cc*(1.0+dd)*nout**dd+deltam)
      vpout = eiflag*pvp(npout,nnout)
c
c
      deodep = uq*dtpdep+vpout*dnpdep
      deoden = uq*dtnden+vnout*dnnden
c
      dmpdep = t+dvpodp*npout/gpo
      dmpden = dvpodn*nnout/gno
      dmndep = dvnodp*npout/gpo
      dmnden = t+dvnodn*nnout/gno
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                 -------------------------------------
c                 !      derivatives of alpha         !
c                 !      particle internal energy &   !
c                 !      entropy densities and the    !
c                 !      chem. pot.  w.r.t. to eta_p, !
c                 !      ate_n & t. these are used in !
c                 !      calculating the derivatives  !
c                 !      w.r.t. the independant vars  !
c                 !      (baryon density, t, and ye)  !
c                 !                                   !
c                 -------------------------------------
c
c
c
c
      s_alfa = alfdns*(2.5-mualfa/t)
c
c
      e_alfa = alfdns*(1.5*t-balpha)
c
c                  derivative of pressure potential w.r.t. t
      dv_dt = dv_dpo*dnpodt+dv_dno*dnnodt
c
c                  derivative of outside pressure w.r.t. t
      dpodt = ovr23*uq*2.5*(tau_po+tau_no)/t+dv_dt
c
c
      dmuadt = 2.0*dmpodt+2.0*dmnodt-v_alfa*dpodt
c
c                  derivative of alpha particle density w.r.t. t
      dnadt = 1.5*alfdns/t-alfdns*mualfa/(t**2)+alfdns*dmuadt/t
c
c
      dsadt = dnadt*(2.5-mualfa/t)-alfdns*dmuadt/t+alfdns*mualfa/t**2
c
      deadt = dnadt*(1.5*t-balpha)+1.5*alfdns
c
c
      dv_dep = dv_dpo*npout/gpo
      dv_den = dv_dno*nnout/gno
c
      dpodep = ovr23*uq*dtpdep+dv_dep
      dpoden = ovr23*uq*dtnden+dv_den
c
      dmadep = 2.0*dmpdep+2.0*dmndep-v_alfa*dpodep
      dmaden = 2.0*dmpden+2.0*dmnden-v_alfa*dpoden
c
      dnadep = alfdns*dmadep/t
      dnaden = alfdns*dmaden/t
c
      dsadep = dnadep*(2.5-mualfa/t)-alfdns*dmadep/t
      dsaden = dnaden*(2.5-mualfa/t)-alfdns*dmaden/t
c
      deadep = dnadep*(1.5*t-balpha)
      deaden = dnaden*(1.5*t-balpha)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c23456789012345678901234567890123456789012345678901234567890123456789012
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      s_dens = exalfa*s_out+s_alfa
c
      e_dens = exalfa*e_out+e_alfa
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
c                 ------------------------------------
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 !      temperature derivatives     !
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 ------------------------------------
c
      dna_dt = dnadt+dnadep*dep_dt+dnaden*den_dt
c
c
      dbsdt = (-v_alfa*dna_dt*s_out+
     1    exalfa*(dsodt+dsodep*dep_dt+dsoden*den_dt)+
     2    (dsadt+dsadep*dep_dt+dsaden*den_dt) )/brydns
c
c~~~~~~~~~~~~~~~~~~
c
      dbudt = t*dbsdt
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbfdt = dbudt-s_dens/brydns-t*dbsdt
c
c~~~~~~~~~~~~~~~~~~
c
      dbmudt = ye*(dmpodt+dmpdep*dep_dt+dmpden*den_dt)+
     1    (1.0-ye)*(dmnodt+dmndep*dep_dt+dmnden*den_dt)
c
c~~~~~~~~~~~~~~~~~~
c
      dbpdt = brydns*(dbmudt-dbfdt)
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c
c                 ------------------------------------
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 !       density derivatives        !
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 ------------------------------------
c
c
      dna_dn = dnadep*dep_dn+dnaden*den_dn
c
c
      dbsdn = (-v_alfa*dna_dn*s_out+
     1    exalfa*(dsodep*dep_dn+dsoden*den_dn)+
     2   (dsadep*dep_dn+dsaden*den_dn) )/brydns-s_dens/brydns**2
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbudn = (-v_alfa*dna_dn*e_out+
     1    exalfa*(deodep*dep_dn+deoden*den_dn)+
     2   (deadep*dep_dn+deaden*den_dn) )/brydns-e_dens/brydns**2
c
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbfdn = dbudn-t*dbsdn
c
c~~~~~~~~~~~~~~~~~~
c
      dbmudn = ye*(dmpdep*dep_dn+dmpden*den_dn)+
     1    (1.0-ye)*(dmndep*dep_dn+dmnden*den_dn)
c
c~~~~~~~~~~~~~~~~~~
c
      dbpdn = brydns*(dbmudn-dbfdn)+mubary-bftot
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c
c                 ------------------------------------
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 !         ye derivatives           !
c                 !                                  !
c                 !                                  !
c                 !                                  !
c                 ------------------------------------
c
c
c
c
      dna_dy = dnadep*dep_dy+dnaden*den_dy
c
c
      dbsdy = (-v_alfa*dna_dy*s_out+
     1    exalfa*(dsodep*dep_dy+dsoden*den_dy)+
     2   (dsadep*dep_dy+dsaden*den_dy) )/brydns
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbudy = (-v_alfa*dna_dy*e_out+
     1    exalfa*(deodep*dep_dy+deoden*den_dy)+
     2   (deadep*dep_dy+deaden*den_dy) )/brydns
c
c
c~~~~~~~~~~~~~~~~~~
c
c
      dbfdy = dbudy-t*dbsdy
c
c~~~~~~~~~~~~~~~~~~
c
      dbmudy = ye*(dmpdep*dep_dy+dmpden*den_dy)+muprot+
     1    (1.0-ye)*(dmndep*dep_dy+dmnden*den_dy)-mun
c
c~~~~~~~~~~~~~~~~~~
c
      dbpdy = brydns*(dbmudy-dbfdy)
c
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c
c-----------------------------------------------------------------------
c                end of derivatives of thermodynamic variables
c-----------------------------------------------------------------------
c
c                  total derivatives
c                  (baryons+electrons+photons)
c
czzzz


      dbudt = dbudt*ionadd
      deudt = deudt*eleadd
      dpudt = dpudt*radadd
      dudt  = dbudt + deudt + dpudt

      dbudn = dbudn*ionadd
      deudn = deudn*eleadd
      dpudn = dpudn*radadd
      dudn  = dbudn + deudn + dpudn

      dbudy = dbudy*ionadd
      deudy = deudy*eleadd
      dpudy = dpudy*radadd
      dudy  = dbudy + deudy + dpudy


      dbsdt = dbsdt*ionadd
      desdt = desdt*eleadd
      dpsdt = dpsdt*radadd
      dsdt  = dbsdt + desdt + dpsdt

      dbsdn = dbsdn*ionadd
      desdn = desdn*eleadd
      dpsdn = dpsdn*radadd
      dsdn  = dbsdn + desdn + dpsdn

      dbsdy = dbsdy*ionadd
      desdy = desdy*eleadd
      dpsdy = dpsdy*radadd
      dsdy  = dbsdy + desdy + dpsdy


      dbpdt = dbpdt*ionadd
      depdt = depdt*eleadd
      dppdt = dppdt*radadd
      dpdt  = dbpdt + depdt + dppdt

      dbpdn = dbpdn*ionadd
      depdn = depdn*eleadd
      dppdn = dppdn*radadd
      dpdn  = dbpdn + depdn + dppdn

      dbpdy = dbpdy*ionadd
      depdy = depdy*eleadd
      dppdy = dppdy*radadd
      dpdy  = dbpdy + depdy + dppdy

c
      dmudt = dbmudt+ye*demudt
      dmudn = dbmudn+ye*demudn
      dmudy = dbmudy+ye*demudy
c
c
c                  adiabatic index
      gam_s = brydns*dpdn/ptot+t*(dpdt**2)/(brydns*ptot*dudt)
c
c
      inpvar(2) = nsubs
      inpvar(3) = eta_po
      inpvar(4) = eta_no
c
c
c                           approximate the nuclear density
      nsubi = nsubs
c
c                           use 0.45 as the nuclear proton fraction
      x = 0.45
      a = 4.0d0
c
c                           save the proton number density for use
c                           as the initial guess on next call
      p_prev = npout
c
c
 999  return
c
      end
calfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalf
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         maxwel.for
c
c***********************************************************************
c
c    module:       maxwel
c    type:         subroutine
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         7/13/90
c
c                  please report any problems to me at:
c                  bitnet:  swesty@sunysbnp or
c                  internet: fswesty@astro.sunysb.edu
c
c
c    call line:    call maxwel(inpvar,ye,brydns)
c
c    inputs:       inpvar = temp, internal eng, or entropy
c                  ye = electron fraction
c                  brydns = baryon number density
c
c    outputs:
c
c
c
c 
c    include files:  eos_m4c.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine maxwel(inpvar,ye,brydns,xprev,p_prev,ssflag)
c
      implicit none
c
      double precision outvar(4)
      double precision dpn_dt, dpn_dn, dpn_dy
      double precision dsn_dt, dsn_dn, dsn_dy
      double precision dsb_dt, dsb_dn, dsb_dy
      double precision dmu_dt, dmu_dn, dmu_dy
      double precision dphadt, dphady, deldns
      double precision n_xh, n_xa, n_xn, n_xp, b_xa, b_xn, b_xp
c
c
      include 'eos_m4c.inc'
      include 'el_eos.inc'
      include 'maxwel.inc'
c
c
c                   set the temperature
      t = inpvar(1)
c
c
c                   calculate and save chem. pot. and thermodynamic
c                   quantaties from low end of two phase region
      call nuceos(inpvar,ye,lowdns,xprev,p_prev,ssflag)
c
c
c
c                    if the nuclear eos failed and the reset flag is set
c                    then reset the initial guesses and try again
      if((ssflag.ne.1).and.(rsflag.eq.1)) then
        call reset(inpvar,ye,lowdns,outvar)
        outvar(1) = inpvar(1)
        call nuceos(outvar,ye,lowdns,xprev,p_prev,ssflag)
c
c
c                    make a last ditch effort at convergence
        if(ssflag.ne.1) then
          outvar(2) = 0.155
          outvar(3) = -15.0
          outvar(4) = -20.0
          call nuceos(outvar,ye,lowdns,xprev,p_prev,ssflag)
        else
          inpvar(2) = outvar(2)
          inpvar(3) = outvar(3)
          inpvar(4) = outvar(4)
        endif
c
      endif
c
c
c
c
      prlow = ptot-ppress
      s_low = stot-ps
      f_low = ftot-pf
      mutlow = (1.0-ye)*mun+ye*(muprot+musube)
      muelow = musube
      muhlow = muhat
c
      dpn_dt = dpdt
      dpn_dn = dpdn
      dpn_dy = dpdy
c
      dmu_dt = dmudt
      dmu_dn = dmudn
      dmu_dy = dmudy
c
      dsn_dt = dsdt-dpsdt
      dsn_dn = dsdn
      dsn_dy = dsdy
c
      n_xh = xh
      n_xa = xalfa
      n_xp = xprot
      n_xn = xnut
c
c
      if(ssflag.ne.1) then
        write(*,*) 'maxwel:  nuclear eos failed at try:'
        write(*,*) t,lowdns,ye
        write(*,*) inpvar
        goto 999
      endif
c                   calculate and save chem. pot. and thermodynamic
c                   quantaties from high end of two phase region
      call alfeos(inpvar,ye,hidns,p_prev,ssflag)
c
      prhi = ptot-ppress
      s_hi = stot-ps
      f_hi = ftot-pf
      muthi = (1.0-ye)*mun+ye*(muprot+musube)
      muehi = musube
      muhhi = muhat
c
c
      dsb_dt = dsdt-dpsdt
      dsb_dn = dsdn
      dsb_dy = dsdy
c
c
      b_xa = xalfa
      b_xp = xprot
      b_xn = xnut
c
c
      if(ssflag.ne.1) then
        write(*,*) 'maxwel:  alfa eos failed at try:'
        write(*,*) t,hidns,ye
        write(*,*) inpvar
        goto 999
      endif
c
c                   calculate "average" chem. pot. and pressure
c                   in order to avoid numerical problems
      mutild = (mutlow+muthi)/2.0
      prtild = (prlow+prhi)/2.0
c
c                   calculate phase fraction
      phasef = (brydns-lowdns)/(hidns-lowdns)
c
c
c                   electron number density
      nsube = brydns*ye
c
c                   call electron eos to determine the
c                   electron chemical potential

c JCM:
c      write(*,*) 'before anyelectron3'
      call any_electron(t,ye,brydns) 

c      call el_eos(t,ye,brydns)



c
c
      muhat = musube+(1.0-phasef)*(muhlow-muelow)+phasef*(muhhi-muehi)
c
      mun = mutild+ye*(muhat-musube)
c
      muprot = mun-muhat

c JCM: dimensionless chemical potentials for protons and neutrons
      etapls=muprot/t
      etanls=mun/t

c
c                   calculate thermodynamic quantities
c
      stot = ((1.0-phasef)*s_low*lowdns+phasef*s_hi*hidns)/brydns+ps
c
      ftot = (lowdns*f_low+mutild*(brydns-lowdns))/brydns+pf
c
      utot = ftot+t*stot+pu
c
      ptot = prtild+ppress
c
c
      xh = (1.0-phasef)*n_xh
      xalfa = (1.0-phasef)*n_xa
      xnut = (1.0-phasef)*n_xn
      xprot = (1.0-phasef)*n_xp
      xalfa2 = phasef*b_xa
      xnut2 = phasef*b_xn
      xprot2 = phasef*b_xp
c
c
c
c
      deldns = hidns-lowdns
c
c
      dphadt = ((brydns-lowdns)/deldns**2-1.0/deldns)*dnl_dt-
     1    ((brydns-lowdns)/deldns**2)*dnh_dt
c
      dpdt = dpn_dt+dpn_dn*dnl_dt
      dmudt = dmu_dt+dmu_dn*dnl_dt
      dsdt = (1.0-phasef)*lowdns*(dsn_dt+dsn_dn*dnl_dt)/brydns+
     2 (1.0-phasef)*s_low*dnl_dt/brydns-lowdns*s_low*dphadt/brydns+
     3    (dphadt*s_hi*hidns+phasef*dnh_dt*s_hi+
     4    phasef*hidns*(dsb_dt+dsb_dn*dnh_dt))/brydns+dpsdt
      dudt = dmudt-dpdt/brydns+stot+t*dsdt
c 
c
      dpdn = 0.0
      dmudn = 0.0
      dsdn = -dpdt/brydns**2
      dudn = (lowdns*(mutild-ftot)/brydns**2)+t*dsdn
c
c
      dphady = ((brydns-lowdns)/deldns**2-1.0/deldns)*dnl_dy-
     1    ((brydns-lowdns)/deldns**2)*dnh_dy
c
      dpdy = dpn_dy+dpn_dn*dnl_dy
      dmudy = dmu_dy+dmu_dn*dnl_dy
      dsdy = (1.0-phasef)*lowdns*(dsn_dy+dsn_dn*dnl_dy)/brydns+
     2 (1.0-phasef)*s_low*dnl_dy/brydns-lowdns*s_low*dphady/brydns+
     3    (dphady*s_hi*hidns+phasef*dnh_dy*s_hi+
     4    phasef*hidns*(dsb_dy+dsb_dn*dnh_dy))/brydns
      dudy = dmudy-dpdy/brydns+t*dsdy
c
c
c
c
c             adiabatic index
c             (note that the first term vanishes in this expression)
      gam_s = t*(dpdt**2)/(brydns*ptot*dudt)
c
c
 999  return
c
c
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         eoslog.for
c
c***********************************************************************
c
c    module:       eoslog
c    type:         subroutine
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         12/15/90
c
c                  please report any problems to me at:
c                  bitnet:  swesty@sunysbnp or
c                  internet: fswesty@astro.sunysb.edu or
c                            fswesty@sbast3.sunysb.edu
c
c
c    call line:    call eoslog(inpvar,ye,brydns,eosflg)
c
c
c    inputs:       inpvar = temp, internal eng, or entropy
c                  ye = electron fraction
c                  brydns = baryon number density
c
c
c
c    outputs:      eosflg = 1 --> not implemented in model 4b
c                           2 --> general eos
c                           3 --> bulk eos (includes alpha's)
c
c
c
c 
c    include files:  eos_m4c.inc, maxwel.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine eoslog(inpvar,ye,brydns,eosflg)
c
c
      implicit none
c
c                       this include file contains all variable
c                       declarations.  note:: no implicit typing
c                       scheme is followed in this code; if you
c                       have any doubt as to a variables type check
c                       it!!!!.  also note that all variables are
c                       declared explicitly.
c
      include 'eos_m4c.inc'
      include 'maxwel.inc'
c
      double precision nlow, nhi, n_cut, temp_1, temp_2, t_bndy
c
      double precision lmm, lmp, lpm, lpp
      double precision dndy1, dndy2
c
c
c
 10   continue
c
c
c                         set t equal to the input variable (any calls
c                         with entropy or internal energy should go
c                         the the eos_m4b subroutine)
c
      t = inpvar(1)
c
c-----------------------------------------------------------------------
c         code to figure out the boundaries from the tables
c-----------------------------------------------------------------------
c
c
c
c
      if(ye.gt.y_hi) then
c                         ye is too large for eos
c
        write(*,*) ' eoslog:: cant do ye = ',ye, 'at this time'
        write(*,*) ' eoslog:: assuming ye =',y_hi,' instead'
        ye = y_hi-1.0d-6
        goto 10
c
      elseif(ye.ge.y_low) then
c                         calculate high and low boundary densities
c                         for the maxwell construction
c
c----------------------------------------------------------
c           calc ye index
c----------------------------------------------------------
c
        yfrac = (ye-y_low)/(y_hi-y_low)
        j_mxwl = int(yfrac*(numye-1))+1
        delt_y = (y_hi-y_low)/dble(numye-1)

c
c JCM:
c        if(j_mxwl>numye) then
c           j_mxwl=numye-1
c        end if
c JCM:
c        if(j_mxwl<1) then
c           j_mxwl=1
c        end if
c

c
        yminus = y_low+dble(j_mxwl-1)*delt_y
        yplus = y_low+dble(j_mxwl)*delt_y
        if((ye.ge.yminus).and.(ye.le.yplus)) then
          j_bd = j_mxwl
          j_bndy = j_mxwl
        elseif(ye.gt.yplus) then
          j_mxwl = j_mxwl+1
c JCM:
c          if(j_mxwl>numye) then
c             j_mxwl=numye-1
c          end if
          j_bd = j_mxwl
          j_bndy = j_mxwl
          yminus = y_low+dble(j_mxwl-1)*delt_y
          yplus = y_low+dble(j_mxwl)*delt_y
        else
          j_mxwl = j_mxwl-1
c JCM:
c          if(j_mxwl<1) then
c             j_mxwl=1
c          end if
          j_bd = j_mxwl
          j_bndy = j_mxwl
          yminus = y_low+dble(j_mxwl-1)*delt_y
          yplus = y_low+dble(j_mxwl)*delt_y
        endif
c
c
        if(j_mxwl.gt.(numye-1)) then
          j_mxwl = numye-1
          j_bd = j_mxwl
          j_bndy = j_mxwl
          yminus = y_low+dble(j_mxwl-1)*delt_y
          yplus = y_low+dble(j_mxwl)*delt_y
        endif
c
c
        yintrp = (ye-yminus)/(yplus-yminus)
c
c----------------------------------------------------------
c           calc t index
c----------------------------------------------------------

c
        tfrac = (t-t_low)/(t_hi-t_low)
        i_mxwl = int(tfrac*(numtmp-1))+1
        delt_t = (t_hi-t_low)/dble(numtmp-1)

c
c JCM:
c        if(i_mxwl>numtmp) then
c           i_mxwl=numtmp-1
c        end if
c JCM: (SUPERMARK)
        if(i_mxwl<1) then
           i_mxwl=1
c JCM:
c           write(*,*) 'stuck?'
c           ssflag=0
c           goto 999
        end if
c

c
        tminus = t_low+dble(i_mxwl-1)*delt_t
        tplus = t_low+dble(i_mxwl)*delt_t

        if((t.gt.tminus).and.(t.le.tplus)) then
          tminus = t_low+dble(i_mxwl-1)*delt_t
          tplus = t_low+dble(i_mxwl)*delt_t
        elseif(t.gt.tplus) then   
          i_mxwl = i_mxwl+1
c JCM:
c        if(i_mxwl>numtmp) then
c           i_mxwl=numtmp-1
c        end if
          tminus = t_low+dble(i_mxwl-1)*delt_t
          tplus = t_low+dble(i_mxwl)*delt_t
        else
          i_mxwl = i_mxwl-1
c JCM: (SUPERMARK)
c        if(i_mxwl<1) then
c           i_mxwl=1
c           write(*,*) '2stuck?'
c     Assume failed
           ssflag=0
           goto 999
c        end if
          tminus = t_low+dble(i_mxwl-1)*delt_t
          tplus = t_low+dble(i_mxwl)*delt_t
        endif
c
c
        if(i_mxwl.gt.(numtmp-1)) then
          i_mxwl = numtmp-1
          tminus = t_low+dble(i_mxwl-1)*delt_t
          tplus = t_low+dble(i_mxwl)*delt_t
        endif
c
c
        tintrp = (t-tminus)/(tplus-tminus)
c
c
c
c
c                find the temperature and density at the top of the
c                maxwel construction
c
cc      t_mxwl = yintrp*(t_h(j_mxwl+1)-t_h(j_mxwl))+t_h(j_mxwl)
cc      d_mxwl = yintrp*(d_h(j_mxwl+1)-d_h(j_mxwl))+d_h(j_mxwl)
        t_mxwl = dmin1(t_h(j_mxwl+1),t_h(j_mxwl))
        if(t_h(j_mxwl+1).gt.t_h(j_mxwl)) then
          d_mxwl = d_h(j_mxwl)
        else
          d_mxwl = d_h(j_mxwl+1)
        endif
c
c
c
c--------------------------------------------------------------------
c            interpolate to get maxwell construction densities
c--------------------------------------------------------------------
c
c
c
c JCM: DEBUG:
c        write(*,*) 'before dns_1'
c        write(*,*) i_mxwl,j_mxwl+1
c        write(*,*) i_mxwl,j_mxwl
c        write(*,*) numtmp,numye
        dns_1 = yintrp*(brylow(i_mxwl,j_mxwl+1)-brylow(i_mxwl,j_mxwl))+
     1               brylow(i_mxwl,j_mxwl)
c JCM: DEBUG:
c        write(*,*) 'before dns_2'
        dns_2 = yintrp*
     1        (brylow(i_mxwl+1,j_mxwl+1)-brylow(i_mxwl+1,j_mxwl))+
     2               brylow(i_mxwl+1,j_mxwl)
c
        lowdns = tintrp*(dns_2-dns_1)+dns_1
c
c                derivative of lower density w.r.t. t
        dnl_dt = (dns_2-dns_1)/delt_t
c
c JCM: DEBUG:
c        write(*,*) 'before dndy1'
        dndy1 = (brylow(i_mxwl,j_mxwl+1)-brylow(i_mxwl,j_mxwl))/delt_y
c JCM: DEBUG:
c        write(*,*) 'before dndy2'
        dndy2 = (brylow(i_mxwl+1,j_mxwl+1)-
     1      brylow(i_mxwl+1,j_mxwl))/delt_y
        dnl_dy = tintrp*(dndy2-dndy1)+dndy1
c
c
c
c
        if(ye.gt.y_cut) then
c
          dns_1 = yintrp*
     1        (bryhi(i_mxwl,j_mxwl+1)-bryhi(i_mxwl,j_mxwl))+
     2        bryhi(i_mxwl,j_mxwl)
          dns_2 = yintrp*
     1        (bryhi(i_mxwl+1,j_mxwl+1)-bryhi(i_mxwl+1,j_mxwl))+
     2               bryhi(i_mxwl+1,j_mxwl)
c
          hidns = tintrp*(dns_2-dns_1)+dns_1
c
c                derivative of higher density w.r.t. t
          dnh_dt = (dns_2-dns_1)/delt_t
c
c
        dndy1 = (bryhi(i_mxwl,j_mxwl+1)-
     1      bryhi(i_mxwl,j_mxwl))/delt_y
        dndy2 = (bryhi(i_mxwl+1,j_mxwl+1)-
     1      bryhi(i_mxwl+1,j_mxwl))/delt_y
        dnh_dy = tintrp*(dndy2-dndy1)+dndy1
c
c
        else
          hidns = lowdns
        endif
c
c
c--------------------------------------------------------------------
c--------------------------------------------------------------------
c
c                       ye is too low
      else
        write(*,*) ' eoslog:: cant do ye = ',ye, 'at this time'
        write(*,*) ' eoslog:: assuming ye =',y_low,' instead'
        ye = y_low+1.0d-6
        goto 10
      endif
c
c
c
c
c
      dltln1 = (lncut-lnlow)/dble(numlow-1)
      dltln2 = (lnhi-lncut)/dble(numhi-1)
c
c
      nlow = 10.0**lnlow
      nhi = 10.0**lnhi
      n_cut = 10.0**lncut
      logbry = dlog10(brydns)
      logbch = logbry
c
c
c
c
c----------------------------------------------------------
c           calc t index
c----------------------------------------------------------
c
c
      if(logbry.ge.lnhi) then
        i_bd = nbpnts
        i_bndy = nbpnts


c JCM:
        if( (i_bndy.gt.nbpnts).OR.(j_bndy+1.gt.numye) ) then
           write(*,*) "lbound requested out of bounds",i_bndy+1,nbpnts,j_bndy+1,numye
        end if

c JCM:
c        write(*,*) 'before t_bndy'
        t_bndy = yintrp*
     1           (lbound(i_bndy,j_bndy+1)-lbound(i_bndy,j_bndy))+
     2            lbound(i_bndy,j_bndy)
        tchk_b = 1.01*t_bndy
        tchk_n = 0.95*t_bndy
        goto 70
      elseif((logbry.lt.lnhi).and.(logbry.gt.lncut)) then
c
        i_bd = int((logbry-lncut)/dltln2)+numlow
        lnmins = lncut+dble(i_bd-numlow)*dltln2
        lnplus = lncut+dble(i_bd-numlow+1)*dltln2
        if((logbch.le.lnplus).and.(logbch.ge.lnmins)) then
          i_bndy = i_bd
        elseif(logbch.gt.lnplus) then
          i_bd = i_bd+1
          i_bndy = i_bd
          lnmins = lncut+dble(i_bndy-numlow)*dltln2
          lnplus = lncut+dble(i_bndy-numlow+1)*dltln2
        else
          i_bd = i_bd-1
          i_bndy = i_bd
c JCM:
          if(i_bndy<1) then
c     Assume failed
             ssflag=0
             goto 999
          end if
          lnmins = lncut+dble(i_bndy-numlow)*dltln2
          lnplus = lncut+dble(i_bndy-numlow+1)*dltln2
        endif
c
      elseif((logbry.le.lncut).and.(logbry.gt.lnlow)) then
c
        i_bd = int((logbry-lnlow)/dltln1)+1
        lnmins = lnlow+dble(i_bd-1)*dltln1
        lnplus = lnlow+dble(i_bd)*dltln1
        if((logbch.le.lnplus).and.(logbch.ge.lnmins)) then
          i_bndy = i_bd
        elseif(logbch.gt.lnplus) then
          i_bd = i_bd+1
          i_bndy = i_bd
          lnmins = lnlow+dble(i_bndy-1)*dltln1
          lnplus = lnlow+dble(i_bndy)*dltln1
        else
          i_bd = i_bd-1
          i_bndy = i_bd
c JCM:
          if(i_bndy<1) then
c     Assume failed
             ssflag=0
             goto 999
          end if
          lnmins = lnlow+dble(i_bndy-1)*dltln1
          lnplus = lnlow+dble(i_bndy)*dltln1
        endif
c
      else
c JCM: (hopefully bug fix will find this out)
c         write(*,*) "Never chose bound"
c    
      end if
c
      if(i_bndy.gt.(nbpnts-1)
c JCM:
     1   .OR.(i_bndy+1.gt.nbpnts).OR.(j_bndy+1.gt.numye)
     1     ) then
c JCM: hopefully this is caught
c         write(*,*) "bug-fix",i_bndy,nbpnts
        i_bd = nbpnts-1
        i_bndy = i_bd
        lnmins = lncut+dble(i_bndy-numlow)*dltln2
        lnplus = lncut+dble(i_bndy-numlow+1)*dltln2
      endif
c
c
c
c JCM:
      if( (i_bndy+1.gt.nbpnts).OR.(j_bndy+1.gt.numye) ) then
         write(*,*) "lbound requested out of bounds2",
     1               i_bndy+1,nbpnts,j_bndy+1,numye
         write(*,*) i_bd
      end if
c JCM:
c JCM:
      if(i_bndy<1) then
c     Assume failed
         ssflag=0
         goto 999
      end if
      if(j_bndy<1) then
c     Assume failed
         ssflag=0
         goto 999
      end if
      if(i_bndy+1>nbpnts) then
c     Assume failed
         ssflag=0
         goto 999
      end if
      if(j_bndy+1>numye) then
c     Assume failed
         ssflag=0
         goto 999
      end if
c      write(*,*) 'before lmm'
      lmm = lbound(i_bndy,j_bndy)
      lpm = lbound(i_bndy+1,j_bndy)
      lmp = lbound(i_bndy,j_bndy+1)
      lpp = lbound(i_bndy+1,j_bndy+1)
c
      lnfrac = (logbch-lnmins)/(lnplus-lnmins)
c
c                interpolate in ye first
c
c JCM:
c      write(*,*) 'before temp_1',i_bndy,j_bndy
      temp_1 = yintrp*
     1           (lbound(i_bndy,j_bndy+1)-lbound(i_bndy,j_bndy))+
     2            lbound(i_bndy,j_bndy)
      temp_2 = yintrp*
     1        (lbound(i_bndy+1,j_bndy+1)-lbound(i_bndy+1,j_bndy))+
     2               lbound(i_bndy+1,j_bndy)
c
c                interpolate in density between the two ye
c                interpolated values
c
      t_bndy = lnfrac*(temp_2-temp_1)+temp_1
c
      tchk_b = 1.01*t_bndy
      tchk_n = 0.95*t_bndy
c
      if((lmm.ge.lpm).or.(lmp.gt.lpp)) then
        tchk_n = dmax1(0.0d0,dmin1(0.95*tchk_n,t_bndy-3.0))
      endif
c
 70   continue
c
c----------------------------------------------------------
c----------------------------------------------------------
c
c
c
c-----------------------------------------------------------------------
c               eos logic
c-----------------------------------------------------------------------
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c                     if t is below the maximum for maxwel construction
      if(t.lt.t_mxwl) then
c                       if rho is greater than the upper max. con.
c                       density the use the bulk eos
        if(brydns.gt.hidns) then
          eosflg = 3
c                       else if rho is greater than the lower max. con.
c                       density then
        elseif(brydns.gt.lowdns) then
c                         if ye is large enough to have a signifigant
c                         max con then use the maxwell con. eos
          if(ye.gt.y_cut) then
            eosflg = 4
c                         otherwise use the bulk eos
          else
            eosflg = 3
          endif
c
c                       if density is greater than the minimum
c                       maxwell con. density, then we know that we are
c                       in the nuclear eos density
        elseif(brydns.gt.d_mxwl) then
          eosflg = 2
c
c
c                       otherwise check the boundary table
        else
c
c                         if t is well below the phase boundary curve
c                         then use the nuclear eos
          if(t.lt.tchk_n) then
            eosflg = 2
c                         otherwise if t is near the boundary, first
c                         try the nuclear eos and if not successfull
c                         then use the bulk eos
          elseif(t.lt.tchk_b) then
            eosflg = 1
          else
c                         otherwise t is well above the boundary so
c                         use the bulk eos
            eosflg = 3
          endif
        endif
c
c                     otherwise t is above the maximum for a maxwell
c                     construction
      else
c                       if density is greater than that at the top of
c                       the maxwell construction then use the bulk eos
        if(brydns.gt.d_mxwl) then
          eosflg = 3
c
c                       otherwise density is below the maxwell con.
        else
c
c                         if t is well below the phase boundary curve
c                         then use the nuclear eos
          if(t.lt.tchk_n) then
            eosflg = 2
c
c                         otherwise if t is near the phase boundary
c                         curve then try the nuclear eos and if not
c                         successfull then use the bulk eos
          elseif(t.lt.tchk_b) then
            eosflg = 1
c
c                         otherwise t is well above the phase boundary
c                         curve so use the bulk eos
          else
            eosflg = 3
          endif
        endif
      endif  
c
c
c-----------------------------------------------------------------------
c                         done with eos logic so return eosflg
c-----------------------------------------------------------------------
c
 999  return
c
c
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         reset.for
c
c***********************************************************************
c
c    module:       reset
c    type:         subroutine
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         12/21/90
c
c                  please report any problems to me at:
c                  bitnet:  swesty@sunysbnp or
c                  internet: fswesty@astro.sunysb.edu or
c                            fswesty@sbast3.sunysb.edu
c
c
c    call line:    call reset(inpvar,ye,brydns,outvar)
c
c
c    inputs:       inpvar = temp, nsubi, eta_po, eta_no
c                  ye = electron fraction
c                  brydns = baryon number density
c
c
c
c    outputs:      outvar = array of length 4 containing reset values
c                  for the initial guesses
c
c
c
c 
c    include files: none
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine reset(inpvar,ye,brydns,outvar)
c
c
      implicit none
c
c
c                      subroutine parameters
c
      double precision inpvar(4), outvar(4), ye, brydns
c
c
c                      local variables
c
      double precision zpg, zng, eta_pg, eta_ng, pi, uq, mq, t, efrac
c
c                      functions
c
      double precision finv12
c
c-----------------------------------------------------------------------
c
      t = inpvar(1)
c
c
      pi = 3.1415927
      uq = 20.721
c
      mq = (t/uq)**1.5
c
c
      efrac = 0.5*ye
c
      zng = 2.0*(pi**2)*brydns*(1.0-efrac)/mq
c
      zpg = 2.0*(pi**2)*brydns*efrac/mq
c
      eta_ng = finv12(zng)
c
      eta_pg = finv12(zpg)
c
      outvar(1) = inpvar(1)
      outvar(2) = inpvar(2)
      outvar(3) = eta_pg
      outvar(4) = eta_ng
c
c
c-----------------------------------------------------------------------
c
 999  return
c
c
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         aloadmx.for
c    module:       loadmx
c    type:         loadmx
c
c    purpose:      load the look-up table for the maxwell construction
c
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    date:         7/16/90
c
c    call line:    call loadmx
c
c    inputs:       n/a
c
c    outputs       n/a
c
c    subroutine calls: eos_m4c
c 
c    include files:  eos_m4c.inc, maxwel.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine loadmx
c
      implicit none
c
c
      include 'eos_m4c.inc'
      include 'maxwel.inc'
c
      integer ntmp, nye, nye2, num_bp
      integer lun1, lun2, kk, kmin
      parameter(lun1=54,lun2=55)
c
c
      integer fnml1, fnml2
      character*60 fname1, fname2
c
      double precision n_sm, symm_m, comp_m, bindem, sym_sm, sig_sm
      double precision n_sb, symm_b, comp_b, bindeb, sym_sb, sig_sb
      double precision mscdd3, bscdd3
c
      include 'force.inc'
      include 'const.dek'
c
c
c
cc      call getfnm(fname1,fnml1,'enter ascii maxwell fname:',26)
c JCM:
c      fname1 = 'max220.atb'
      fname1 = 'maxwel.atb'
      fnml1=10
c
cc      call getfnm(fname2,fnml2,'enter ascii boundary fname:',27)
c JCM:
c      fname2 = 'bd220.atb'
      fname2 = 'bound.atb'
      fnml2=9
c
c
c
c-----------------------------------------------------------------------
c        read the file maxwell construction data file
c-----------------------------------------------------------------------
c
c
c
c      write(*,*) "Open maxwel.atb"
      open(unit=lun1,file=fname1(1:fnml1),status='old')
c
c
c
c
c
c
      read(lun1,*) n_sm, symm_m
      read(lun1,*) comp_m,bindem
      read(lun1,*) sym_sm, sig_sm
c
c
c
c
      read(lun1,*) ntmp,nye
      read(lun1,*) t_low,t_hi
      read(lun1,*) y_low,y_hi
c
c
c
      if((ntmp.ne.numtmp).or.(nye.ne.numye)) then
        write(*,*) 'loadmx:  mxwl table is incompatible with arrays'
        stop
      endif
c
c
      do 101 j=1,numye,1
        do 100 i=1,numtmp,3
c JCM:
c           write(*,*) 'brylow',i,j
          kmin = min0(i+2,numtmp)
          
c JCM:
c          if(kk<1 .OR. kk>numtmp) then
c             write(*,*) 'brylow out of bounds tmp',kk,numtmp
c             write(*,*) 'kmin=',kmin
c          end if
c          if(j<1 .OR. j>numye) then
c             write(*,*) 'brylow out of bounds ye',j,numye
c             write(*,*) 'kmin=',kmin
c          end if
c          write(*,*) 'before read'
          read(lun1,*) (brylow(kk,j),kk=i,kmin,1)
 100    continue
 101  continue
c
c
      do 103 j=1,numye,1
        do 102 i=1,numtmp,3
c JCM:
c           write(*,*) 'bryye',i,j
          kmin = min0(i+2,numtmp)
          read(lun1,*) (bryhi(kk,j),kk=i,kmin,1)
 102    continue
 103  continue
c
c
c
      do 104 i=1,numye,3
c JCM:
c           write(*,*) 'bryye2',i
        kmin = min0(i+2,numye)
        read(lun1,*) (t_h(kk),kk=i,kmin,1)
 104  continue
c
c
      do 105 i=1,numye,3
c JCM:
c           write(*,*) 'bryye3',i
        kmin = min0(i+2,numye)
        read(lun1,*) (d_h(kk),kk=i,kmin,1)
 105  continue
c
      read(lun1,*) ycut
      read(lun1,*) mscdd3
c
c
      close(unit=lun1,status='keep')
c      write(*,*) "Close maxwel.atb"

c
c
c..      write(*,*)
c..      write(*,*) '<<loadmx:  maxwell con. table is initialized>>'
c..      write(*,*)
c
c
c
c-----------------------------------------------------------------------
c        read the file boundary data file
c-----------------------------------------------------------------------
c
c
c
c      write(*,*) "Open bound.atb"
      open(unit=lun2,file=fname2(1:fnml2),status='old')
c
c
c
c
      read(lun2,*) n_sb,symm_b
      read(lun2,*) comp_b,bindeb
      read(lun2,*) sym_sb,sig_sb
c
c
c
c
c
c
      read(lun2,*) num_bp,nye2
      read(lun2,*) lnl,lnh,lnc
      read(lun2,*) y_low2,y_hi2
c
c
      if((nbpnts.ne.num_bp).or.(nye2.ne.numye)) then
        write(*,*) 'loadmx:  bndy table is incompatible with arrays'
        stop
      endif
c
      if(abs(lnl-lnlow).gt.1.0d-10) then
        write(*,*) 'loadmx:  lower end of phase bndy is inconsist.'
        stop
      endif
c
c
      if(abs(lnh-lnhi).gt.1.0d-10) then
        write(*,*) 'loadmx:  upper end of phase bndy is inconsist.'
        stop
      endif
c
c
      if(abs(lnc-lncut).gt.1.0d-10) then
        write(*,*) 'loadmx:  mid cut of phase bndy is inconsist.'
        stop
      endif
c
      if(abs(y_low-y_low2).gt.1.0d-10) then
        write(*,*) 'loadmx:  lower ye limits are inconsist.'
        stop
      endif
c
      if(abs(y_hi-y_hi2).gt.1.0d-10) then
        write(*,*) 'loadmx:  upper ye limits are inconsist.'
        stop
      endif
c
c
      do 201 j=1,numye,1
        do 200 i=1,nbpnts,3
c JCM:
c           write(*,*) 'bryye4',i,j
          kmin = min0(i+2,nbpnts)
          read(lun2,*) (lbound(kk,j),kk=i,kmin,1)
 200    continue
 201  continue
c
c
      do 203 j=1,numye,1
        do 202 i=1,nbpnts,3
c JCM:
c           write(*,*) 'bryye5',i,j
          kmin = min0(i+2,nbpnts)
          read(lun2,*) (ubound(kk,j),kk=i,kmin,1)
 202    continue
 203  continue
c
      read(lun2,*) bscdd3
c
      if(abs(mscdd3-bscdd3).gt.1.0d-10) then
        write(*,*) 'loadmx:  scrdd3 values are inconsist.'
        stop
      endif
c
c
c
      close(unit=lun2,status='keep')
c      write(*,*) "Close bound.atb"

c
c
c..      write(*,*)
c..      write(*,*) '<<loadmx:  boundary table is initialized>>'
c..      write(*,*)
c
c
c
c-----------------------------------------------------------------------
c                  all arrays are now loaded so return
c-----------------------------------------------------------------------
c
      scrdd3 = bscdd3
      n_s = n_sm
      nsubs = n_sm
      symm = symm_m
      comp = comp_m
      bind_e = bindem
      sym_s = sym_sm
      sig_s = sig_sm
c
c20      skyrmc=(.3*((hbar*c)**2)/massn)*(1.5*n_s*(pi**2))**ovr23
c20      dd = (comp+2.0*skyrmc)/(3.0*skyrmc+9.0*bind_e)
c20      bb = (skyrmc*(2.0**ovr23-1.0)-symm)/n_s
c20      aa = (ovr23*skyrmc-dd*(skyrmc+bind_e))/(n_s*(dd-1.0))-bb
c20      cc = (comp+2.0*skyrmc)/(9.0*dd*(dd-1.0)*n_s**dd)
      skyrmc=(.3*((hbar*c)**2)/massn)*(1.5*n_s*(pi**2))**ovr23
      dd = (comp+2.0*skyrmc+scrdd3*(comp-4.0*skyrmc-18.0*bind_e))/
     1    ((1.0-scrdd3)*(3.0*skyrmc+9.0*bind_e))
      bb = (skyrmc*(2.0**ovr23-1.0)-symm)/n_s
      cc = ((ovr3*skyrmc+bind_e)*(1.0+scrdd3)**2)/((n_s**dd)*(dd-1.0))
      aa = ((ovr23*skyrmc-dd*(skyrmc+bind_e)-scrdd3*(ovr3*skyrmc+bind_e)
     1    )/(n_s*(dd-1.0)) )-bb
      dd3 = scrdd3/(n_s**(dd-1.0))
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c..      write(*,*)
c..      write(*,*) '<<loadmx:  skyrme parameters for this run are:>>'
c..      write(*,*) 'abcd: ',aa,bb,cc,dd,scrdd3
c..      write(*,*) ' satur. density, symmetry engy, & compression mod.:'
c..      write(*,*) n_sm, symm_m, comp_m
c..      write(*,*) n_sb, symm_b, comp_b
c..      write(*,*) ' binding engy, surf. symm. engy, & surface tension:'
c..      write(*,*) bindem,sym_sm,sig_sm
c..      write(*,*) bindeb,sym_sb,sig_sb
c..c
c..      write(*,*)
c
c
c
c..      write(*,*)
c..      write(*,*) '<<loadmx: fermi integral tables are initialized>>'
c..      write(*,*)
c
c
 999  return
c
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       getfnm.for
c    type:         subroutine
c    author:       f. douglas swesty
c    date:         8/5/89
c
c    purpose:      obtains a file name from user
c
c    call line:    call getfnm(fname,fnml,prompt,promtl)
c
c    inputs:       prompt = string to pompt user with (c*60)
c                  promtl = length of prompt string (i)
c
c    outputs:      fname = file name (c*60)
c                  fnml = file name length (i)
c*************************************************************************
c
      subroutine getfnm(fname,fnml,prompt,promtl)
c
      implicit none
c
      integer fnml, promtl
      character*60 fname, prompt
c
c                       local variables
c
      integer stdout, stdin, i
      data stdout/6/, stdin/5/
c
c                       prompt user for file name
c
      write(stdout,'(t2,a,$)') prompt(1:promtl)
      read(stdin,'(a)') fname
c
c                        figure out input file name length
      do 10 i=1,20,1
c JCM:
c           write(*,*) 'file',i
        if(fname(i:i).eq.' ') goto 20
 10   continue
c
 20   continue
      fnml = i-1
c
c
 999  return
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       f_1_2
c    type:         double precision function
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    call line:    f_1_2(y)      (1/2th fermi integral)
c
c    inputs:       y (double precision)   (argument)
c
c    return:       1/2th fermi integral (double precision)
c
c***********************************************************************
      double precision function f_1_2(y)
      implicit none
c                     
      integer nlow,nhigh,n
      parameter(n=201)
      double precision y, a(7), eta(n), f12(n), f12a(n), th
      double precision f0, f1, x2
c
c                retain values between calls
      save nlow,nhigh
c
      data a/6.16850274d0,1.77568655d0,6.92965606d0,
     $   .176776695d0,6.41500299d-02,.4d0,1.32934039d0/
      data th,nlow,nhigh/0.33333333333d0,1,n/
c
c
c                    cubic spline data
      data eta(  1),f12(  1),f12a(  1)
     $/-1.00000000e+01, 4.02339983e-05, 0.000000000000e+00/
      data eta(  2),f12(  2),f12a(  2)
     $/-9.80000000e+00, 4.91417786e-05, 5.971803206329e-05/
      data eta(  3),f12(  3),f12a(  3)
     $/-9.60000000e+00, 6.00216308e-05, 5.693865674682e-05/
      data eta(  4),f12(  4),f12a(  4)
     $/-9.40000000e+00, 7.33101788e-05, 7.383171094946e-05/
      data eta(  5),f12(  5),f12a(  5)
     $/-9.20000000e+00, 8.95406574e-05, 8.902408945532e-05/
      data eta(  6),f12(  6),f12a(  6)
     $/-9.00000000e+00, 1.09364317e-04, 1.090490812293e-04/
      data eta(  7),f12(  7),f12a(  7)
     $/-8.80000000e+00, 1.33576705e-04, 1.330888456275e-04/
      data eta(  8),f12(  8),f12a(  8)
     $/-8.60000000e+00, 1.63148990e-04, 1.625800862606e-04/
      data eta(  9),f12(  9),f12a(  9)
     $/-8.40000000e+00, 1.99267717e-04, 1.985571093302e-04/
      data eta( 10),f12( 10),f12a( 10)
     $/-8.20000000e+00, 2.43381819e-04, 2.424977264186e-04/
      data eta( 11),f12( 11),f12a( 11)
     $/-8.00000000e+00, 2.97260769e-04, 2.961791849954e-04/
      data eta( 12),f12( 12),f12a( 12)
     $/-7.80000000e+00, 3.63065723e-04, 3.616861335996e-04/
      data eta( 13),f12( 13),f12a( 13)
     $/-7.60000000e+00, 4.43435136e-04, 4.417451306063e-04/
      data eta( 14),f12( 14),f12a( 14)
     $/-7.40000000e+00, 5.41591857e-04, 5.394295439752e-04/
      data eta( 15),f12( 15),f12a( 15)
     $/-7.20000000e+00, 6.61470011e-04, 6.587516434928e-04/
      data eta( 16),f12( 16),f12a( 16)
     $/-7.00000000e+00, 8.07873978e-04, 8.044358320536e-04/
      data eta( 17),f12( 17),f12a( 17)
     $/-6.80000000e+00, 9.86669358e-04, 9.822169782926e-04/
      data eta( 18),f12( 18),f12a( 18)
     $/-6.60000000e+00, 1.20501561e-03, 1.199327054776e-03/
      data eta( 19),f12( 19),f12a( 19)
     $/-6.40000000e+00, 1.47165314e-03, 1.464166502604e-03/
      data eta( 20),f12( 20),f12a( 20)
     $/-6.20000000e+00, 1.79724747e-03, 1.787526934808e-03/
      data eta( 21),f12( 21),f12a( 21)
     $/-6.00000000e+00, 2.19481561e-03, 2.181797258165e-03/
      data eta( 22),f12( 22),f12a( 22)
     $/-5.80000000e+00, 2.68023420e-03, 2.662851532532e-03/
      data eta( 23),f12( 23),f12a( 23)
     $/-5.60000000e+00, 3.27287085e-03, 3.249505611708e-03/
      data eta( 24),f12( 24),f12a( 24)
     $/-5.40000000e+00, 3.99634114e-03, 3.964172020638e-03/
      data eta( 25),f12( 25),f12a( 25)
     $/-5.20000000e+00, 4.87942214e-03, 4.835412805738e-03/
      data eta( 26),f12( 26),f12a( 26)
     $/-5.00000000e+00, 5.95717982e-03, 5.895678756409e-03/
      data eta( 27),f12( 27),f12a( 27)
     $/-4.80000000e+00, 7.27229903e-03, 7.186101668628e-03/
      data eta( 28),f12( 28),f12a( 28)
     $/-4.60000000e+00, 8.87672161e-03, 8.755420069081e-03/
      data eta( 29),f12( 29),f12a( 29)
     $/-4.40000000e+00, 1.08335980e-02, 1.066028955505e-02/
      data eta( 30),f12( 30),f12a( 30)
     $/-4.20000000e+00, 1.32195998e-02, 1.297223321072e-02/
      data eta( 31),f12( 31),f12a( 31)
     $/-4.00000000e+00, 1.61277414e-02, 1.577174760206e-02/
      data eta( 32),f12( 32),f12a( 32)
     $/-3.80000000e+00, 1.96706583e-02, 1.915707138105e-02/
      data eta( 33),f12( 33),f12a( 33)
     $/-3.60000000e+00, 2.39845128e-02, 2.324060687373e-02/
      data eta( 34),f12( 34),f12a( 34)
     $/-3.40000000e+00, 2.92335232e-02, 2.815388612401e-02/
      data eta( 35),f12( 35),f12a( 35)
     $/-3.20000000e+00, 3.56152142e-02, 3.404593863024e-02/
      data eta( 36),f12( 36),f12a( 36)
     $/-3.00000000e+00, 4.33663810e-02, 4.108372935503e-02/
      data eta( 37),f12( 37),f12a( 37)
     $/-2.80000000e+00, 5.27697736e-02, 4.945301394966e-02/
      data eta( 38),f12( 38),f12a( 38)
     $/-2.60000000e+00, 6.41614366e-02, 5.934477484635e-02/
      data eta( 39),f12( 39),f12a( 39)
     $/-2.40000000e+00, 7.79383797e-02, 7.095990166492e-02/
      data eta( 40),f12( 40),f12a( 40)
     $/-2.20000000e+00, 9.45664588e-02, 8.448601849395e-02/
      data eta( 41),f12( 41),f12a( 41)
     $/-2.00000000e+00, 1.14587830e-01, 1.000898393593e-01/
      data eta( 42),f12( 42),f12a( 42)
     $/-1.80000000e+00, 1.38627354e-01, 1.178775440690e-01/
      data eta( 43),f12( 43),f12a( 43)
     $/-1.60000000e+00, 1.67396817e-01, 1.378908343646e-01/
      data eta( 44),f12( 44),f12a( 44)
     $/-1.40000000e+00, 2.01696221e-01, 1.600502684727e-01/
      data eta( 45),f12( 45),f12a( 45)
     $/-1.20000000e+00, 2.42410529e-01, 1.841436917447e-01/
      data eta( 46),f12( 46),f12a( 46)
     $/-1.00000000e+00, 2.90500917e-01, 2.097869645484e-01/
      data eta( 47),f12( 47),f12a( 47)
     $/-8.00000000e-01, 3.46989460e-01, 2.364317000617e-01/
      data eta( 48),f12( 48),f12a( 48)
     $/-6.00000000e-01, 4.12937023e-01, 2.633392352050e-01/
      data eta( 49),f12( 49),f12a( 49)
     $/-4.00000000e-01, 4.89414580e-01, 2.897104591184e-01/
      data eta( 50),f12( 50),f12a( 50)
     $/-2.00000000e-01, 5.77470496e-01, 3.145727783213e-01/
      data eta( 51),f12( 51),f12a( 51)
     $/ 0.00000000e+00, 6.78093925e-01, 3.371253775963e-01/
      data eta( 52),f12( 52),f12a( 52)
     $/ 2.00000000e-01, 7.92181447e-01, 3.565396612936e-01/
      data eta( 53),f12( 53),f12a( 53)
     $/ 4.00000000e-01, 9.20506015e-01, 3.722728772295e-01/
      data eta( 54),f12( 54),f12a( 54)
     $/ 6.00000000e-01, 1.06369475e+00, 3.839938797886e-01/
      data eta( 55),f12( 55),f12a( 55)
     $/ 8.00000000e-01, 1.22221592e+00, 3.916168536161e-01/
      data eta( 56),f12( 56),f12a( 56)
     $/ 1.00000000e+00, 1.39637545e+00, 3.952927057472e-01/
      data eta( 57),f12( 57),f12a( 57)
     $/ 1.20000000e+00, 1.58632329e+00, 3.954588233952e-01/
      data eta( 58),f12( 58),f12a( 58)
     $/ 1.40000000e+00, 1.79206851e+00, 3.924790006721e-01/
      data eta( 59),f12( 59),f12a( 59)
     $/ 1.60000000e+00, 2.01349622e+00, 3.869986739163e-01/
      data eta( 60),f12( 60),f12a( 60)
     $/ 1.80000000e+00, 2.25039083e+00, 3.795613036627e-01/
      data eta( 61),f12( 61),f12a( 61)
     $/ 2.00000000e+00, 2.50245792e+00, 3.706281114328e-01/
      data eta( 62),f12( 62),f12a( 62)
     $/ 2.20000000e+00, 2.76934439e+00, 3.608332506061e-01/
      data eta( 63),f12( 63),f12a( 63)
     $/ 2.40000000e+00, 3.05065972e+00, 3.503678861428e-01/
      data eta( 64),f12( 64),f12a( 64)
     $/ 2.60000000e+00, 3.34598833e+00, 3.396872048225e-01/
      data eta( 65),f12( 65),f12a( 65)
     $/ 2.80000000e+00, 3.65490490e+00, 3.290772945673e-01/
      data eta( 66),f12( 66),f12a( 66)
     $/ 3.00000000e+00, 3.97698528e+00, 3.185751169082e-01/
      data eta( 67),f12( 67),f12a( 67)
     $/ 3.20000000e+00, 4.31181109e+00, 3.084367378000e-01/
      data eta( 68),f12( 68),f12a( 68)
     $/ 3.40000000e+00, 4.65897715e+00, 2.987154318921e-01/
      data eta( 69),f12( 69),f12a( 69)
     $/ 3.60000000e+00, 5.01809514e+00, 2.894910346313e-01/
      data eta( 70),f12( 70),f12a( 70)
     $/ 3.80000000e+00, 5.38879550e+00, 2.806759295830e-01/
      data eta( 71),f12( 71),f12a( 71)
     $/ 4.00000000e+00, 5.77072680e+00, 2.724462470369e-01/
      data eta( 72),f12( 72),f12a( 72)
     $/ 4.20000000e+00, 6.16355908e+00, 2.646860822693e-01/
      data eta( 73),f12( 73),f12a( 73)
     $/ 4.40000000e+00, 6.56698239e+00, 2.574639238859e-01/
      data eta( 74),f12( 74),f12a( 74)
     $/ 4.60000000e+00, 6.98070586e+00, 2.504822221872e-01/
      data eta( 75),f12( 75),f12a( 75)
     $/ 4.80000000e+00, 7.40445435e+00, 2.443601873650e-01/
      data eta( 76),f12( 76),f12a( 76)
     $/ 5.00000000e+00, 7.83797658e+00, 2.381380283526e-01/
      data eta( 77),f12( 77),f12a( 77)
     $/ 5.20000000e+00, 8.28102899e+00, 2.326146992243e-01/
      data eta( 78),f12( 78),f12a( 78)
     $/ 5.40000000e+00, 8.73338914e+00, 2.275641747505e-01/
      data eta( 79),f12( 79),f12a( 79)
     $/ 5.60000000e+00, 9.19485021e+00, 2.222666017741e-01/
      data eta( 80),f12( 80),f12a( 80)
     $/ 5.80000000e+00, 9.66520906e+00, 2.180364181530e-01/
      data eta( 81),f12( 81),f12a( 81)
     $/ 6.00000000e+00, 1.01442864e+01, 2.133612256141e-01/
      data eta( 82),f12( 82),f12a( 82)
     $/ 6.20000000e+00, 1.06319029e+01, 2.093926793908e-01/
      data eta( 83),f12( 83),f12a( 83)
     $/ 6.40000000e+00, 1.11278965e+01, 2.056330568230e-01/
      data eta( 84),f12( 84),f12a( 84)
     $/ 6.60000000e+00, 1.16321146e+01, 2.017500933175e-01/
      data eta( 85),f12( 85),f12a( 85)
     $/ 6.80000000e+00, 1.21444066e+01, 1.984515699067e-01/
      data eta( 86),f12( 86),f12a( 86)
     $/ 7.00000000e+00, 1.26646369e+01, 1.951886270556e-01/
      data eta( 87),f12( 87),f12a( 87)
     $/ 7.20000000e+00, 1.31926749e+01, 1.919489218711e-01/
      data eta( 88),f12( 88),f12a( 88)
     $/ 7.40000000e+00, 1.37283938e+01, 1.891506854598e-01/
      data eta( 89),f12( 89),f12a( 89)
     $/ 7.60000000e+00, 1.42716773e+01, 1.861383362905e-01/
      data eta( 90),f12( 90),f12a( 90)
     $/ 7.80000000e+00, 1.48224099e+01, 1.836609693776e-01/
      data eta( 91),f12( 91),f12a( 91)
     $/ 8.00000000e+00, 1.53804867e+01, 1.808477861992e-01/
      data eta( 92),f12( 92),f12a( 92)
     $/ 8.20000000e+00, 1.59458020e+01, 1.787228858260e-01/
      data eta( 93),f12( 93),f12a( 93)
     $/ 8.40000000e+00, 1.65182614e+01, 1.758756704959e-01/
      data eta( 94),f12( 94),f12a( 94)
     $/ 8.60000000e+00, 1.70977635e+01, 1.741794321910e-01/
      data eta( 95),f12( 95),f12a( 95)
     $/ 8.80000000e+00, 1.76842275e+01, 1.716916007395e-01/
      data eta( 96),f12( 96),f12a( 96)
     $/ 9.00000000e+00, 1.82775617e+01, 1.695841648522e-01/
      data eta( 97),f12( 97),f12a( 97)
     $/ 9.20000000e+00, 1.88776803e+01, 1.676317398519e-01/
      data eta( 98),f12( 98),f12a( 98)
     $/ 9.40000000e+00, 1.94845052e+01, 1.658338757391e-01/
      data eta( 99),f12( 99),f12a( 99)
     $/ 9.60000000e+00, 2.00979619e+01, 1.638027571928e-01/
      data eta(100),f12(100),f12a(100)
     $/ 9.80000000e+00, 2.07179742e+01, 1.622950954886e-01/
      data eta(101),f12(101),f12a(101)
     $/ 1.00000000e+01, 2.13444734e+01, 1.600518608534e-01/
      data eta(102),f12(102),f12a(102)
     $/ 1.02000000e+01, 2.19773812e+01, 1.587874610986e-01/
      data eta(103),f12(103),f12a(103)
     $/ 1.04000000e+01, 2.26166368e+01, 1.569682947511e-01/
      data eta(104),f12(104),f12a(104)
     $/ 1.06000000e+01, 2.32621732e+01, 1.554593598982e-01/
      data eta(105),f12(105),f12a(105)
     $/ 1.08000000e+01, 2.39139295e+01, 1.541792656549e-01/
      data eta(106),f12(106),f12a(106)
     $/ 1.10000000e+01, 2.45718484e+01, 1.522135774835e-01/
      data eta(107),f12(107),f12a(107)
     $/ 1.12000000e+01, 2.52358613e+01, 1.510664244108e-01/
      data eta(108),f12(108),f12a(108)
     $/ 1.14000000e+01, 2.59059148e+01, 1.496107248727e-01/
      data eta(109),f12(109),f12a(109)
     $/ 1.16000000e+01, 2.65819535e+01, 1.482706760991e-01/
      data eta(110),f12(110),f12a(110)
     $/ 1.18000000e+01, 2.72639241e+01, 1.470915707300e-01/
      data eta(111),f12(111),f12a(111)
     $/ 1.20000000e+01, 2.79517770e+01, 1.457080409816e-01/
      data eta(112),f12(112),f12a(112)
     $/ 1.22000000e+01, 2.86454568e+01, 1.441112653439e-01/
      data eta(113),f12(113),f12a(113)
     $/ 1.24000000e+01, 2.93449082e+01, 1.435868976411e-01/
      data eta(114),f12(114),f12a(114)
     $/ 1.26000000e+01, 3.00500951e+01, 1.418661440931e-01/
      data eta(115),f12(115),f12a(115)
     $/ 1.28000000e+01, 3.07609620e+01, 1.409485259852e-01/
      data eta(116),f12(116),f12a(116)
     $/ 1.30000000e+01, 3.14774652e+01, 1.397847519674e-01/
      data eta(117),f12(117),f12a(117)
     $/ 1.32000000e+01, 3.21995592e+01, 1.385324661454e-01/
      data eta(118),f12(118),f12a(118)
     $/ 1.34000000e+01, 3.29271975e+01, 1.377303834490e-01/
      data eta(119),f12(119),f12a(119)
     $/ 1.36000000e+01, 3.36603403e+01, 1.362210000606e-01/
      data eta(120),f12(120),f12a(120)
     $/ 1.38000000e+01, 3.43989420e+01, 1.362206163067e-01/
      data eta(121),f12(121),f12a(121)
     $/ 1.40000000e+01, 3.51429758e+01, 1.337115347148e-01/
      data eta(122),f12(122),f12a(122)
     $/ 1.42000000e+01, 3.58923769e+01, 1.340282448335e-01/
      data eta(123),f12(123),f12a(123)
     $/ 1.44000000e+01, 3.66471262e+01, 1.324054859504e-01/
      data eta(124),f12(124),f12a(124)
     $/ 1.46000000e+01, 3.74071779e+01, 1.317098113657e-01/
      data eta(125),f12(125),f12a(125)
     $/ 1.48000000e+01, 3.81724978e+01, 1.309852685862e-01/
      data eta(126),f12(126),f12a(126)
     $/ 1.50000000e+01, 3.89430513e+01, 1.293891142900e-01/
      data eta(127),f12(127),f12a(127)
     $/ 1.52000000e+01, 3.97187891e+01, 1.291032742542e-01/
      data eta(128),f12(128),f12a(128)
     $/ 1.54000000e+01, 4.04996882e+01, 1.283927886917e-01/
      data eta(129),f12(129),f12a(129)
     $/ 1.56000000e+01, 4.12857180e+01, 1.269305709802e-01/
      data eta(130),f12(130),f12a(130)
     $/ 1.58000000e+01, 4.20768328e+01, 1.266349273854e-01/
      data eta(131),f12(131),f12a(131)
     $/ 1.60000000e+01, 4.28730059e+01, 1.252747194805e-01/
      data eta(132),f12(132),f12a(132)
     $/ 1.62000000e+01, 4.36741991e+01, 1.252811946920e-01/
      data eta(133),f12(133),f12a(133)
     $/ 1.64000000e+01, 4.44803934e+01, 1.237655017505e-01/
      data eta(134),f12(134),f12a(134)
     $/ 1.66000000e+01, 4.52915468e+01, 1.235217983048e-01/
      data eta(135),f12(135),f12a(135)
     $/ 1.68000000e+01, 4.61076365e+01, 1.225923050320e-01/
      data eta(136),f12(136),f12a(136)
     $/ 1.70000000e+01, 4.69286280e+01, 1.213789815681e-01/
      data eta(137),f12(137),f12a(137)
     $/ 1.72000000e+01, 4.77544832e+01, 1.214467686948e-01/
      data eta(138),f12(138),f12a(138)
     $/ 1.74000000e+01, 4.85851870e+01, 1.201239436519e-01/
      data eta(139),f12(139),f12a(139)
     $/ 1.76000000e+01, 4.94207048e+01, 1.201574566968e-01/
      data eta(140),f12(140),f12a(140)
     $/ 1.78000000e+01, 5.02610178e+01, 1.185262295614e-01/
      data eta(141),f12(141),f12a(141)
     $/ 1.80000000e+01, 5.11060801e+01, 1.181326250585e-01/
      data eta(142),f12(142),f12a(142)
     $/ 1.82000000e+01, 5.19558725e+01, 1.184582702039e-01/
      data eta(143),f12(143),f12a(143)
     $/ 1.84000000e+01, 5.28103838e+01, 1.158692941259e-01/
      data eta(144),f12(144),f12a(144)
     $/ 1.86000000e+01, 5.36695566e+01, 1.172895532908e-01/
      data eta(145),f12(145),f12a(145)
     $/ 1.88000000e+01, 5.45334024e+01, 1.159224927128e-01/
      data eta(146),f12(146),f12a(146)
     $/ 1.90000000e+01, 5.54018793e+01, 1.136854758577e-01/
      data eta(147),f12(147),f12a(147)
     $/ 1.92000000e+01, 5.62749338e+01, 1.159756038573e-01/
      data eta(148),f12(148),f12a(148)
     $/ 1.94000000e+01, 5.71525927e+01, 1.130721087129e-01/
      data eta(149),f12(149),f12a(149)
     $/ 1.96000000e+01, 5.80347986e+01, 1.137859612890e-01/
      data eta(150),f12(150),f12a(150)
     $/ 1.98000000e+01, 5.89215441e+01, 1.127240461327e-01/
      data eta(151),f12(151),f12a(151)
     $/ 2.00000000e+01, 5.98127985e+01, 1.116528541808e-01/
      data eta(152),f12(152),f12a(152)
     $/ 2.02000000e+01, 6.07085314e+01, 1.124395371436e-01/
      data eta(153),f12(153),f12a(153)
     $/ 2.04000000e+01, 6.16087389e+01, 1.097789972451e-01/
      data eta(154),f12(154),f12a(154)
     $/ 2.06000000e+01, 6.25133677e+01, 1.116394738743e-01/
      data eta(155),f12(155),f12a(155)
     $/ 2.08000000e+01, 6.34224367e+01, 1.096931072587e-01/
      data eta(156),f12(156),f12a(156)
     $/ 2.10000000e+01, 6.43359013e+01, 1.089280970924e-01/
      data eta(157),f12(157),f12a(157)
     $/ 2.12000000e+01, 6.52537327e+01, 1.096145043705e-01/
      data eta(158),f12(158),f12a(158)
     $/ 2.14000000e+01, 6.61759281e+01, 1.072138854249e-01/
      data eta(159),f12(159),f12a(159)
     $/ 2.16000000e+01, 6.71024418e+01, 1.092749539287e-01/
      data eta(160),f12(160),f12a(160)
     $/ 2.18000000e+01, 6.80332966e+01, 1.068512988636e-01/
      data eta(161),f12(161),f12a(161)
     $/ 2.20000000e+01, 6.89684391e+01, 1.064748506159e-01/
      data eta(162),f12(162),f12a(162)
     $/ 2.22000000e+01, 6.99078465e+01, 1.069842986745e-01/
      data eta(163),f12(163),f12a(163)
     $/ 2.24000000e+01, 7.08515186e+01, 1.052929546852e-01/
      data eta(164),f12(164),f12a(164)
     $/ 2.26000000e+01, 7.17994175e+01, 1.058638825818e-01/
      data eta(165),f12(165),f12a(165)
     $/ 2.28000000e+01, 7.27515431e+01, 1.052565149906e-01/
      data eta(166),f12(166),f12a(166)
     $/ 2.30000000e+01, 7.37078724e+01, 1.036650574545e-01/
      data eta(167),f12(167),f12a(167)
     $/ 2.32000000e+01, 7.46683674e+01, 1.049382551915e-01/
      data eta(168),f12(168),f12a(168)
     $/ 2.34000000e+01, 7.56330357e+01, 1.025769217803e-01/
      data eta(169),f12(169),f12a(169)
     $/ 2.36000000e+01, 7.66018314e+01, 1.038640576851e-01/
      data eta(170),f12(170),f12a(170)
     $/ 2.38000000e+01, 7.75747700e+01, 1.034018474827e-01/
      data eta(171),f12(171),f12a(171)
     $/ 2.40000000e+01, 7.85518284e+01, 1.004985523838e-01/
      data eta(172),f12(172),f12a(172)
     $/ 2.42000000e+01, 7.95329456e+01, 1.034239429815e-01/
      data eta(173),f12(173),f12a(173)
     $/ 2.44000000e+01, 8.05181599e+01, 1.003706756894e-01/
      data eta(174),f12(174),f12a(174)
     $/ 2.46000000e+01, 8.15074177e+01, 1.016183542594e-01/
      data eta(175),f12(175),f12a(175)
     $/ 2.48000000e+01, 8.25007267e+01, 1.008359072745e-01/
      data eta(176),f12(176),f12a(176)
     $/ 2.50000000e+01, 8.34980640e+01, 9.928301664272e-02/
      data eta(177),f12(177),f12a(177)
     $/ 2.52000000e+01, 8.44993916e+01, 1.005770261543e-01/
      data eta(178),f12(178),f12a(178)
     $/ 2.54000000e+01, 8.55047245e+01, 9.920387874048e-02/
      data eta(179),f12(179),f12a(179)
     $/ 2.56000000e+01, 8.65140324e+01, 9.885745888275e-02/
      data eta(180),f12(180),f12a(180)
     $/ 2.58000000e+01, 8.75272999e+01, 9.930628572935e-02/
      data eta(181),f12(181),f12a(181)
     $/ 2.60000000e+01, 8.85445194e+01, 9.671739819958e-02/
      data eta(182),f12(182),f12a(182)
     $/ 2.62000000e+01, 8.95656452e+01, 9.976912147403e-02/
      data eta(183),f12(183),f12a(183)
     $/ 2.64000000e+01, 9.05907154e+01, 9.586611590140e-02/
      data eta(184),f12(184),f12a(184)
     $/ 2.66000000e+01, 9.16196613e+01, 9.812141492108e-02/
      data eta(185),f12(185),f12a(185)
     $/ 2.68000000e+01, 9.26525059e+01, 9.645322441392e-02/
      data eta(186),f12(186),f12a(186)
     $/ 2.70000000e+01, 9.36892185e+01, 9.626568742494e-02/
      data eta(187),f12(187),f12a(187)
     $/ 2.72000000e+01, 9.47297840e+01, 9.641902588607e-02/
      data eta(188),f12(188),f12a(188)
     $/ 2.74000000e+01, 9.57741947e+01, 9.483820903090e-02/
      data eta(189),f12(189),f12a(189)
     $/ 2.76000000e+01, 9.68224201e+01, 9.643313798854e-02/
      data eta(190),f12(190),f12a(190)
     $/ 2.78000000e+01, 9.78744755e+01, 9.392923901678e-02/
      data eta(191),f12(191),f12a(191)
     $/ 2.80000000e+01, 9.89303150e+01, 9.546490594284e-02/
      data eta(192),f12(192),f12a(192)
     $/ 2.82000000e+01, 9.99899540e+01, 9.413613721391e-02/
      data eta(193),f12(193),f12a(193)
     $/ 2.84000000e+01, 1.01053362e+02, 9.334054520234e-02/
      data eta(194),f12(194),f12a(194)
     $/ 2.86000000e+01, 1.02120516e+02, 9.440168197065e-02/
      data eta(195),f12(195),f12a(195)
     $/ 2.88000000e+01, 1.03191431e+02, 9.320272691961e-02/
      data eta(196),f12(196),f12a(196)
     $/ 2.90000000e+01, 1.04266077e+02, 9.243741035120e-02/
      data eta(197),f12(197),f12a(197)
     $/ 2.92000000e+01, 1.05344431e+02, 9.324763167608e-02/
      data eta(198),f12(198),f12a(198)
     $/ 2.94000000e+01, 1.06426508e+02, 9.302206294354e-02/
      data eta(199),f12(199),f12a(199)
     $/ 2.96000000e+01, 1.07512262e+02, 8.621411654839e-02/
      data eta(200),f12(200),f12a(200)
     $/ 2.98000000e+01, 1.08601709e+02, 1.160714708634e-01/
      data eta(201),f12(201),f12a(201)
     $/ 3.00000000e+01, 1.09694826e+02, 0.000000000000e+00/

c
c
c          if y is between -10 and 30 then use spline table
      if((y.le.30.0d0).and.(y.ge.-10.0d0)) then
        call intrp(y,f1,eta,f12,f12a,n,nlow,nhigh)
c          else if y is greater than 30 use degenerate approximation
      elseif(y.gt.30.0d0) then
        x2=y**(-2)
        f1=a(6)*y*sqrt(y)*th*(5.0d0+(a(1)+(3*a(2)+7.0d0*x2*a(3))*x2)*x2)  
c          else if y is less than -10 use nondegenerate approximation
      else
        f0=dexp(y)   
        f1=a(7)*th*f0*(2.0d0-(4.0d0*a(4)-(6.0d0*a(5)-0.25d0*f0)*f0)*f0) 
      endif
c
      f_1_2=f1  
 999  return
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       f_3_2
c    type:         double precision function
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    call line:    f_3_2(y)      (3/2th fermi integral)
c
c    inputs:       y (double precision)   (argument)
c
c    return:       3/2th fermi integral (double precision)
c
c***********************************************************************
      double precision function f_3_2(y)
      implicit none 
      integer n, nlow, nhigh
      parameter(n=201)
      double precision y, a(7), eta(n), f32(n), f32a(n), x2, f0, f1
c
c                retain values between calls
      save nlow,nhigh
c
      data a/6.16850274d0,1.77568655d0,6.92965606d0,.176776695d0,
     $ 6.41500299d-02,.4d0,1.32934039d0/
      data nlow, nhigh/1,n/
c
c
c                    cubic spline data
      data eta(  1),f32(  1),f32a(  1)
     $/-1.00000000e+01, 6.03514791e-05, 0.000000000000e+00/
      data eta(  2),f32(  2),f32a(  2)
     $/-9.80000000e+00, 7.37133833e-05, 8.958056137096e-05/
      data eta(  3),f32(  3),f32a(  3)
     $/-9.60000000e+00, 9.00335163e-05, 8.541207451613e-05/
      data eta(  4),f32(  4),f32a(  4)
     $/-9.40000000e+00, 1.09966868e-04, 1.107539455645e-04/
      data eta(  5),f32(  5),f32a(  5)
     $/-9.20000000e+00, 1.34313381e-04, 1.335463382257e-04/
      data eta(  6),f32(  6),f32a(  6)
     $/-9.00000000e+00, 1.64050060e-04, 1.635856015326e-04/
      data eta(  7),f32(  7),f32a(  7)
     $/-8.80000000e+00, 2.00370374e-04, 1.996565056440e-04/
      data eta(  8),f32(  8),f32a(  8)
     $/-8.60000000e+00, 2.44731440e-04, 2.439011758913e-04/
      data eta(  9),f32(  9),f32a(  9)
     $/-8.40000000e+00, 2.98913459e-04, 2.978817407908e-04/
      data eta( 10),f32( 10),f32a( 10)
     $/-8.20000000e+00, 3.65090447e-04, 3.638172109455e-04/
      data eta( 11),f32( 11),f32a( 11)
     $/-8.00000000e+00, 4.45917576e-04, 4.443705654272e-04/
      data eta( 12),f32( 12),f32a( 12)
     $/-7.80000000e+00, 5.44637980e-04, 5.426917773457e-04/
      data eta( 13),f32( 13),f32a( 13)
     $/-7.60000000e+00, 6.65211541e-04, 6.628358751901e-04/
      data eta( 14),f32( 14),f32a( 14)
     $/-7.40000000e+00, 8.12475410e-04, 8.095109218942e-04/
      data eta( 15),f32( 15),f32a( 15)
     $/-7.20000000e+00, 9.92335991e-04, 9.886272372330e-04/
      data eta( 16),f32( 16),f32a( 16)
     $/-7.00000000e+00, 1.21200623e-03, 1.207428829174e-03/
      data eta( 17),f32( 17),f32a( 17)
     $/-6.80000000e+00, 1.48029524e-03, 1.474473096071e-03/
      data eta( 18),f32( 18),f32a( 18)
     $/-6.60000000e+00, 1.80795780e-03, 1.800711286541e-03/
      data eta( 19),f32( 19),f32a( 19)
     $/-6.40000000e+00, 2.20812770e-03, 2.198782757766e-03/
      data eta( 20),f32( 20),f32a( 20)
     $/-6.20000000e+00, 2.69683736e-03, 2.685121682395e-03/
      data eta( 21),f32( 21),f32a( 21)
     $/-6.00000000e+00, 3.29366449e-03, 3.278351012652e-03/
      data eta( 22),f32( 22),f32a( 22)
     $/-5.80000000e+00, 4.02250013e-03, 4.002750766996e-03/
      data eta( 23),f32( 23),f32a( 23)
     $/-5.60000000e+00, 4.91251063e-03, 4.886874919365e-03/
      data eta( 24),f32( 24),f32a( 24)
     $/-5.40000000e+00, 5.99928957e-03, 5.965015555545e-03/
      data eta( 25),f32( 25),f32a( 25)
     $/-5.20000000e+00, 7.32625567e-03, 7.281136858453e-03/
      data eta( 26),f32( 26),f32a( 26)
     $/-5.00000000e+00, 8.94638640e-03, 8.885131510642e-03/
      data eta( 27),f32( 27),f32a( 27)
     $/-4.80000000e+00, 1.09242697e-02, 1.084122259898e-02/
      data eta( 28),f32( 28),f32a( 28)
     $/-4.60000000e+00, 1.33386545e-02, 1.322520309344e-02/
      data eta( 29),f32( 29),f32a( 29)
     $/-4.40000000e+00, 1.62855038e-02, 1.612764002725e-02/
      data eta( 30),f32( 30),f32a( 30)
     $/-4.20000000e+00, 1.98816718e-02, 1.966204179756e-02/
      data eta( 31),f32( 31),f32a( 31)
     $/-4.00000000e+00, 2.42694099e-02, 2.395970778252e-02/
      data eta( 32),f32( 32),f32a( 32)
     $/-3.80000000e+00, 2.96217115e-02, 2.918365207237e-02/
      data eta( 33),f32( 33),f32a( 33)
     $/-3.60000000e+00, 3.61488024e-02, 3.552407892800e-02/
      data eta( 34),f32( 34),f32a( 34)
     $/-3.40000000e+00, 4.41058287e-02, 4.321034221562e-02/
      data eta( 35),f32( 35),f32a( 35)
     $/-3.20000000e+00, 5.38020553e-02, 5.251459720952e-02/
      data eta( 36),f32( 36),f32a( 36)
     $/-3.00000000e+00, 6.56117518e-02, 6.375175394630e-02/
      data eta( 37),f32( 37),f32a( 37)
     $/-2.80000000e+00, 7.99869169e-02, 7.729867700529e-02/
      data eta( 38),f32( 38),f32a( 38)
     $/-2.60000000e+00, 9.74722300e-02, 9.357573803256e-02/
      data eta( 39),f32( 39),f32a( 39)
     $/-2.40000000e+00, 1.18722076e-01, 1.130783058645e-01/
      data eta( 40),f32( 40),f32a( 40)
     $/-2.20000000e+00, 1.44520123e-01, 1.363411885096e-01/
      data eta( 41),f32( 41),f32a( 41)
     $/-2.00000000e+00, 1.75800983e-01, 1.639788900970e-01/
      data eta( 42),f32( 42),f32a( 42)
     $/-1.80000000e+00, 2.13674326e-01, 1.966157011023e-01/
      data eta( 43),f32( 43),f32a( 43)
     $/-1.60000000e+00, 2.59450115e-01, 2.349252054937e-01/
      data eta( 44),f32( 44),f32a( 44)
     $/-1.40000000e+00, 3.14665116e-01, 2.795652769229e-01/
      data eta( 45),f32( 45),f32a( 45)
     $/-1.20000000e+00, 3.81109037e-01, 3.311516868147e-01/
      data eta( 46),f32( 46),f32a( 46)
     $/-1.00000000e+00, 4.60848816e-01, 3.902066758184e-01/
      data eta( 47),f32( 47),f32a( 47)
     $/-8.00000000e-01, 5.56249276e-01, 4.571237599118e-01/
      data eta( 48),f32( 48),f32a( 48)
     $/-6.00000000e-01, 6.69988349e-01, 5.320902345345e-01/
      data eta( 49),f32( 49),f32a( 49)
     $/-4.00000000e-01, 8.05064514e-01, 6.150791019500e-01/
      data eta( 50),f32( 50),f32a( 50)
     $/-2.00000000e-01, 9.64795128e-01, 7.057607076653e-01/
      data eta( 51),f32( 51),f32a( 51)
     $/ 0.00000000e+00, 1.15280381e+00, 8.035882673888e-01/
      data eta( 52),f32( 52),f32a( 52)
     $/ 2.00000000e-01, 1.37299815e+00, 9.077349227796e-01/
      data eta( 53),f32( 53),f32a( 53)
     $/ 4.00000000e-01, 1.62953690e+00, 1.017133541493e+00/
      data eta( 54),f32( 54),f32a( 54)
     $/ 6.00000000e-01, 1.92678872e+00, 1.130691411249e+00/
      data eta( 55),f32( 55),f32a( 55)
     $/ 8.00000000e-01, 2.26928741e+00, 1.247131313511e+00/
      data eta( 56),f32( 56),f32a( 56)
     $/ 1.00000000e+00, 2.66168267e+00, 1.365268834706e+00/
      data eta( 57),f32( 57),f32a( 57)
     $/ 1.20000000e+00, 3.10869223e+00, 1.483938347666e+00/
      data eta( 58),f32( 58),f32a( 58)
     $/ 1.40000000e+00, 3.61505681e+00, 1.602230774630e+00/
      data eta( 59),f32( 59),f32a( 59)
     $/ 1.60000000e+00, 4.18550170e+00, 1.719185053814e+00/
      data eta( 60),f32( 60),f32a( 60)
     $/ 1.80000000e+00, 4.82470143e+00, 1.834255010113e+00/
      data eta( 61),f32( 61),f32a( 61)
     $/ 2.00000000e+00, 5.53725398e+00, 1.946717905735e+00/
      data eta( 62),f32( 62),f32a( 62)
     $/ 2.20000000e+00, 6.32765782e+00, 2.056566866947e+00/
      data eta( 63),f32( 63),f32a( 63)
     $/ 2.40000000e+00, 7.20030272e+00, 2.163173626476e+00/
      data eta( 64),f32( 64),f32a( 64)
     $/ 2.60000000e+00, 8.15945459e+00, 2.266784127147e+00/
      data eta( 65),f32( 65),f32a( 65)
     $/ 2.80000000e+00, 9.20925546e+00, 2.367039864934e+00/
      data eta( 66),f32( 66),f32a( 66)
     $/ 3.00000000e+00, 1.03537161e+01, 2.464021913115e+00/
      data eta( 67),f32( 67),f32a( 67)
     $/ 3.20000000e+00, 1.15967200e+01, 2.558361482607e+00/
      data eta( 68),f32( 68),f32a( 68)
     $/ 3.40000000e+00, 1.29420350e+01, 2.649197156458e+00/
      data eta( 69),f32( 69),f32a( 69)
     $/ 3.60000000e+00, 1.43933013e+01, 2.737544891559e+00/
      data eta( 70),f32( 70),f32a( 70)
     $/ 3.80000000e+00, 1.59540503e+01, 2.823028277308e+00/
      data eta( 71),f32( 71),f32a( 71)
     $/ 4.00000000e+00, 1.76277032e+01, 2.905926999207e+00/
      data eta( 72),f32( 72),f32a( 72)
     $/ 4.20000000e+00, 1.94175764e+01, 2.986308725864e+00/
      data eta( 73),f32( 73),f32a( 73)
     $/ 4.40000000e+00, 2.13268934e+01, 3.065408097336e+00/
      data eta( 74),f32( 74),f32a( 74)
     $/ 4.60000000e+00, 2.33587976e+01, 3.140138884795e+00/
      data eta( 75),f32( 75),f32a( 75)
     $/ 4.80000000e+00, 2.55163160e+01, 3.216166363482e+00/
      data eta( 76),f32( 76),f32a( 76)
     $/ 5.00000000e+00, 2.78024469e+01, 3.287070661277e+00/
      data eta( 77),f32( 77),f32a( 77)
     $/ 5.20000000e+00, 3.02200609e+01, 3.358015991410e+00/
      data eta( 78),f32( 78),f32a( 78)
     $/ 5.40000000e+00, 3.27719889e+01, 3.427965373083e+00/
      data eta( 79),f32( 79),f32a( 79)
     $/ 5.60000000e+00, 3.54610072e+01, 3.493667516261e+00/
      data eta( 80),f32( 80),f32a( 80)
     $/ 5.80000000e+00, 3.82897883e+01, 3.561784561870e+00/
      data eta( 81),f32( 81),f32a( 81)
     $/ 6.00000000e+00, 4.12610064e+01, 3.624744236258e+00/
      data eta( 82),f32( 82),f32a( 82)
     $/ 6.20000000e+00, 4.43772212e+01, 3.688743493097e+00/
      data eta( 83),f32( 83),f32a( 83)
     $/ 6.40000000e+00, 4.76409770e+01, 3.751431791354e+00/
      data eta( 84),f32( 84),f32a( 84)
     $/ 6.60000000e+00, 5.10547763e+01, 3.812054341488e+00/
      data eta( 85),f32( 85),f32a( 85)
     $/ 6.80000000e+00, 5.46210566e+01, 3.872500842692e+00/
      data eta( 86),f32( 86),f32a( 86)
     $/ 7.00000000e+00, 5.83422213e+01, 3.930602287744e+00/
      data eta( 87),f32( 87),f32a( 87)
     $/ 7.20000000e+00, 6.22206164e+01, 3.989650006333e+00/
      data eta( 88),f32( 88),f32a( 88)
     $/ 7.40000000e+00, 6.62585851e+01, 4.046837686926e+00/
      data eta( 89),f32( 89),f32a( 89)
     $/ 7.60000000e+00, 7.04584142e+01, 4.102059245965e+00/
      data eta( 90),f32( 90),f32a( 90)
     $/ 7.80000000e+00, 7.48223363e+01, 4.158875329212e+00/
      data eta( 91),f32( 91),f32a( 91)
     $/ 8.00000000e+00, 7.93525945e+01, 4.212854437187e+00/
      data eta( 92),f32( 92),f32a( 92)
     $/ 8.20000000e+00, 8.40513631e+01, 4.266266922047e+00/
      data eta( 93),f32( 93),f32a( 93)
     $/ 8.40000000e+00, 8.89207860e+01, 4.320222874619e+00/
      data eta( 94),f32( 94),f32a( 94)
     $/ 8.60000000e+00, 9.39630071e+01, 4.372571579482e+00/
      data eta( 95),f32( 95),f32a( 95)
     $/ 8.80000000e+00, 9.91801320e+01, 4.425060807447e+00/
      data eta( 96),f32( 96),f32a( 96)
     $/ 9.00000000e+00, 1.04574244e+02, 4.475250190736e+00/
      data eta( 97),f32( 97),f32a( 97)
     $/ 9.20000000e+00, 1.10147364e+02, 4.525138429613e+00/
      data eta( 98),f32( 98),f32a( 98)
     $/ 9.40000000e+00, 1.15901507e+02, 4.577646090805e+00/
      data eta( 99),f32( 99),f32a( 99)
     $/ 9.60000000e+00, 1.21838717e+02, 4.624327207175e+00/
      data eta(100),f32(100),f32a(100)
     $/ 9.80000000e+00, 1.27960940e+02, 4.676995080485e+00/
      data eta(101),f32(101),f32a(101)
     $/ 1.00000000e+01, 1.34270176e+02, 4.719642470896e+00/
      data eta(102),f32(102),f32a(102)
     $/ 1.02000000e+01, 1.40768269e+02, 4.772985035935e+00/
      data eta(103),f32(103),f32a(103)
     $/ 1.04000000e+01, 1.47457218e+02, 4.816817385355e+00/
      data eta(104),f32(104),f32a(104)
     $/ 1.06000000e+01, 1.54338871e+02, 4.865345422653e+00/
      data eta(105),f32(105),f32a(105)
     $/ 1.08000000e+01, 1.61415135e+02, 4.913450924021e+00/
      data eta(106),f32(106),f32a(106)
     $/ 1.10000000e+01, 1.68687886e+02, 4.953900881277e+00/
      data eta(107),f32(107),f32a(107)
     $/ 1.12000000e+01, 1.76158863e+02, 5.004845550873e+00/
      data eta(108),f32(108),f32a(108)
     $/ 1.14000000e+01, 1.83829975e+02, 5.046966915220e+00/
      data eta(109),f32(109),f32a(109)
     $/ 1.16000000e+01, 1.91702992e+02, 5.093036788259e+00/
      data eta(110),f32(110),f32a(110)
     $/ 1.18000000e+01, 1.99779728e+02, 5.138735931733e+00/
      data eta(111),f32(111),f32a(111)
     $/ 1.20000000e+01, 2.08061970e+02, 5.177919484819e+00/
      data eta(112),f32(112),f32a(112)
     $/ 1.22000000e+01, 2.16551396e+02, 5.227186128993e+00/
      data eta(113),f32(113),f32a(113)
     $/ 1.24000000e+01, 2.25249821e+02, 5.263185999195e+00/
      data eta(114),f32(114),f32a(114)
     $/ 1.26000000e+01, 2.34158879e+02, 5.315019874243e+00/
      data eta(115),f32(115),f32a(115)
     $/ 1.28000000e+01, 2.43280430e+02, 5.350684503817e+00/
      data eta(116),f32(116),f32a(116)
     $/ 1.30000000e+01, 2.52616062e+02, 5.394392110504e+00/
      data eta(117),f32(117),f32a(117)
     $/ 1.32000000e+01, 2.62167458e+02, 5.436347054167e+00/
      data eta(118),f32(118),f32a(118)
     $/ 1.34000000e+01, 2.71936318e+02, 5.479819672817e+00/
      data eta(119),f32(119),f32a(119)
     $/ 1.36000000e+01, 2.81924325e+02, 5.516424254575e+00/
      data eta(120),f32(120),f32a(120)
     $/ 1.38000000e+01, 2.92133065e+02, 5.564433308870e+00/
      data eta(121),f32(121),f32a(121)
     $/ 1.40000000e+01, 3.02564278e+02, 5.596792509966e+00/
      data eta(122),f32(122),f32a(122)
     $/ 1.42000000e+01, 3.13219429e+02, 5.639096651265e+00/
      data eta(123),f32(123),f32a(123)
     $/ 1.44000000e+01, 3.24100167e+02, 5.684870884964e+00/
      data eta(124),f32(124),f32a(124)
     $/ 1.46000000e+01, 3.35208199e+02, 5.715519808889e+00/
      data eta(125),f32(125),f32a(125)
     $/ 1.48000000e+01, 3.46544991e+02, 5.767049879467e+00/
      data eta(126),f32(126),f32a(126)
     $/ 1.50000000e+01, 3.58112282e+02, 5.791130673258e+00/
      data eta(127),f32(127),f32a(127)
     $/ 1.52000000e+01, 3.69911385e+02, 5.840227427499e+00/
      data eta(128),f32(128),f32a(128)
     $/ 1.54000000e+01, 3.81944008e+02, 5.875959616731e+00/
      data eta(129),f32(129),f32a(129)
     $/ 1.56000000e+01, 3.94211678e+02, 5.912984105596e+00/
      data eta(130),f32(130),f32a(130)
     $/ 1.58000000e+01, 4.06715920e+02, 5.957903960865e+00/
      data eta(131),f32(131),f32a(131)
     $/ 1.60000000e+01, 4.19458352e+02, 5.983900050964e+00/
      data eta(132),f32(132),f32a(132)
     $/ 1.62000000e+01, 4.32440285e+02, 6.031645835280e+00/
      data eta(133),f32(133),f32a(133)
     $/ 1.64000000e+01, 4.45663338e+02, 6.057516607911e+00/
      data eta(134),f32(134),f32a(134)
     $/ 1.66000000e+01, 4.59128884e+02, 6.112237733043e+00/
      data eta(135),f32(135),f32a(135)
     $/ 1.68000000e+01, 4.72838723e+02, 6.137482459949e+00/
      data eta(136),f32(136),f32a(136)
     $/ 1.70000000e+01, 4.86794106e+02, 6.169432427160e+00/
      data eta(137),f32(137),f32a(137)
     $/ 1.72000000e+01, 5.00996376e+02, 6.217837831413e+00/
      data eta(138),f32(138),f32a(138)
     $/ 1.74000000e+01, 5.15447221e+02, 6.245466247192e+00/
      data eta(139),f32(139),f32a(139)
     $/ 1.76000000e+01, 5.30147965e+02, 6.285147179775e+00/
      data eta(140),f32(140),f32a(140)
     $/ 1.78000000e+01, 5.45100114e+02, 6.324695033745e+00/
      data eta(141),f32(141),f32a(141)
     $/ 1.80000000e+01, 5.60305131e+02, 6.346272685250e+00/
      data eta(142),f32(142),f32a(142)
     $/ 1.82000000e+01, 5.75764237e+02, 6.403564225258e+00/
      data eta(143),f32(143),f32a(143)
     $/ 1.84000000e+01, 5.91479142e+02, 6.409320413722e+00/
      data eta(144),f32(144),f32a(144)
     $/ 1.86000000e+01, 6.07450822e+02, 6.475404119805e+00/
      data eta(145),f32(145),f32a(145)
     $/ 1.88000000e+01, 6.23681230e+02, 6.498263107104e+00/
      data eta(146),f32(146),f32a(146)
     $/ 1.90000000e+01, 6.40171525e+02, 6.514593451770e+00/
      data eta(147),f32(147),f32a(147)
     $/ 1.92000000e+01, 6.56922746e+02, 6.582263085821e+00/
      data eta(148),f32(148),f32a(148)
     $/ 1.94000000e+01, 6.73936845e+02, 6.588054204946e+00/
      data eta(149),f32(149),f32a(149)
     $/ 1.96000000e+01, 6.91214799e+02, 6.643770094354e+00/
      data eta(150),f32(150),f32a(150)
     $/ 1.98000000e+01, 7.08758256e+02, 6.662315417678e+00/
      data eta(151),f32(151),f32a(151)
     $/ 2.00000000e+01, 7.26568315e+02, 6.697268234941e+00/
      data eta(152),f32(152),f32a(152)
     $/ 2.02000000e+01, 7.44646317e+02, 6.740061642552e+00/
      data eta(153),f32(153),f32a(153)
     $/ 2.04000000e+01, 7.62993791e+02, 6.763285194860e+00/
      data eta(154),f32(154),f32a(154)
     $/ 2.06000000e+01, 7.81611955e+02, 6.810297577950e+00/
      data eta(155),f32(155),f32a(155)
     $/ 2.08000000e+01, 8.00502335e+02, 6.827924493401e+00/
      data eta(156),f32(156),f32a(156)
     $/ 2.10000000e+01, 8.19665971e+02, 6.866404448437e+00/
      data eta(157),f32(157),f32a(157)
     $/ 2.12000000e+01, 8.39104264e+02, 6.905007712841e+00/
      data eta(158),f32(158),f32a(158)
     $/ 2.14000000e+01, 8.58818620e+02, 6.923014700220e+00/
      data eta(159),f32(159),f32a(159)
     $/ 2.16000000e+01, 8.78810136e+02, 6.976933486221e+00/
      data eta(160),f32(160),f32a(160)
     $/ 2.18000000e+01, 8.99080522e+02, 6.999751354928e+00/
      data eta(161),f32(161),f32a(161)
     $/ 2.20000000e+01, 9.19630815e+02, 7.010111094081e+00/
      data eta(162),f32(162),f32a(162)
     $/ 2.22000000e+01, 9.40461930e+02, 7.083104268755e+00/
      data eta(163),f32(163),f32a(163)
     $/ 2.24000000e+01, 9.61575822e+02, 7.074021830882e+00/
      data eta(164),f32(164),f32a(164)
     $/ 2.26000000e+01, 9.82973161e+02, 7.137858407662e+00/
      data eta(165),f32(165),f32a(165)
     $/ 2.28000000e+01, 1.00465578e+03, 7.166544538534e+00/
      data eta(166),f32(166),f32a(166)
     $/ 2.30000000e+01, 1.02662479e+03, 7.154613438199e+00/
      data eta(167),f32(167),f32a(167)
     $/ 2.32000000e+01, 1.04888077e+03, 7.260501708652e+00/
      data eta(168),f32(168),f32a(168)
     $/ 2.34000000e+01, 1.07142618e+03, 7.217879727200e+00/
      data eta(169),f32(169),f32a(169)
     $/ 2.36000000e+01, 1.09426114e+03, 7.300479382516e+00/
      data eta(170),f32(170),f32a(170)
     $/ 2.38000000e+01, 1.11738761e+03, 7.306702742765e+00/
      data eta(171),f32(171),f32a(171)
     $/ 2.40000000e+01, 1.14080655e+03, 7.343209646437e+00/
      data eta(172),f32(172),f32a(172)
     $/ 2.42000000e+01, 1.16451920e+03, 7.376958671475e+00/
      data eta(173),f32(173),f32a(173)
     $/ 2.44000000e+01, 1.18852677e+03, 7.386955667671e+00/
      data eta(174),f32(174),f32a(174)
     $/ 2.46000000e+01, 1.21283023e+03, 7.458718657792e+00/
      data eta(175),f32(175),f32a(175)
     $/ 2.48000000e+01, 1.23743155e+03, 7.457169701212e+00/
      data eta(176),f32(176),f32a(176)
     $/ 2.50000000e+01, 1.26233133e+03, 7.481602537356e+00/
      data eta(177),f32(177),f32a(177)
     $/ 2.52000000e+01, 1.28753067e+03, 7.550420149389e+00/
      data eta(178),f32(178),f32a(178)
     $/ 2.54000000e+01, 1.31303140e+03, 7.525216865068e+00/
      data eta(179),f32(179),f32a(179)
     $/ 2.56000000e+01, 1.33883389e+03, 7.612712390278e+00/
      data eta(180),f32(180),f32a(180)
     $/ 2.58000000e+01, 1.36494022e+03, 7.599933573862e+00/
      data eta(181),f32(181),f32a(181)
     $/ 2.60000000e+01, 1.39135086e+03, 7.634053314309e+00/
      data eta(182),f32(182),f32a(182)
     $/ 2.62000000e+01, 1.41806705e+03, 7.696353168910e+00/
      data eta(183),f32(183),f32a(183)
     $/ 2.64000000e+01, 1.44509061e+03, 7.686034010016e+00/
      data eta(184),f32(184),f32a(184)
     $/ 2.66000000e+01, 1.47242203e+03, 7.738510790971e+00/
      data eta(185),f32(185),f32a(185)
     $/ 2.68000000e+01, 1.50006278e+03, 7.759422826178e+00/
      data eta(186),f32(186),f32a(186)
     $/ 2.70000000e+01, 1.52801395e+03, 7.786797904321e+00/
      data eta(187),f32(187),f32a(187)
     $/ 2.72000000e+01, 1.55627664e+03, 7.821385556534e+00/
      data eta(188),f32(188),f32a(188)
     $/ 2.74000000e+01, 1.58485208e+03, 7.840159869515e+00/
      data eta(189),f32(189),f32a(189)
     $/ 2.76000000e+01, 1.61374137e+03, 7.895474965363e+00/
      data eta(190),f32(190),f32a(190)
     $/ 2.78000000e+01, 1.64294608e+03, 7.890940269096e+00/
      data eta(191),f32(191),f32a(191)
     $/ 2.80000000e+01, 1.67246671e+03, 7.928763958255e+00/
      data eta(192),f32(192),f32a(192)
     $/ 2.82000000e+01, 1.70230460e+03, 7.983003897891e+00/
      data eta(193),f32(193),f32a(193)
     $/ 2.84000000e+01, 1.73246121e+03, 7.947220450177e+00/
      data eta(194),f32(194),f32a(194)
     $/ 2.86000000e+01, 1.76293668e+03, 8.057114301331e+00/
      data eta(195),f32(195),f32a(195)
     $/ 2.88000000e+01, 1.79373355e+03, 8.034322344553e+00/
      data eta(196),f32(196),f32a(196)
     $/ 2.90000000e+01, 1.82485221e+03, 8.074096320488e+00/
      data eta(197),f32(197),f32a(197)
     $/ 2.92000000e+01, 1.85629361e+03, 8.080292373472e+00/
      data eta(198),f32(198),f32a(198)
     $/ 2.94000000e+01, 1.88805936e+03, 8.257234185640e+00/
      data eta(199),f32(199),f32a(199)
     $/ 2.96000000e+01, 1.92014981e+03, 7.595770883867e+00/
      data eta(200),f32(200),f32a(200)
     $/ 2.98000000e+01, 1.95256705e+03, 1.037818227902e+01/
      data eta(201),f32(201),f32a(201)
     $/ 3.00000000e+01, 1.98531168e+03, 0.000000000000e+00/
c
c          if y is between -10 and 30 then use spline table
      if((y.le.30.0d0).and.(y.ge.-10.0d0)) then
        call intrp(y,f1,eta,f32,f32a,n,nlow,nhigh)
c          else if y is greater than 30 use degenerate approximation
      elseif(y.gt.30.0d0) then
        x2=y**(-2)
        f1=a(6)*sqrt(y)*(1.0d0/x2+a(1)-(a(2)+x2*a(3))*x2)  
      else
c          else if y is less than -10 use nondegenerate approximation
        f0=dexp(y)   
        f1=a(7)*f0*(1.0d0-(a(4)-(a(5)-0.03125d0*f0)*f0)*f0) 
      endif
c
      f_3_2=f1  
 999  return
      end   
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       finv12
c    type:         double precision function
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    call line:    finv12(y)      (inverse of the 1/2th fermi integral)
c
c    inputs:       y (double precision)   (argument)
c
c    return:       inverse of fermi integral (double precision)
c
c***********************************************************************
      double precision function finv12(y)
      implicit none
      integer n, nlow, nhigh
      parameter(n=201)
      double precision y, ai(8),f12(n),fia(n),eta(n), x2, x4, f1
c
c                retain values between calls
      save nlow,nhigh
c
      data ai/-.822467032d0,-1.21761363d0,-9.16138616d0,
     $ 0.398942281d0,.0732748216d0,-1.310707d0,1.12837917d0,  
     $ 8.2810645d-3/ 
      data nlow,nhigh/1,n/
c
c
c                    cubic spline data
      data eta(  1),f12(  1),fia(  1)
     $/-1.00000000e+01, 4.02339983e-05, 0.000000000000e+00/
      data eta(  2),f12(  2),fia(  2)
     $/-9.80000000e+00, 4.91417786e-05,-5.518729757333e+08/
      data eta(  3),f12(  3),fia(  3)
     $/-9.60000000e+00, 6.00216308e-05,-2.369113017381e+08/
      data eta(  4),f12(  4),fia(  4)
     $/-9.40000000e+00, 7.33101788e-05,-1.908760407552e+08/
      data eta(  5),f12(  5),fia(  5)
     $/-9.20000000e+00, 8.95406574e-05,-1.202176248026e+08/
      data eta(  6),f12(  6),fia(  6)
     $/-9.00000000e+00, 1.09364317e-04,-8.245454770317e+07/
      data eta(  7),f12(  7),fia(  7)
     $/-8.80000000e+00, 1.33576705e-04,-5.481585297269e+07/
      data eta(  8),f12(  8),fia(  8)
     $/-8.60000000e+00, 1.63148990e-04,-3.685695721234e+07/
      data eta(  9),f12(  9),fia(  9)
     $/-8.40000000e+00, 1.99267717e-04,-2.467974532472e+07/
      data eta( 10),f12( 10),fia( 10)
     $/-8.20000000e+00, 2.43381819e-04,-1.655022961155e+07/
      data eta( 11),f12( 11),fia( 11)
     $/-8.00000000e+00, 2.97260769e-04,-1.109333129137e+07/
      data eta( 12),f12( 12),fia( 12)
     $/-7.80000000e+00, 3.63065723e-04,-7.436285795690e+06/
      data eta( 13),f12( 13),fia( 13)
     $/-7.60000000e+00, 4.43435136e-04,-4.985364210851e+06/
      data eta( 14),f12( 14),fia( 14)
     $/-7.40000000e+00, 5.41591857e-04,-3.341805327658e+06/
      data eta( 15),f12( 15),fia( 15)
     $/-7.20000000e+00, 6.61470011e-04,-2.240407265052e+06/
      data eta( 16),f12( 16),fia( 16)
     $/-7.00000000e+00, 8.07873978e-04,-1.501972226436e+06/
      data eta( 17),f12( 17),fia( 17)
     $/-6.80000000e+00, 9.86669358e-04,-1.006911415258e+06/
      data eta( 18),f12( 18),fia( 18)
     $/-6.60000000e+00, 1.20501561e-03,-6.751030308338e+05/
      data eta( 19),f12( 19),fia( 19)
     $/-6.40000000e+00, 1.47165314e-03,-4.526103593217e+05/
      data eta( 20),f12( 20),fia( 20)
     $/-6.20000000e+00, 1.79724747e-03,-3.034927371978e+05/
      data eta( 21),f12( 21),fia( 21)
     $/-6.00000000e+00, 2.19481561e-03,-2.034893180359e+05/
      data eta( 22),f12( 22),fia( 22)
     $/-5.80000000e+00, 2.68023420e-03,-1.364626916017e+05/
      data eta( 23),f12( 23),fia( 23)
     $/-5.60000000e+00, 3.27287085e-03,-9.151916808167e+04/
      data eta( 24),f12( 24),fia( 24)
     $/-5.40000000e+00, 3.99634114e-03,-6.137988251686e+04/
      data eta( 25),f12( 25),fia( 25)
     $/-5.20000000e+00, 4.87942214e-03,-4.117642510051e+04/
      data eta( 26),f12( 26),fia( 26)
     $/-5.00000000e+00, 5.95717982e-03,-2.762388797406e+04/
      data eta( 27),f12( 27),fia( 27)
     $/-4.80000000e+00, 7.27229903e-03,-1.853723945326e+04/
      data eta( 28),f12( 28),fia( 28)
     $/-4.60000000e+00, 8.87672161e-03,-1.244247237626e+04/
      data eta( 29),f12( 29),fia( 29)
     $/-4.40000000e+00, 1.08335980e-02,-8.353189544372e+03/
      data eta( 30),f12( 30),fia( 30)
     $/-4.20000000e+00, 1.32195998e-02,-5.610484545822e+03/
      data eta( 31),f12( 31),fia( 31)
     $/-4.00000000e+00, 1.61277414e-02,-3.769624980571e+03/
      data eta( 32),f12( 32),fia( 32)
     $/-3.80000000e+00, 1.96706583e-02,-2.534200684298e+03/
      data eta( 33),f12( 33),fia( 33)
     $/-3.60000000e+00, 2.39845128e-02,-1.704676561770e+03/
      data eta( 34),f12( 34),fia( 34)
     $/-3.40000000e+00, 2.92335232e-02,-1.147572737681e+03/
      data eta( 35),f12( 35),fia( 35)
     $/-3.20000000e+00, 3.56152142e-02,-7.732385076806e+02/
      data eta( 36),f12( 36),fia( 36)
     $/-3.00000000e+00, 4.33663810e-02,-5.215910937804e+02/
      data eta( 37),f12( 37),fia( 37)
     $/-2.80000000e+00, 5.27697736e-02,-3.523253054057e+02/
      data eta( 38),f12( 38),fia( 38)
     $/-2.60000000e+00, 6.41614366e-02,-2.383626995336e+02/
      data eta( 39),f12( 39),fia( 39)
     $/-2.40000000e+00, 7.79383797e-02,-1.615785832573e+02/
      data eta( 40),f12( 40),fia( 40)
     $/-2.20000000e+00, 9.45664588e-02,-1.097815116306e+02/
      data eta( 41),f12( 41),fia( 41)
     $/-2.00000000e+00, 1.14587830e-01,-7.479631187407e+01/
      data eta( 42),f12( 42),fia( 42)
     $/-1.80000000e+00, 1.38627354e-01,-5.112413248282e+01/
      data eta( 43),f12( 43),fia( 43)
     $/-1.60000000e+00, 1.67396817e-01,-3.507905015366e+01/
      data eta( 44),f12( 44),fia( 44)
     $/-1.40000000e+00, 2.01696221e-01,-2.417700470772e+01/
      data eta( 45),f12( 45),fia( 45)
     $/-1.20000000e+00, 2.42410529e-01,-1.674984464907e+01/
      data eta( 46),f12( 46),fia( 46)
     $/-1.00000000e+00, 2.90500917e-01,-1.167337092503e+01/
      data eta( 47),f12( 47),fia( 47)
     $/-8.00000000e-01, 3.46989460e-01,-8.190719927811e+00/
      data eta( 48),f12( 48),fia( 48)
     $/-6.00000000e-01, 4.12937023e-01,-5.790648969028e+00/
      data eta( 49),f12( 49),fia( 49)
     $/-4.00000000e-01, 4.89414580e-01,-4.128944296030e+00/
      data eta( 50),f12( 50),fia( 50)
     $/-2.00000000e-01, 5.77470496e-01,-2.971061447648e+00/
      data eta( 51),f12( 51),fia( 51)
     $/ 0.00000000e+00, 6.78093925e-01,-2.159721811273e+00/
      data eta( 52),f12( 52),fia( 52)
     $/ 2.00000000e-01, 7.92181447e-01,-1.586687451373e+00/
      data eta( 53),f12( 53),fia( 53)
     $/ 4.00000000e-01, 9.20506015e-01,-1.178969413999e+00/
      data eta( 54),f12( 54),fia( 54)
     $/ 6.00000000e-01, 1.06369475e+00,-8.863670774363e-01/
      data eta( 55),f12( 55),fia( 55)
     $/ 8.00000000e-01, 1.22221592e+00,-6.744468494356e-01/
      data eta( 56),f12( 56),fia( 56)
     $/ 1.00000000e+00, 1.39637545e+00,-5.194863803958e-01/
      data eta( 57),f12( 57),fia( 57)
     $/ 1.20000000e+00, 1.58632329e+00,-4.051201736579e-01/
      data eta( 58),f12( 58),fia( 58)
     $/ 1.40000000e+00, 1.79206851e+00,-3.197436337379e-01/
      data eta( 59),f12( 59),fia( 59)
     $/ 1.60000000e+00, 2.01349622e+00,-2.554204570394e-01/
      data eta( 60),f12( 60),fia( 60)
     $/ 1.80000000e+00, 2.25039083e+00,-2.064307210907e-01/
      data eta( 61),f12( 61),fia( 61)
     $/ 2.00000000e+00, 2.50245792e+00,-1.687058961462e-01/
      data eta( 62),f12( 62),fia( 62)
     $/ 2.20000000e+00, 2.76934439e+00,-1.394151646274e-01/
      data eta( 63),f12( 63),fia( 63)
     $/ 2.40000000e+00, 3.05065972e+00,-1.163727819558e-01/
      data eta( 64),f12( 64),fia( 64)
     $/ 2.60000000e+00, 3.34598833e+00,-9.810926234532e-02/
      data eta( 65),f12( 65),fia( 65)
     $/ 2.80000000e+00, 3.65490490e+00,-8.349866180907e-02/
      data eta( 66),f12( 66),fia( 66)
     $/ 3.00000000e+00, 3.97698528e+00,-7.167023124865e-02/
      data eta( 67),f12( 67),fia( 67)
     $/ 3.20000000e+00, 4.31181109e+00,-6.203344907539e-02/
      data eta( 68),f12( 68),fia( 68)
     $/ 3.40000000e+00, 4.65897715e+00,-5.410768842926e-02/
      data eta( 69),f12( 69),fia( 69)
     $/ 3.60000000e+00, 5.01809514e+00,-4.753934908271e-02/
      data eta( 70),f12( 70),fia( 70)
     $/ 3.80000000e+00, 5.38879550e+00,-4.203693557621e-02/
      data eta( 71),f12( 71),fia( 71)
     $/ 4.00000000e+00, 5.77072680e+00,-3.741510484228e-02/
      data eta( 72),f12( 72),fia( 72)
     $/ 4.20000000e+00, 6.16355908e+00,-3.349161263478e-02/
      data eta( 73),f12( 73),fia( 73)
     $/ 4.40000000e+00, 6.56698239e+00,-3.014726514264e-02/
      data eta( 74),f12( 74),fia( 74)
     $/ 4.60000000e+00, 6.98070586e+00,-2.725049073768e-02/
      data eta( 75),f12( 75),fia( 75)
     $/ 4.80000000e+00, 7.40445435e+00,-2.478810817325e-02/
      data eta( 76),f12( 76),fia( 76)
     $/ 5.00000000e+00, 7.83797658e+00,-2.259805758188e-02/
      data eta( 77),f12( 77),fia( 77)
     $/ 5.20000000e+00, 8.28102899e+00,-2.071312823782e-02/
      data eta( 78),f12( 78),fia( 78)
     $/ 5.40000000e+00, 8.73338914e+00,-1.906423734281e-02/
      data eta( 79),f12( 79),fia( 79)
     $/ 5.60000000e+00, 9.19485021e+00,-1.756402815298e-02/
      data eta( 80),f12( 80),fia( 80)
     $/ 5.80000000e+00, 9.66520906e+00,-1.628975776444e-02/
      data eta( 81),f12( 81),fia( 81)
     $/ 6.00000000e+00, 1.01442864e+01,-1.510245658209e-02/
      data eta( 82),f12( 82),fia( 82)
     $/ 6.20000000e+00, 1.06319029e+01,-1.407142912242e-02/
      data eta( 83),f12( 83),fia( 83)
     $/ 6.40000000e+00, 1.11278965e+01,-1.314229513765e-02/
      data eta( 84),f12( 84),fia( 84)
     $/ 6.60000000e+00, 1.16321146e+01,-1.228449605398e-02/
      data eta( 85),f12( 85),fia( 85)
     $/ 6.80000000e+00, 1.21444066e+01,-1.153091235501e-02/
      data eta( 86),f12( 86),fia( 86)
     $/ 7.00000000e+00, 1.26646369e+01,-1.083803827012e-02/
      data eta( 87),f12( 87),fia( 87)
     $/ 7.20000000e+00, 1.31926749e+01,-1.019988677668e-02/
      data eta( 88),f12( 88),fia( 88)
     $/ 7.40000000e+00, 1.37283938e+01,-9.631384626898e-03/
      data eta( 89),f12( 89),fia( 89)
     $/ 7.60000000e+00, 1.42716773e+01,-9.093442710361e-03/
      data eta( 90),f12( 90),fia( 90)
     $/ 7.80000000e+00, 1.48224099e+01,-8.618263823991e-03/
      data eta( 91),f12( 91),fia( 91)
     $/ 8.00000000e+00, 1.53804867e+01,-8.160345480975e-03/
      data eta( 92),f12( 92),fia( 92)
     $/ 8.20000000e+00, 1.59458020e+01,-7.762531673437e-03/
      data eta( 93),f12( 93),fia( 93)
     $/ 8.40000000e+00, 1.65182614e+01,-7.360350480125e-03/
      data eta( 94),f12( 94),fia( 94)
     $/ 8.60000000e+00, 1.70977635e+01,-7.030118584822e-03/
      data eta( 95),f12( 95),fia( 95)
     $/ 8.80000000e+00, 1.76842275e+01,-6.688634687327e-03/
      data eta( 96),f12( 96),fia( 96)
     $/ 9.00000000e+00, 1.82775617e+01,-6.382668439251e-03/
      data eta( 97),f12( 97),fia( 97)
     $/ 9.20000000e+00, 1.88776803e+01,-6.100108206762e-03/
      data eta( 98),f12( 98),fia( 98)
     $/ 9.40000000e+00, 1.94845052e+01,-5.838946705110e-03/
      data eta( 99),f12( 99),fia( 99)
     $/ 9.60000000e+00, 2.00979619e+01,-5.584582666939e-03/
      data eta(100),f12(100),fia(100)
     $/ 9.80000000e+00, 2.07179742e+01,-5.361259022162e-03/
      data eta(101),f12(101),fia(101)
     $/ 1.00000000e+01, 2.13444734e+01,-5.126492837303e-03/
      data eta(102),f12(102),fia(102)
     $/ 1.02000000e+01, 2.19773812e+01,-4.934700434242e-03/
      data eta(103),f12(103),fia(103)
     $/ 1.04000000e+01, 2.26166368e+01,-4.735625230311e-03/
      data eta(104),f12(104),fia(104)
     $/ 1.06000000e+01, 2.32621732e+01,-4.556040567237e-03/
      data eta(105),f12(105),fia(105)
     $/ 1.08000000e+01, 2.39139295e+01,-4.391444865609e-03/
      data eta(106),f12(106),fia(106)
     $/ 1.10000000e+01, 2.45718484e+01,-4.216034760485e-03/
      data eta(107),f12(107),fia(107)
     $/ 1.12000000e+01, 2.52358613e+01,-4.071259985244e-03/
      data eta(108),f12(108),fia(108)
     $/ 1.14000000e+01, 2.59059148e+01,-3.924867157171e-03/
      data eta(109),f12(109),fia(109)
     $/ 1.16000000e+01, 2.65819535e+01,-3.788293366519e-03/
      data eta(110),f12(110),fia(110)
     $/ 1.18000000e+01, 2.72639241e+01,-3.661781142180e-03/
      data eta(111),f12(111),fia(111)
     $/ 1.20000000e+01, 2.79517770e+01,-3.535793689118e-03/
      data eta(112),f12(112),fia(112)
     $/ 1.22000000e+01, 2.86454568e+01,-3.410614036319e-03/
      data eta(113),f12(113),fia(113)
     $/ 1.24000000e+01, 2.93449082e+01,-3.315348054643e-03/
      data eta(114),f12(114),fia(114)
     $/ 1.26000000e+01, 3.00500951e+01,-3.196961694843e-03/
      data eta(115),f12(115),fia(115)
     $/ 1.28000000e+01, 3.07609620e+01,-3.101468403870e-03/
      data eta(116),f12(116),fia(116)
     $/ 1.30000000e+01, 3.14774652e+01,-3.004323966790e-03/
      data eta(117),f12(117),fia(117)
     $/ 1.32000000e+01, 3.21995592e+01,-2.909424474011e-03/
      data eta(118),f12(118),fia(118)
     $/ 1.34000000e+01, 3.29271975e+01,-2.827365948037e-03/
      data eta(119),f12(119),fia(119)
     $/ 1.36000000e+01, 3.36603403e+01,-2.734508928049e-03/
      data eta(120),f12(120),fia(120)
     $/ 1.38000000e+01, 3.43989420e+01,-2.674531887656e-03/
      data eta(121),f12(121),fia(121)
     $/ 1.40000000e+01, 3.51429758e+01,-2.568768499103e-03/
      data eta(122),f12(122),fia(122)
     $/ 1.42000000e+01, 3.58923769e+01,-2.520360850761e-03/
      data eta(123),f12(123),fia(123)
     $/ 1.44000000e+01, 3.66471262e+01,-2.437521643781e-03/
      data eta(124),f12(124),fia(124)
     $/ 1.46000000e+01, 3.74071779e+01,-2.374823121522e-03/
      data eta(125),f12(125),fia(125)
     $/ 1.48000000e+01, 3.81724978e+01,-2.313510297832e-03/
      data eta(126),f12(126),fia(126)
     $/ 1.50000000e+01, 3.89430513e+01,-2.239491881527e-03/
      data eta(127),f12(127),fia(127)
     $/ 1.52000000e+01, 3.97187891e+01,-2.190399198766e-03/
      data eta(128),f12(128),fia(128)
     $/ 1.54000000e+01, 4.04996882e+01,-2.135558627693e-03/
      data eta(129),f12(129),fia(129)
     $/ 1.56000000e+01, 4.12857180e+01,-2.070549224394e-03/
      data eta(130),f12(130),fia(130)
     $/ 1.58000000e+01, 4.20768328e+01,-2.026371946726e-03/
      data eta(131),f12(131),fia(131)
     $/ 1.60000000e+01, 4.28730059e+01,-1.966932919251e-03/
      data eta(132),f12(132),fia(132)
     $/ 1.62000000e+01, 4.36741991e+01,-1.930467763422e-03/
      data eta(133),f12(133),fia(133)
     $/ 1.64000000e+01, 4.44803934e+01,-1.872060355250e-03/
      data eta(134),f12(134),fia(134)
     $/ 1.66000000e+01, 4.52915468e+01,-1.834576476007e-03/
      data eta(135),f12(135),fia(135)
     $/ 1.68000000e+01, 4.61076365e+01,-1.788034097328e-03/
      data eta(136),f12(136),fia(136)
     $/ 1.70000000e+01, 4.69286280e+01,-1.739134736556e-03/
      data eta(137),f12(137),fia(137)
     $/ 1.72000000e+01, 4.77544832e+01,-1.709653214469e-03/
      data eta(138),f12(138),fia(138)
     $/ 1.74000000e+01, 4.85851870e+01,-1.661786300207e-03/
      data eta(139),f12(139),fia(139)
     $/ 1.76000000e+01, 4.94207048e+01,-1.633815786441e-03/
      data eta(140),f12(140),fia(140)
     $/ 1.78000000e+01, 5.02610178e+01,-1.584346715242e-03/
      data eta(141),f12(141),fia(141)
     $/ 1.80000000e+01, 5.11060801e+01,-1.552916203350e-03/
      data eta(142),f12(142),fia(142)
     $/ 1.82000000e+01, 5.19558725e+01,-1.531255731001e-03/
      data eta(143),f12(143),fia(143)
     $/ 1.84000000e+01, 5.28103838e+01,-1.473403678003e-03/
      data eta(144),f12(144),fia(144)
     $/ 1.86000000e+01, 5.36695566e+01,-1.467474616122e-03/
      data eta(145),f12(145),fia(145)
     $/ 1.88000000e+01, 5.45334024e+01,-1.426838647689e-03/
      data eta(146),f12(146),fia(146)
     $/ 1.90000000e+01, 5.54018793e+01,-1.377498261708e-03/
      data eta(147),f12(147),fia(147)
     $/ 1.92000000e+01, 5.62749338e+01,-1.383183816947e-03/
      data eta(148),f12(148),fia(148)
     $/ 1.94000000e+01, 5.71525927e+01,-1.327524740849e-03/
      data eta(149),f12(149),fia(149)
     $/ 1.96000000e+01, 5.80347986e+01,-1.315576157266e-03/
      data eta(150),f12(150),fia(150)
     $/ 1.98000000e+01, 5.89215441e+01,-1.283351521326e-03/
      data eta(151),f12(151),fia(151)
     $/ 2.00000000e+01, 5.98127985e+01,-1.252219378154e-03/
      data eta(152),f12(152),fia(152)
     $/ 2.02000000e+01, 6.07085314e+01,-1.242160880550e-03/
      data eta(153),f12(153),fia(153)
     $/ 2.04000000e+01, 6.16087389e+01,-1.194983711474e-03/
      data eta(154),f12(154),fia(154)
     $/ 2.06000000e+01, 6.25133677e+01,-1.197567559790e-03/
      data eta(155),f12(155),fia(155)
     $/ 2.08000000e+01, 6.34224367e+01,-1.159493320310e-03/
      data eta(156),f12(156),fia(156)
     $/ 2.10000000e+01, 6.43359013e+01,-1.135131594020e-03/
      data eta(157),f12(157),fia(157)
     $/ 2.12000000e+01, 6.52537327e+01,-1.125982877246e-03/
      data eta(158),f12(158),fia(158)
     $/ 2.14000000e+01, 6.61759281e+01,-1.085955718327e-03/
      data eta(159),f12(159),fia(159)
     $/ 2.16000000e+01, 6.71024418e+01,-1.091435662477e-03/
      data eta(160),f12(160),fia(160)
     $/ 2.18000000e+01, 6.80332966e+01,-1.052359718121e-03/
      data eta(161),f12(161),fia(161)
     $/ 2.20000000e+01, 6.89684391e+01,-1.034522675630e-03/
      data eta(162),f12(162),fia(162)
     $/ 2.22000000e+01, 6.99078465e+01,-1.025329275552e-03/
      data eta(163),f12(163),fia(163)
     $/ 2.24000000e+01, 7.08515186e+01,-9.955965775667e-04/
      data eta(164),f12(164),fia(164)
     $/ 2.26000000e+01, 7.17994175e+01,-9.877545065490e-04/
      data eta(165),f12(165),fia(165)
     $/ 2.28000000e+01, 7.27515431e+01,-9.690227539283e-04/
      data eta(166),f12(166),fia(166)
     $/ 2.30000000e+01, 7.37078724e+01,-9.420188557201e-04/
      data eta(167),f12(167),fia(167)
     $/ 2.32000000e+01, 7.46683674e+01,-9.412102489942e-04/
      data eta(168),f12(168),fia(168)
     $/ 2.34000000e+01, 7.56330357e+01,-9.082185476753e-04/
      data eta(169),f12(169),fia(169)
     $/ 2.36000000e+01, 7.66018314e+01,-9.080175988176e-04/
      data eta(170),f12(170),fia(170)
     $/ 2.38000000e+01, 7.75747700e+01,-8.923514704106e-04/
      data eta(171),f12(171),fia(171)
     $/ 2.40000000e+01, 7.85518284e+01,-8.566171195270e-04/
      data eta(172),f12(172),fia(172)
     $/ 2.42000000e+01, 7.95329456e+01,-8.706095710928e-04/
      data eta(173),f12(173),fia(173)
     $/ 2.44000000e+01, 8.05181599e+01,-8.344314949249e-04/
      data eta(174),f12(174),fia(174)
     $/ 2.46000000e+01, 8.15074177e+01,-8.346049903250e-04/
      data eta(175),f12(175),fia(175)
     $/ 2.48000000e+01, 8.25007267e+01,-8.180174375162e-04/
      data eta(176),f12(176),fia(176)
     $/ 2.50000000e+01, 8.34980640e+01,-7.958462623517e-04/
      data eta(177),f12(177),fia(177)
     $/ 2.52000000e+01, 8.44993916e+01,-7.966215802492e-04/
      data eta(178),f12(178),fia(178)
     $/ 2.54000000e+01, 8.55047245e+01,-7.763841457450e-04/
      data eta(179),f12(179),fia(179)
     $/ 2.56000000e+01, 8.65140324e+01,-7.646832835317e-04/
      data eta(180),f12(180),fia(180)
     $/ 2.58000000e+01, 8.75272999e+01,-7.591181446306e-04/
      data eta(181),f12(181),fia(181)
     $/ 2.60000000e+01, 8.85445194e+01,-7.309023009314e-04/
      data eta(182),f12(182),fia(182)
     $/ 2.62000000e+01, 8.95656452e+01,-7.452741062596e-04/
      data eta(183),f12(183),fia(183)
     $/ 2.64000000e+01, 9.05907154e+01,-7.079420670479e-04/
      data eta(184),f12(184),fia(184)
     $/ 2.66000000e+01, 9.16196613e+01,-7.165112531131e-04/
      data eta(185),f12(185),fia(185)
     $/ 2.68000000e+01, 9.26525059e+01,-6.963409034097e-04/
      data eta(186),f12(186),fia(186)
     $/ 2.70000000e+01, 9.36892185e+01,-6.873393148553e-04/
      data eta(187),f12(187),fia(187)
     $/ 2.72000000e+01, 9.47297840e+01,-6.807802206753e-04/
      data eta(188),f12(188),fia(188)
     $/ 2.74000000e+01, 9.57741947e+01,-6.623413987856e-04/
      data eta(189),f12(189),fia(189)
     $/ 2.76000000e+01, 9.68224201e+01,-6.661242492164e-04/
      data eta(190),f12(190),fia(190)
     $/ 2.78000000e+01, 9.78744755e+01,-6.418176402524e-04/
      data eta(191),f12(191),fia(191)
     $/ 2.80000000e+01, 9.89303150e+01,-6.453620029091e-04/
      data eta(192),f12(192),fia(192)
     $/ 2.82000000e+01, 9.99899540e+01,-6.295277664925e-04/
      data eta(193),f12(193),fia(193)
     $/ 2.84000000e+01, 1.01053362e+02,-6.176912881028e-04/
      data eta(194),f12(194),fia(194)
     $/ 2.86000000e+01, 1.02120516e+02,-6.181417564908e-04/
      data eta(195),f12(195),fia(195)
     $/ 2.88000000e+01, 1.03191431e+02,-6.038839429712e-04/
      data eta(196),f12(196),fia(196)
     $/ 2.90000000e+01, 1.04266077e+02,-5.927623476398e-04/
      data eta(197),f12(197),fia(197)
     $/ 2.92000000e+01, 1.05344431e+02,-5.918987997347e-04/
      data eta(198),f12(198),fia(198)
     $/ 2.94000000e+01, 1.06426508e+02,-5.840707347615e-04/
      data eta(199),f12(199),fia(199)
     $/ 2.96000000e+01, 1.07512262e+02,-5.368048493408e-04/
      data eta(200),f12(200),fia(200)
     $/ 2.98000000e+01, 1.08601709e+02,-7.132081093603e-04/
      data eta(201),f12(201),fia(201)
     $/ 3.00000000e+01, 1.09694826e+02, 0.000000000000e+00/
c
c
c          if y is between 4.0234d-5 and 109.695 then use spline table
      if((y.le.109.695d0).and.(y.ge.4.0234d-5)) then
        call intrp(y,f1,f12,eta,fia,n,nlow,nhigh)
c          else if y is greater than 109.695 use degen. approximation
      elseif(y.gt.109.695d0) then
        x2=(1.5d0*y)**(0.666666667d0)  
        x4=1.0d0/(x2*x2) 
        f1=x2*(1.0d0+(ai(1)+(ai(2)+ai(3)*x4)*x4)*x4)
c          else if y is less than 4.0234d-5 use nondegen. approximation
      else
        f1=log(ai(7)*max(y,1.0d-20)*(1.d0+(ai(4)+(ai(5)+ai(8)*y)*y)*y))
      endif
c
      finv12=f1  
 999  return
      end   
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       fhalf
c    type:         double precision function
c    author:       f. douglas swesty, dpt of physics suny @ stony brook
c
c    call line:    fhalf(y)
c
c    inputs:       y (double precision)   (argument)
c
c    return:       ratio of 1/2th fermi integral to the -1/2th fermi
c                  integral (double precision)
c
c***********************************************************************
      double precision function fhalf(y)
      implicit none
c
      integer n, nlow, nhigh
      parameter(n=201)
      double precision y, a(7), eta(n), fr(n), fra(n), th, x2, f0, f1
c
c
c                retain values between calls
      save nlow,nhigh
c
      data a/6.16850274d0,1.77568655d0,6.92965606d0,.176776695d0,
     $ 6.41500299d-02,.4d0,1.32934039d0/
      data th,nlow,nhigh/0.3333333333333d0,1,n/
c
c
c                  cubic spline data
      data eta(  1),fr(  1),fra(  1)
     $/-1.00000000e+01, 5.00008047e-01, 0.000000000000e+00/
      data eta(  2),fr(  2),fra(  2)
     $/-9.80000000e+00, 5.00009809e-01, 1.160669214176e-05/
      data eta(  3),fr(  3),fra(  3)
     $/-9.60000000e+00, 5.00011969e-01, 1.327318909372e-05/
      data eta(  4),fr(  4),fra(  4)
     $/-9.40000000e+00, 5.00014639e-01, 1.170462092124e-05/
      data eta(  5),fr(  5),fra(  5)
     $/-9.20000000e+00, 5.00017836e-01, 1.896244399075e-05/
      data eta(  6),fr(  6),fra(  6)
     $/-9.00000000e+00, 5.00021791e-01, 2.619339170395e-05/
      data eta(  7),fr(  7),fra(  7)
     $/-8.80000000e+00, 5.00026692e-01, 1.825697414678e-05/
      data eta(  8),fr(  8),fra(  8)
     $/-8.60000000e+00, 5.00032512e-01, 3.857146281817e-05/
      data eta(  9),fr(  9),fra(  9)
     $/-8.40000000e+00, 5.00039730e-01, 3.706621670383e-05/
      data eta( 10),fr( 10),fra( 10)
     $/-8.20000000e+00, 5.00048524e-01, 4.966876577432e-05/
      data eta( 11),fr( 11),fra( 11)
     $/-8.00000000e+00, 5.00059289e-01, 5.986611168059e-05/
      data eta( 12),fr( 12),fra( 12)
     $/-7.80000000e+00, 5.00072436e-01, 6.815575901197e-05/
      data eta( 13),fr( 13),fra( 13)
     $/-7.60000000e+00, 5.00088423e-01, 9.353156612777e-05/
      data eta( 14),fr( 14),fra( 14)
     $/-7.40000000e+00, 5.00108035e-01, 1.015077511200e-04/
      data eta( 15),fr( 15),fra( 15)
     $/-7.20000000e+00, 5.00131901e-01, 1.384040424797e-04/
      data eta( 16),fr( 16),fra( 16)
     $/-7.00000000e+00, 5.00161156e-01, 1.532660352027e-04/
      data eta( 17),fr( 17),fra( 17)
     $/-6.80000000e+00, 5.00196763e-01, 2.014816386790e-04/
      data eta( 18),fr( 18),fra( 18)
     $/-6.60000000e+00, 5.00240334e-01, 2.353001371116e-04/
      data eta( 19),fr( 19),fra( 19)
     $/-6.40000000e+00, 5.00293492e-01, 2.953124554528e-04/
      data eta( 20),fr( 20),fra( 20)
     $/-6.20000000e+00, 5.00358457e-01, 3.546103132566e-04/
      data eta( 21),fr( 21),fra( 21)
     $/-6.00000000e+00, 5.00437760e-01, 4.368757029944e-04/
      data eta( 22),fr( 22),fra( 22)
     $/-5.80000000e+00, 5.00534601e-01, 5.285417635068e-04/
      data eta( 23),fr( 23),fra( 23)
     $/-5.60000000e+00, 5.00652775e-01, 6.489897161003e-04/
      data eta( 24),fr( 24),fra( 24)
     $/-5.40000000e+00, 5.00797055e-01, 7.913279951832e-04/
      data eta( 25),fr( 25),fra( 25)
     $/-5.20000000e+00, 5.00973167e-01, 9.605790184029e-04/
      data eta( 26),fr( 26),fra( 26)
     $/-5.00000000e+00, 5.01188039e-01, 1.180342230219e-03/
      data eta( 27),fr( 27),fra( 27)
     $/-4.80000000e+00, 5.01450297e-01, 1.425959520518e-03/
      data eta( 28),fr( 28),fra( 28)
     $/-4.60000000e+00, 5.01770117e-01, 1.750084178541e-03/
      data eta( 29),fr( 29),fra( 29)
     $/-4.40000000e+00, 5.02160252e-01, 2.120996973156e-03/
      data eta( 30),fr( 30),fra( 30)
     $/-4.20000000e+00, 5.02635834e-01, 2.582978496851e-03/
      data eta( 31),fr( 31),fra( 31)
     $/-4.00000000e+00, 5.03215376e-01, 3.140940050287e-03/
      data eta( 32),fr( 32),fra( 32)
     $/-3.80000000e+00, 5.03921295e-01, 3.809961012388e-03/
      data eta( 33),fr( 33),fra( 33)
     $/-3.60000000e+00, 5.04780575e-01, 4.623412578668e-03/
      data eta( 34),fr( 34),fra( 34)
     $/-3.40000000e+00, 5.05825854e-01, 5.596028708118e-03/
      data eta( 35),fr( 35),fra( 35)
     $/-3.20000000e+00, 5.07096250e-01, 6.760222613033e-03/
      data eta( 36),fr( 36),fra( 36)
     $/-3.00000000e+00, 5.08638542e-01, 8.147420929689e-03/
      data eta( 37),fr( 37),fra( 37)
     $/-2.80000000e+00, 5.10508521e-01, 9.803125317419e-03/
      data eta( 38),fr( 38),fra( 38)
     $/-2.60000000e+00, 5.12772505e-01, 1.174090923418e-02/
      data eta( 39),fr( 39),fra( 39)
     $/-2.40000000e+00, 5.15508351e-01, 1.401245248421e-02/
      data eta( 40),fr( 40),fra( 40)
     $/-2.20000000e+00, 5.18807170e-01, 1.665511284146e-02/
      data eta( 41),fr( 41),fra( 41)
     $/-2.00000000e+00, 5.22774696e-01, 1.967341586786e-02/
      data eta( 42),fr( 42),fra( 42)
     $/-1.80000000e+00, 5.27531885e-01, 2.310044392523e-02/
      data eta( 43),fr( 43),fra( 43)
     $/-1.60000000e+00, 5.33215712e-01, 2.692057340338e-02/
      data eta( 44),fr( 44),fra( 44)
     $/-1.40000000e+00, 5.39978788e-01, 3.110469072503e-02/
      data eta( 45),fr( 45),fra( 45)
     $/-1.20000000e+00, 5.47988069e-01, 3.559136204354e-02/
      data eta( 46),fr( 46),fra( 46)
     $/-1.00000000e+00, 5.57422393e-01, 4.028629703731e-02/
      data eta( 47),fr( 47),fra( 47)
     $/-8.00000000e-01, 5.68468657e-01, 4.505445242759e-02/
      data eta( 48),fr( 48),fra( 48)
     $/-6.00000000e-01, 5.81316396e-01, 4.971717382570e-02/
      data eta( 49),fr( 49),fra( 49)
     $/-4.00000000e-01, 5.96151066e-01, 5.411642817285e-02/
      data eta( 50),fr( 50),fra( 50)
     $/-2.00000000e-01, 6.13146980e-01, 5.800386974609e-02/
      data eta( 51),fr( 51),fra( 51)
     $/ 0.00000000e+00, 6.32458839e-01, 6.125961642856e-02/
      data eta( 52),fr( 52),fra( 52)
     $/ 2.00000000e-01, 6.54215519e-01, 6.368105824286e-02/
      data eta( 53),fr( 53),fra( 53)
     $/ 4.00000000e-01, 6.78513328e-01, 6.518526072255e-02/
      data eta( 54),fr( 54),fra( 54)
     $/ 6.00000000e-01, 7.05412292e-01, 6.575130141701e-02/
      data eta( 55),fr( 55),fra( 55)
     $/ 8.00000000e-01, 7.34934891e-01, 6.535476854033e-02/
      data eta( 56),fr( 56),fra( 56)
     $/ 1.00000000e+00, 7.67065889e-01, 6.408943410837e-02/
      data eta( 57),fr( 57),fra( 57)
     $/ 1.20000000e+00, 8.01755749e-01, 6.211687075523e-02/
      data eta( 58),fr( 58),fra( 58)
     $/ 1.40000000e+00, 8.38925803e-01, 5.947211293743e-02/
      data eta( 59),fr( 59),fra( 59)
     $/ 1.60000000e+00, 8.78471892e-01, 5.639992416069e-02/
      data eta( 60),fr( 60),fra( 60)
     $/ 1.80000000e+00, 9.20271987e-01, 5.302910851018e-02/
      data eta( 61),fr( 61),fra( 61)
     $/ 2.00000000e+00, 9.64191670e-01, 4.942180766085e-02/
      data eta( 62),fr( 62),fra( 62)
     $/ 2.20000000e+00, 1.01008812e+00, 4.579899792645e-02/
      data eta( 63),fr( 63),fra( 63)
     $/ 2.40000000e+00, 1.05781651e+00, 4.217259322424e-02/
      data eta( 64),fr( 64),fra( 64)
     $/ 2.60000000e+00, 1.10723239e+00, 3.863512430722e-02/
      data eta( 65),fr( 65),fra( 65)
     $/ 2.80000000e+00, 1.15819499e+00, 3.529307034029e-02/
      data eta( 66),fr( 66),fra( 66)
     $/ 3.00000000e+00, 1.21057018e+00, 3.208270748942e-02/
      data eta( 67),fr( 67),fra( 67)
     $/ 3.20000000e+00, 1.26423024e+00, 2.910594156410e-02/
      data eta( 68),fr( 68),fra( 68)
     $/ 3.40000000e+00, 1.31905603e+00, 2.635426039667e-02/
      data eta( 69),fr( 69),fra( 69)
     $/ 3.60000000e+00, 1.37493737e+00, 2.380801215631e-02/
      data eta( 70),fr( 70),fra( 70)
     $/ 3.80000000e+00, 1.43177253e+00, 2.148777375912e-02/
      data eta( 71),fr( 71),fra( 71)
     $/ 4.00000000e+00, 1.48946852e+00, 1.936381862584e-02/
      data eta( 72),fr( 72),fra( 72)
     $/ 4.20000000e+00, 1.54794050e+00, 1.745753296313e-02/
      data eta( 73),fr( 73),fra( 73)
     $/ 4.40000000e+00, 1.60711200e+00, 1.573191622765e-02/
      data eta( 74),fr( 74),fra( 74)
     $/ 4.60000000e+00, 1.66691366e+00, 1.413981428697e-02/
      data eta( 75),fr( 75),fra( 75)
     $/ 4.80000000e+00, 1.72728245e+00, 1.277958760511e-02/
      data eta( 76),fr( 76),fra( 76)
     $/ 5.00000000e+00, 1.78816286e+00, 1.148276396969e-02/
      data eta( 77),fr( 77),fra( 77)
     $/ 5.20000000e+00, 1.84950378e+00, 1.036724436250e-02/
      data eta( 78),fr( 78),fra( 78)
     $/ 5.40000000e+00, 1.91126027e+00, 9.382681356530e-03/
      data eta( 79),fr( 79),fra( 79)
     $/ 5.60000000e+00, 1.97339221e+00, 8.419828946720e-03/
      data eta( 80),fr( 80),fra( 80)
     $/ 5.80000000e+00, 2.03586235e+00, 7.668390275724e-03/
      data eta( 81),fr( 81),fra( 81)
     $/ 6.00000000e+00, 2.09863908e+00, 6.894960215311e-03/
      data eta( 82),fr( 82),fra( 82)
     $/ 6.20000000e+00, 2.16169262e+00, 6.274318363239e-03/
      data eta( 83),fr( 83),fra( 83)
     $/ 6.40000000e+00, 2.22499741e+00, 5.694183220240e-03/
      data eta( 84),fr( 84),fra( 84)
     $/ 6.60000000e+00, 2.28853030e+00, 5.163412244369e-03/
      data eta( 85),fr( 85),fra( 85)
     $/ 6.80000000e+00, 2.35227011e+00, 4.691055087970e-03/
      data eta( 86),fr( 86),fra( 86)
     $/ 7.00000000e+00, 2.41619827e+00, 4.323617268449e-03/
      data eta( 87),fr( 87),fra( 87)
     $/ 7.20000000e+00, 2.48029883e+00, 3.875965713102e-03/
      data eta( 88),fr( 88),fra( 88)
     $/ 7.40000000e+00, 2.54455567e+00, 3.614743308902e-03/
      data eta( 89),fr( 89),fra( 89)
     $/ 7.60000000e+00, 2.60895650e+00, 3.262132860466e-03/
      data eta( 90),fr( 90),fra( 90)
     $/ 7.80000000e+00, 2.67348849e+00, 3.011856988495e-03/
      data eta( 91),fr( 91),fra( 91)
     $/ 8.00000000e+00, 2.73814091e+00, 2.753550229065e-03/
      data eta( 92),fr( 92),fra( 92)
     $/ 8.20000000e+00, 2.80290404e+00, 2.582380523235e-03/
      data eta( 93),fr( 93),fra( 93)
     $/ 8.40000000e+00, 2.86776968e+00, 2.292083501014e-03/
      data eta( 94),fr( 94),fra( 94)
     $/ 8.60000000e+00, 2.93272843e+00, 2.215304435771e-03/
      data eta( 95),fr( 95),fra( 95)
     $/ 8.80000000e+00, 2.99777480e+00, 1.991168746407e-03/
      data eta( 96),fr( 96),fra( 96)
     $/ 9.00000000e+00, 3.06290146e+00, 1.862794325246e-03/
      data eta( 97),fr( 97),fra( 97)
     $/ 9.20000000e+00, 3.12810241e+00, 1.701422191793e-03/
      data eta( 98),fr( 98),fra( 98)
     $/ 9.40000000e+00, 3.19337213e+00, 1.646713595709e-03/
      data eta( 99),fr( 99),fra( 99)
     $/ 9.60000000e+00, 3.25870691e+00, 1.470861367159e-03/
      data eta(100),fr(100),fra(100)
     $/ 9.80000000e+00, 3.32410133e+00, 1.415936859359e-03/
      data eta(101),fr(101),fra(101)
     $/ 1.00000000e+01, 3.38955181e+00, 1.273721641867e-03/
      data eta(102),fr(102),fra(102)
     $/ 1.02000000e+01, 3.45505393e+00, 1.235581482298e-03/
      data eta(103),fr(103),fra(103)
     $/ 1.04000000e+01, 3.52060505e+00, 1.134274717011e-03/
      data eta(104),fr(104),fra(104)
     $/ 1.06000000e+01, 3.58620184e+00, 1.077069548628e-03/
      data eta(105),fr(105),fra(105)
     $/ 1.08000000e+01, 3.65184173e+00, 1.022991960613e-03/
      data eta(106),fr(106),fra(106)
     $/ 1.10000000e+01, 3.71752230e+00, 9.341272710622e-04/
      data eta(107),fr(107),fra(107)
     $/ 1.12000000e+01, 3.78324041e+00, 8.692110436287e-04/
      data eta(108),fr(108),fra(108)
     $/ 1.14000000e+01, 3.84899382e+00, 8.850650648729e-04/
      data eta(109),fr(109),fra(109)
     $/ 1.16000000e+01, 3.91478159e+00, 7.459549148300e-04/
      data eta(110),fr(110),fra(110)
     $/ 1.18000000e+01, 3.98060041e+00, 7.864194213587e-04/
      data eta(111),fr(111),fra(111)
     $/ 1.20000000e+01, 4.04644996e+00, 7.192672054582e-04/
      data eta(112),fr(112),fra(112)
     $/ 1.22000000e+01, 4.11232811e+00, 6.253942766733e-04/
      data eta(113),fr(113),fra(113)
     $/ 1.24000000e+01, 4.17823227e+00, 6.811530337942e-04/
      data eta(114),fr(114),fra(114)
     $/ 1.26000000e+01, 4.24416274e+00, 5.969604189414e-04/
      data eta(115),fr(115),fra(115)
     $/ 1.28000000e+01, 4.31011748e+00, 5.709059055339e-04/
      data eta(116),fr(116),fra(116)
     $/ 1.30000000e+01, 4.37609516e+00, 5.615397884449e-04/
      data eta(117),fr(117),fra(117)
     $/ 1.32000000e+01, 4.44209499e+00, 5.044814122883e-04/
      data eta(118),fr(118),fra(118)
     $/ 1.34000000e+01, 4.50811538e+00, 5.045695992625e-04/
      data eta(119),fr(119),fra(119)
     $/ 1.36000000e+01, 4.57415563e+00, 4.552324352849e-04/
      data eta(120),fr(120),fra(120)
     $/ 1.38000000e+01, 4.64021477e+00, 5.087014064024e-04/
      data eta(121),fr(121),fra(121)
     $/ 1.40000000e+01, 4.70629289e+00, 3.581955061927e-04/
      data eta(122),fr(122),fra(122)
     $/ 1.42000000e+01, 4.77238707e+00, 4.664301602782e-04/
      data eta(123),fr(123),fra(123)
     $/ 1.44000000e+01, 4.83849842e+00, 3.514279312654e-04/
      data eta(124),fr(124),fra(124)
     $/ 1.46000000e+01, 4.90462510e+00, 4.274941163720e-04/
      data eta(125),fr(125),fra(125)
     $/ 1.48000000e+01, 4.97076776e+00, 3.365813609877e-04/
      data eta(126),fr(126),fra(126)
     $/ 1.50000000e+01, 5.03692462e+00, 3.549707958916e-04/
      data eta(127),fr(127),fra(127)
     $/ 1.52000000e+01, 5.10309513e+00, 2.910214326093e-04/
      data eta(128),fr(128),fra(128)
     $/ 1.54000000e+01, 5.16927849e+00, 4.096124892769e-04/
      data eta(129),fr(129),fra(129)
     $/ 1.56000000e+01, 5.23547613e+00, 2.118142438439e-04/
      data eta(130),fr(130),fra(130)
     $/ 1.58000000e+01, 5.30168451e+00, 3.539701323030e-04/
      data eta(131),fr(131),fra(131)
     $/ 1.60000000e+01, 5.36790514e+00, 2.099303736157e-04/
      data eta(132),fr(132),fra(132)
     $/ 1.62000000e+01, 5.43413621e+00, 3.726817539688e-04/
      data eta(133),fr(133),fra(133)
     $/ 1.64000000e+01, 5.50037998e+00, 2.035344437886e-04/
      data eta(134),fr(134),fra(134)
     $/ 1.66000000e+01, 5.56663346e+00, 2.697503364877e-04/
      data eta(135),fr(135),fra(135)
     $/ 1.68000000e+01, 5.63289722e+00, 2.603197709474e-04/
      data eta(136),fr(136),fra(136)
     $/ 1.70000000e+01, 5.69917089e+00, 1.752530698713e-04/
      data eta(137),fr(137),fra(137)
     $/ 1.72000000e+01, 5.76545287e+00, 2.848424287755e-04/
      data eta(138),fr(138),fra(138)
     $/ 1.74000000e+01, 5.83174481e+00, 1.787223011498e-04/
      data eta(139),fr(139),fra(139)
     $/ 1.76000000e+01, 5.89804520e+00, 2.684866184639e-04/
      data eta(140),fr(140),fra(140)
     $/ 1.78000000e+01, 5.96435493e+00, 1.488933760423e-04/
      data eta(141),fr(141),fra(141)
     $/ 1.80000000e+01, 6.03067153e+00, 1.660321432524e-04/
      data eta(142),fr(142),fra(142)
     $/ 1.82000000e+01, 6.09699556e+00, 3.007372019870e-04/
      data eta(143),fr(143),fra(143)
     $/ 1.84000000e+01, 6.16332905e+00, 5.039185382223e-05/
      data eta(144),fr(144),fra(144)
     $/ 1.86000000e+01, 6.22966780e+00, 2.880169444751e-04/
      data eta(145),fr(145),fra(145)
     $/ 1.88000000e+01, 6.29601555e+00, 1.459654809555e-04/
      data eta(146),fr(146),fra(146)
     $/ 1.90000000e+01, 6.36236965e+00, 8.022559364040e-05/
      data eta(147),fr(147),fra(147)
     $/ 1.92000000e+01, 6.42872883e+00, 2.967073127773e-04/
      data eta(148),fr(148),fra(148)
     $/ 1.94000000e+01, 6.49509677e+00, 4.547040400385e-05/
      data eta(149),fr(149),fra(149)
     $/ 1.96000000e+01, 6.56146941e+00, 2.275888933824e-04/
      data eta(150),fr(150),fra(150)
     $/ 1.98000000e+01, 6.62784906e+00, 9.594149535815e-05/
      data eta(151),fr(151),fra(151)
     $/ 2.00000000e+01, 6.69423382e+00, 1.539878947448e-04/
      data eta(152),fr(152),fra(152)
     $/ 2.02000000e+01, 6.76062461e+00, 1.936619445003e-04/
      data eta(153),fr(153),fra(153)
     $/ 2.04000000e+01, 6.82702164e+00, 6.161813072681e-06/
      data eta(154),fr(154),fra(154)
     $/ 2.06000000e+01, 6.89342166e+00, 2.302557839252e-04/
      data eta(155),fr(155),fra(155)
     $/ 2.08000000e+01, 6.95982860e+00, 1.118128460033e-04/
      data eta(156),fr(156),fra(156)
     $/ 2.10000000e+01, 7.02624073e+00, 1.013932536890e-04/
      data eta(157),fr(157),fra(157)
     $/ 2.12000000e+01, 7.09265730e+00, 1.470875659220e-04/
      data eta(158),fr(158),fra(158)
     $/ 2.14000000e+01, 7.15907881e+00, 5.141533931037e-05/
      data eta(159),fr(159),fra(159)
     $/ 2.16000000e+01, 7.22550396e+00, 1.950632132465e-04/
      data eta(160),fr(160),fra(160)
     $/ 2.18000000e+01, 7.29193509e+00, 6.331854420200e-05/
      data eta(161),fr(161),fra(161)
     $/ 2.20000000e+01, 7.35836991e+00, 1.054146921341e-04/
      data eta(162),fr(162),fra(162)
     $/ 2.22000000e+01, 7.42480867e+00, 1.072633188645e-04/
      data eta(163),fr(163),fra(163)
     $/ 2.24000000e+01, 7.49125168e+00, 1.027449618964e-04/
      data eta(164),fr(164),fra(164)
     $/ 2.26000000e+01, 7.55769894e+00, 1.183006063983e-04/
      data eta(165),fr(165),fra(165)
     $/ 2.28000000e+01, 7.62415059e+00, 8.253924684004e-05/
      data eta(166),fr(166),fra(166)
     $/ 2.30000000e+01, 7.69060579e+00, 8.451551828560e-05/
      data eta(167),fr(167),fra(167)
     $/ 2.32000000e+01, 7.75706437e+00, 8.685419149568e-05/
      data eta(168),fr(168),fra(168)
     $/ 2.34000000e+01, 7.82352627e+00, 6.485775084105e-05/
      data eta(169),fr(169),fra(169)
     $/ 2.36000000e+01, 7.88999116e+00, 1.043033206627e-04/
      data eta(170),fr(170),fra(170)
     $/ 2.38000000e+01, 7.95646040e+00, 1.690790286878e-04/
      data eta(171),fr(171),fra(171)
     $/ 2.40000000e+01, 8.02293437e+00,-7.229290158689e-05/
      data eta(172),fr(172),fra(172)
     $/ 2.42000000e+01, 8.08940876e+00, 1.853851230301e-04/
      data eta(173),fr(173),fra(173)
     $/ 2.44000000e+01, 8.15588780e+00, 2.682857763816e-05/
      data eta(174),fr(174),fra(174)
     $/ 2.46000000e+01, 8.22236941e+00, 9.315857978144e-05/
      data eta(175),fr(175),fra(175)
     $/ 2.48000000e+01, 8.28885454e+00, 1.279132160594e-04/
      data eta(176),fr(176),fra(176)
     $/ 2.50000000e+01, 8.35534335e+00,-5.151623010689e-05/
      data eta(177),fr(177),fra(177)
     $/ 2.52000000e+01, 8.42183289e+00, 1.870374120185e-04/
      data eta(178),fr(178),fra(178)
     $/ 2.54000000e+01, 8.48832695e+00,-1.961965051994e-05/
      data eta(179),fr(179),fra(179)
     $/ 2.56000000e+01, 8.55482254e+00, 1.229752535397e-04/
      data eta(180),fr(180),fra(180)
     $/ 2.58000000e+01, 8.62132190e+00, 9.217292731240e-05/
      data eta(181),fr(181),fra(181)
     $/ 2.60000000e+01, 8.68782416e+00,-5.687752100593e-05/
      data eta(182),fr(182),fra(182)
     $/ 2.62000000e+01, 8.75432688e+00, 2.031850501482e-04/
      data eta(183),fr(183),fra(183)
     $/ 2.64000000e+01, 8.82083439e+00,-3.579723391714e-05/
      data eta(184),fr(184),fra(184)
     $/ 2.66000000e+01, 8.88734279e+00, 7.422832091008e-05/
      data eta(185),fr(185),fra(185)
     $/ 2.68000000e+01, 8.95385391e+00, 1.447580583760e-04/
      data eta(186),fr(186),fra(186)
     $/ 2.70000000e+01, 9.02036919e+00,-2.735083889056e-05/
      data eta(187),fr(187),fra(187)
     $/ 2.72000000e+01, 9.08688502e+00, 4.542808910866e-05/
      data eta(188),fr(188),fra(188)
     $/ 2.74000000e+01, 9.15340256e+00, 1.035480068597e-04/
      data eta(189),fr(189),fra(189)
     $/ 2.76000000e+01, 9.21992316e+00,-6.679452567842e-07/
      data eta(190),fr(190),fra(190)
     $/ 2.78000000e+01, 9.28644514e+00, 1.046883666256e-04/
      data eta(191),fr(191),fra(191)
     $/ 2.80000000e+01, 9.35297006e+00, 2.438337814088e-05/
      data eta(192),fr(192),fra(192)
     $/ 2.82000000e+01, 9.41949693e+00, 8.980690836273e-05/
      data eta(193),fr(193),fra(193)
     $/ 2.84000000e+01, 9.48602609e+00,-3.931880271592e-05/
      data eta(194),fr(194),fra(194)
     $/ 2.86000000e+01, 9.55255558e+00, 1.158160177556e-04/
      data eta(195),fr(195),fra(195)
     $/ 2.88000000e+01, 9.61908818e+00, 4.256278618984e-05/
      data eta(196),fr(196),fra(196)
     $/ 2.90000000e+01, 9.68562279e+00, 1.531962767634e-05/
      data eta(197),fr(197),fra(197)
     $/ 2.92000000e+01, 9.75215876e+00, 1.018906708333e-04/
      data eta(198),fr(198),fra(198)
     $/ 2.94000000e+01, 9.81869734e+00,-3.387023032781e-05/
      data eta(199),fr(199),fra(199)
     $/ 2.96000000e+01, 9.88523639e+00, 1.063057625796e-04/
      data eta(200),fr(200),fra(200)
     $/ 2.98000000e+01, 9.95177816e+00, 1.564710994435e-05/
      data eta(201),fr(201),fra(201)
     $/ 3.00000000e+01, 1.00183211e+01, 0.000000000000e+00/
c
c          if y is between -10 and 30 then use the spline tables
      if((y.ge.-10.0d0).and.(y.le.30.0d0)) then
        call intrp(y,f1,eta,fr,fra,n,nlow,nhigh)
c          else if y is greater than 30 use degenerate approximation
      elseif(y.gt.30.0d0) then
        x2=y**(-2)
        f1=y*th*(1.0d0+(0.2d0*a(1)+(0.6d0*a(2)+1.4d0*x2*a(3))*x2)*x2)
     >  /(1.0d0-(0.2d0*th*a(1)-(a(2)-4.2d0*x2*a(3))*x2)*x2)  
      else
c          else if y is less than -10 use nondegenerate approximation
        f0=exp(y)   
        f1=(1.0d0-(2.0d0*a(4)-(3*a(5)-0.125d0*f0)*f0)*f0)/
     >  (2.0d0-(8.0d0*a(4)-(18.0d0*a(5)-f0)*f0)*f0)
      endif
c
      fhalf=f1  
 999  return
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
      subroutine intrp(x,y,xtable,ytable,sc,npts,nlow,nhigh)
c
      implicit none
c
      integer npts
      double precision x, y, xtable(npts), ytable(npts), sc(npts)
c
c
      double precision c1, c2, s1, s2, deltx, sixth
      parameter(sixth=1.0d0/6.0d0)
      integer nlow, nhigh, nmid
c
      nhigh = max(min(npts-1,nhigh),2)
      nlow = max(min(npts-1,nlow),2)
c
c             if old values of bracketing pointers work then use them...
      if((nhigh-nlow).eq.1) then
        if((xtable(nhigh).gt.x).and.(xtable(nlow).le.x)) then
          goto 20
c                 try the next lower pair
        elseif((xtable(nhigh-1).gt.x).and.(xtable(nlow-1).le.x)) then
          nhigh = max(nhigh-1,1)
          nlow = max(nlow-1,1)
          goto 20
c                 try the next higher pair
        elseif((xtable(nhigh+1).gt.x).and.(xtable(nlow+1).le.x)) then
          nhigh = min(nhigh+1,npts)
          nlow = min(nlow+1,npts)
          goto 20
        endif          
      endif
c
c       otherwise find the ordinates that bracket the x value by
c       bisecting the table...
c
c                         set the high & low pointers to their extrema
      nlow = 1
      nhigh = npts
c
 10   nmid = (nhigh+nlow)/2
      if(x.gt.xtable(nmid)) then
        nlow = nmid
      else
        nhigh = nmid
      endif
      if((nhigh-nlow).ne.1) goto 10 
c
c
c                  now that the bracketing ordinates have been found,
c                  interpolate to get the y value
c
 20   continue
      if((nhigh.gt.npts).or.(nhigh.lt.1).or.
     *   (nlow.ge.nhigh).or.(nlow.lt.1) ) then
        write(*,*) ' intrp: nh,nl = ',nhigh,nlow,npts
      endif
c                  distance between nearest x entries
      deltx = xtable(nhigh)-xtable(nlow)
c
c                  now evaluate the polynomial...
c
c                        calculate linear interp. coefficients
      c1 = (xtable(nhigh)-x)/deltx
      c2 = (x-xtable(nlow))/deltx
c
c                        calculate spline interp. coefficients
      s1 = (c1**3)-c1
      s2 = (c2**3)-c2
c
c                        calculate spline interp. value
      y = c1*ytable(nlow)+c2*ytable(nhigh)+
     1  (s1*sc(nlow)+s2*sc(nhigh))*(deltx**2)*sixth
c
c                        if interpolated value of y is not monotonic
c                        then linearly interpolate
      if((y.lt.ytable(nlow)).or.(y.gt.ytable(nhigh))) then
        y = c1*ytable(nlow)+c2*ytable(nhigh)
      endif
c
      return
      end





cc JCM: moved any_electron() to jon_lsbox.f
      subroutine any_electron_old(tin,yein,dnsin) 
      implicit none

      include 'eos_m4c.inc'
      include 'el_eos.inc'

      include 'vector_eos.dek'

c..declare the pass
      double precision tin,yein,dnsin



c..local variables
      double precision ytot1,zbarxx,abar,zbar,yelocal

c..for unit conversions
      double precision fm,fm3,avo,kerg,kev,k2mev,mev2erg,
     1                 clight,light2,me,mecc
      parameter        (fm      = 1.0d-13,
     1                  fm3     = fm*fm*fm,
     2                  avo     = 6.0221367d23,
     3                  kerg    = 1.380658d-16,
     4                  kev     = 8.617385d-5,
     5                  k2mev   = kev * 1.0d-6,
     6                  mev2erg = 1.602177d-6,
     7                  clight  = 2.99792458d10,
     8                  light2  = clight*clight,
     9                  me      = 9.1093897d-28,
     &                  mecc    = me*light2)


c..debug formats
 111  format(1x,1p3e24.16)
 112  format(1x,1p6e14.6)



c..the original e+e- eos in ls
c      call el_eos(tin,yein,dnsin)



c..some better ones

c..set the input vector in cgs units
      temp_row(1) = tin / k2mev
      den_row(1)  = dnsin / (avo*fm3)

c..abar and zbar of the whole mixture
      ytot1   = xnut  + xprot + 0.25d0*xalfa 
      if (a .ne. 0.0) ytot1 = ytot1 + xh/a
      zbarxx  = xprot + 0.5d0*xalfa + x*xh
      abar    = 1.0d0/ytot1
c JCM:
c      abarnum=abar
      zbar    = zbarxx * abar
      yelocal = zbar/abar

      abar_row(1) = abar
c JCM:
      abarnum_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1

c..do it
      call eosfxt


c..set the output vector
      nsube  = xnem_row(1) * fm3
      neplus = xnp_row(1) * fm3
      musube = (etaele_row(1) * kerg*temp_row(1) + mecc)/mev2erg

c..electron thermodynamics
      epress = (pele_row(1) + ppos_row(1)) * fm3/mev2erg
      eu     = (eele_row(1) + epos_row(1) + mecc*avo*yelocal)
     1          / (avo * mev2erg)
      es     = (sele_row(1) + spos_row(1)) / (kerg*avo)
      fsube  = eu - tin * es

c..photon thermodynamics
      ppress = prad_row(1) * fm3/mev2erg
      pu     = erad_row(1) / (avo * mev2erg)
      ps     = srad_row(1) / (kerg*avo)
      pf     = pu - tin * ps


c..chemical potential derivaties
      demudt = (detat_row(1)*temp_row(1) + etaele_row(1)) * kerg
     1         / (mev2erg*k2mev)
      demudn = detad_row(1)*kerg*temp_row(1) / mev2erg / (avo*fm3)
      demudy = (detaz_row(1) - detaa_row(1)/yelocal)
     1         * abar_row(1)*kerg*temp_row(1) * yelocal / mev2erg

c..electron derivatives
      depdt  = dpept_row(1) * fm3/mev2erg / k2mev
      depdn  = dpepd_row(1) * fm3/mev2erg / (avo*fm3)
      depdy  = (dpepz_row(1) - dpepa_row(1)/yelocal)
     1          * abar_row(1) * yelocal * fm3/mev2erg

      deudt  = deept_row(1) / (avo * mev2erg) / k2mev
      deudn  = deepd_row(1) /(avo * mev2erg) / (avo*fm3)
      deudy  = ((deepz_row(1) - deepa_row(1)/yelocal)
     1           *abar_row(1)*yelocal  + mecc*avo)
     2          /(avo * mev2erg)

      desdt  = dsept_row(1) / (kerg*avo) / k2mev
      desdn  = dsepd_row(1) / (kerg*avo) / (avo*fm3)
      desdy  = (dsepz_row(1) -dsepa_row(1)/yelocal)*abar_row(1)*yelocal 
     1         / (kerg*avo) 

c..photon derivatives
      dppdt  = dpradt_row(1) * fm3/mev2erg / k2mev
      dppdn  = dpradd_row(1) * fm3/mev2erg / (avo*fm3)
      dppdy  = (dpradz_row(1) - dprada_row(1)/yelocal)
     1          *abar_row(1)*yelocal * fm3/mev2erg 

      dpudt  = deradt_row(1) / (avo * mev2erg) / k2mev
      dpudn  = deradd_row(1) / (avo * mev2erg) / (avo*fm3)
      dpudy  = (deradz_row(1) - derada_row(1)/yelocal)
     1         *abar_row(1)*yelocal/ (avo * mev2erg)

      dpsdt  = dsradt_row(1) / (kerg*avo) / k2mev
      dpsdn  = dsradd_row(1) / (kerg*avo) / (avo*fm3)
      dpsdy  = (dsradz_row(1) - dsrada_row(1)/yelocal)
     1         *abar_row(1)*yelocal / (kerg*avo) 

      return
      end






c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c                  this file contains the complete electron eos
c                  package el_eos_pak.f.  also needed is the include
c                  file el_eos.inc.  this complete set of routines
c                  is included in the ls eos code as of version 2.7
c
c                  the routines contained in this package are:
c
c                  el_eos, el_rel, el_imt, imtrule, gl16, & f2
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         el_eos.for
c
c***********************************************************************
c
c    module:       el_eos
c    type:         subroutine
c    author:       f. douglas swesty, 
c                  laboratory for computational astrophysics
c                  dept. of astronomy & ncsa
c                  university of illinois at urbana-champaign
c
c    email:        dswesty@ncsa.uiuc.edu
c
c    date:         (original) 2/12/91
c                  (v2.7)     9/15/95
c
c    purpose:      the electron and photon equation of state
c
c
c    call line:    call el_eos(t,ye,brydns)
c
c    inputs:       t = temperature
c                  ye = electron fraction
c                  brydns = baryon number density
c
c    outputs:      none
c
c
c 
c    include files:  el_eos.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine el_eos(t,ye,brydns)
c
      implicit none
c
      double precision t, ye, brydns
c
      include 'el_eos.inc'
c
c
c
c                           plancks constant & speed of light
      double precision hbar, c
      parameter (hbar=6.58217317d-22,c=2.997924581d23)
c
c                           pi and 1/3
      double precision pi, pi2, ovr3, movr3, ovr23
      parameter(pi=3.1415927, pi2=9.8696044)
      parameter(ovr3=0.33333333, movr3=-0.33333333, ovr23=0.66666667)
c
c
      if( ((brydns*ye).gt.1.0d-4).or.(t.gt.5.0d0) ) then
c                      assume that the electrons are relativistic
        call el_rel(t,ye,brydns)
      else
c                      otherwise call the gauss-laguerre version
        call el_imt(t,ye,brydns)
c                      if the stuff is really degenerate call the
c                      relativistic version anyway
        if((musube/t).gt.2.0d3) then
          call el_rel(t,ye,brydns)
        endif
      endif
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
 999  return
c
c
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         el_rel.for
c
c***********************************************************************
c
c    module:       el_rel
c
c    type:         subroutine
c
c    author:       f. douglas swesty
c                  laboratory for computational astrophysics
c                  dept. of astronomy & ncsa
c                  university of illinois at urbana-champaign
c
c    email:        dswesty@ncsa.uiuc.edu
c
c    date:         2/12/91
c
c
c    purpose:      the relativistic electron and photon eos
c
c
c    call line:    call el_rel(t,ye,brydns)
c
c    inputs:       t = temperature
c                  ye = electron fraction
c                  brydns = baryon number density
c
c    outputs:      none
c
c
c 
c    include files:  el_rel.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine el_rel(t,ye,brydns)
c
      implicit none
c
      double precision t, ye, brydns
c
      include 'el_eos.inc'
c
c
c
c                           plancks constant & speed of light
      double precision hbar, c
      parameter (hbar=6.58217317d-22,c=2.997924581d23)
c
c                           pi and 1/3
      double precision pi, pi2, ovr3, movr3, ovr23
      parameter(pi=3.1415927, pi2=9.8696044)
      parameter(ovr3=0.33333333, movr3=-0.33333333, ovr23=0.66666667)
c
c                           2nd fermi integral
      double precision f_2
c
c                           positron degeneracy parameter
      double precision elpeta
c
c
c                    leptons
c
c                    electron number density
      nsube = brydns*ye
c
c                    coefficants for chemical potential
c                    and thermodynamics quantities
      qsube = 1.0/( 3.0*(pi**2)*((hbar*c)**3) )
c
      acoef = 0.5*nsube/qsube
c
      bcoef = (acoef**2+((pi**6)*t**6)/27.0)**0.5
c
      dbdt = (pi**6)*(t**5)/(9.0*bcoef)
c
      ccoef = (acoef+bcoef)**ovr3
c
c
c                    electron chemical potential
      musube = ccoef-ovr3*((pi*t)**2)/ccoef
c
c                    positron degeneracy parameter
      elpeta = -musube/t
c
c                    positron number density
      neplus =  3.0*qsube*(t**3)*f_2(elpeta)
c
c
c                    electron pressure for rel. case
      epress = 0.25*qsube*(musube**4+2.0*(pi*t*musube)**2+
     1 7.0*((pi*t)**4)/15.0)
c
c
c                    electron internal energy per baryon
      eu = 0.75*qsube*(musube**4+2.0*(pi*musube*t)**2+
     1 7.0*((pi*t)**4)/15.0)/brydns
c
c
c                    electron free energy per baryon
      fsube = ((musube*nsube)-epress)/brydns
c
c                    electron entropy per baryon
      es = qsube*(((pi*musube)**2)*t+7.0*(pi**4)*(t**3)/
     1 15.0)/brydns
c
c                    photons
c
c                    photon pressure
      ppress = (pi**2)*(t**4)/(45.0*((hbar*c)**3))
c                    photon entropy per baryon
      ps = 4.0*ppress/(t*brydns)
c
c                    photon internal energy per baryon
      pu = 3.0*ppress/brydns
c
c                    photon free energy per baryon
      pf = pu-t*ps
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c                    derivatives of chem. potential w.r.t. t,
c                    brydns, ye
c
      demudt = dbdt/(3.0*ccoef**2)-ovr23*(pi**2)*t/ccoef+
     1         dbdt*((pi*t)**2)/(9.0*ccoef**4)
c
      demudn = (ye*pi2*(hbar*c)**3)/(musube**2+ovr3*pi2*t**2)
c
      demudy = brydns*demudn/ye
c
c
c                    derivatives of pressure w.r.t. brydns,ye,t
c
      depdn = brydns*ye*demudn
c
      depdy = brydns*depdn/ye
c
      depdt = brydns*(es+ye*demudt)
c
c
c                    derivatives of entropy w.r.t. t,brydns,ye
c
      desdt = es/t+ovr23*(7.0*pi2*(t**2)/15.0+musube*t*demudt)/
     1        (brydns*(hbar*c)**3)
c
      desdn = -1.0*depdt/(brydns**2)
c
      desdy = 2.0*t*qsube*pi2*musube*demudy/brydns
c
c
c                    derivatives of internal energy w.r.t.
c                    t,brydns,ye
      deudt = t*desdt
c
      deudn = (ye*(musube-t*demudt)-eu)/brydns
c
      deudy = 3.0*qsube*((musube**3)+pi2*(t**2)*musube)*
     1        demudy/brydns
c
c
c                               photons
c
c                    derivatives of photon pressure
      dppdn = 0.0
      dppdt = brydns*ps
      dppdy = 0.0
c
c                    derivatives of photon entropy
      dpsdn = -ps/brydns
      dpsdt = 3.0*ps/t
      dpsdy = 0.0
c
c                    derivatives of internal energy
      dpudn = -0.75*t*ps/brydns
      dpudt = 3.0*ps
      dpudy = 0.0
c
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
 999  return
c
c
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    file:         el_imt.for
c
c***********************************************************************
c
c    module:       el_imt
c
c    type:         subroutine
c
c    author:       f. douglas swesty
c                  laboratory for computational astrophysics
c                  dept. of astronomy & ncsa
c                  university of illinois at urbana-champaign
c
c    email:        dswesty@ncsa.uiuc.edu
c
c    date:         9/15/95
c
c    purpose:      the elctron and photon equation of state
c
c
c    call line:    call el_imt(t,ye,brydns)
c
c    inputs:       t = temperature
c                  ye = electron fraction
c                  brydns = baryon number density
c
c    outputs:      none
c
c
c 
c    include files:  el_eos.inc
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine el_imt(t,ye,brydns)
c
      implicit none
c
      double precision t, ye, brydns
c
      include 'el_eos.inc'
c
c
c
c                           plancks constant & speed of light
      double precision hbar, c
      parameter (hbar=6.58217317d-22,c=2.997924581d23)
      double precision mass_e
      parameter(mass_e=0.51100d0)
c
c                           pi and 1/3
      double precision pi, pi2, ovr3, movr3, ovr23
      parameter(pi=3.1415927, pi2=9.8696044)
      parameter(ovr3=0.33333333, movr3=-0.33333333, ovr23=0.66666667)
c
c                           2nd fermi integral
      double precision f_2
c                           positron degeneracy parameter
      double precision elpeta
c                           multiplicative factors in the g-l
c                           quadrature
      double precision tfac, dtfdt, efac, pfac
c                           electron eta and beta
      double precision eta_e, beta, dbtdt, detdt, detdne, d_eta
c                           number density and it's derivatives
      double precision ne_chk, nen, dnendb, dnende
      double precision nem, dnemdb, dnemde
      double precision nep, dnepdb, dnepde
c                           energy density & derivatives
      double precision e_eng, dedb, dede
      double precision em_eng, demdb, demde
      double precision ep_eng, depdb, depde
c                           pressure & derivatives
      double precision e_pr, dpdb, dpde
      double precision pm, dpmdb, dpmde
      double precision pp, dppdb, dppde
c
c                           loop variables
      integer j, maxitr
      parameter(maxitr=60)
c
c                           convergence tolerance for n-r iteration
      double precision epsil
      parameter(epsil=1.0d-10)
c
      double precision x_lo, x_up, eta
c
c-----------------------------------------------------------------------
c
c
c
c
c                    electron number density
      nsube = brydns*ye

c
c                    coefficants for chemical potential
c                    and thermodynamics quantities
      qsube = 1.0d0/( 3.0d0*(pi**2)*((hbar*c)**3) )
c
      acoef = 0.5d0*nsube/qsube
c
      bcoef = (acoef**2+((pi**6)*t**6)/27.0d0)**0.5d0
c
      dbdt = (pi**6)*(t**5)/(9.0d0*bcoef)
c
      ccoef = (acoef+bcoef)**ovr3
c
c
c                    electron chemical potential
      musube = ccoef-ovr3*((pi*t)**2)/ccoef
c
c
c
c                            initial guess at electron degeneracy 
c                            parameter
      eta_e = musube/t
c
      beta = mass_e/t
      dbtdt = -beta/t
c                            multiplicative factor in fermi integrals
      tfac = (t**3)/((pi**2)*(hbar*c)**3)
      dtfdt = 3.0d0*tfac/t
      pfac = tfac*t*ovr3
      efac = tfac*t/brydns
c
c                            loop until the n-r iteration for the
c                            chemical potential converges
      do 20 j=1,maxitr,1

c JCM:
c         write(*,*) 'loopj',j,maxitr
c
c                            zero the electron number density integral
        nem = 0.0d0
        dnemde = 0.0d0
        dnemdb = 0.0d0
        nep = 0.0d0
        dnepde = 0.0d0
        dnepdb = 0.0d0
c
c                            zero the electron energy density integral
        em_eng = 0.0d0
        demde = 0.0d0
        demdb = 0.0d0
        ep_eng = 0.0d0
        depde = 0.0d0
        depdb = 0.0d0
c
c                            zero the electron pressure integral
        pm = 0.0d0
        dpmde = 0.0d0
        dpmdb = 0.0d0
        pp = 0.0d0
        dppde = 0.0d0
        dppdb = 0.0d0
c
c                            zero the positron number density
        neplus = 0.0d0
c
c       ----------------------------------------------------------------
c       |  do the quadrature for electrons first                       |
c       ----------------------------------------------------------------
        eta = eta_e
        if(eta.gt.5.0d-1) then
c       |                    do the imt quadrature in two parts:       |
c       |                    (0,eta) & (eta,eta+60)                    |
          x_lo = 0.0d0
          x_up = eta
          call imtrule(eta,beta,x_lo,x_up,nem,dnemdb,dnemde,
     *                 em_eng,demdb,demde,pm,dpmdb,dpmde)
          x_lo = eta
          x_up = eta+60.0d0
          call imtrule(eta,beta,x_lo,x_up,nem,dnemdb,dnemde,
     *                 em_eng,demdb,demde,pm,dpmdb,dpmde)
c       |                                                              |
        else                                                           !
c       ----------------------------------------------------------------
c       |                    do the imt quadrature in 1 part           |
c       ----------------------------------------------------------------
          x_lo = 0.0d0
          x_up = 60.0d0
          call imtrule(eta,beta,x_lo,x_up,nem,dnemdb,dnemde,
     *                 em_eng,demdb,demde,pm,dpmdb,dpmde)
cc          call gl16(eta,beta,x_lo,x_up,nem,dnemdb,dnemde,
cc     *                 em_eng,demdb,demde,pm,dpmdb,dpmde)
c       |                                                              |
        endif
c       ----------------------------------------------------------------
c       |               now do the quadrature for positrons            |
c       ----------------------------------------------------------------
        eta = -eta_e
        if(eta.gt.5.0d-1) then
c       |                    do the imt quadrature in two parts:       |
c       |                    (0,eta) & (eta,eta+60)                    |
c       |                                                              |
          x_lo = 0.0d0
          x_up = eta
          call imtrule(eta,beta,x_lo,x_up,nep,dnepdb,dnepde,
     *                 ep_eng,depdb,depde,pp,dppdb,dppde)
          x_lo = eta
          x_up = eta+60.0d0
          call imtrule(eta,beta,x_lo,x_up,nep,dnepdb,dnepde,
     *                 ep_eng,depdb,depde,pp,dppdb,dppde)
c       |                                                              |
        else
          x_lo = 0.0d0
          x_up = 60.0d0
          call imtrule(eta,beta,x_lo,x_up,nep,dnepdb,dnepde,
     *                 ep_eng,depdb,depde,pp,dppdb,dppde)
cc          call gl16(eta,beta,x_lo,x_up,nep,dnepdb,dnepde,
cc     *                 ep_eng,depdb,depde,pp,dppdb,dppde)
        endif
c       |                                                              |
c       |                                                              |
c       |                                                              |
c       |                                                              |
c       ----------------------------------------------------------------
c
c                    now sum the electron & positron contributions
        nen = nem-nep
        dnende = dnemde+dnepde
        dnendb = dnemdb-dnepdb
        e_eng = em_eng+ep_eng
        dede = demde-depde
        dedb = demdb+depdb
        e_pr = pm+pp
        dpde = dpmde-dppde
        dpdb = dpmdb+dppdb
        neplus = nep
c
c                            multiply by the temperature factor to get
c                            number density
        ne_chk = tfac*nen
c                            multiply by the temperature factor to get
c                            number density
        neplus = tfac*neplus
c
c                            calculate the new change in eta
        d_eta = -(ne_chk-nsube)/(dnende*tfac)
c
c                   if we've met the convergence criterion...
        if(dabs(d_eta).lt.(epsil*eta_e)) then
c                          then break out the chemical potential
c                          loop
          goto 30
        else
c                          otherwise update the chemical potential
          eta_e = eta_e+d_eta
        endif
c
 20   continue
c
c                                      if we reached this point the
c                                      n-r iteration didn't converge
      write(*,*) ' el_eos: n-r iteration didnt converge!'
      write(*,*) ' el_eos: try = ',t,brydns,ye
      stop
c
c               if we reached this point then the n-r iteration has
c               converged.
 30   continue
c
c                     is the result consistent with electron number
c                     density?
      if(dabs((ne_chk-nsube)/nsube).gt.1.0d-5) then
        write(*,*) ' el_eos: gaussian quadrature converged badly!'
        stop
      endif
c
c             calculate thermodynamic quantities...
c
c                          electron chemical potential
      musube = t*eta_e
c                          electron pressure
      epress = pfac*e_pr
c                          electron internal energy per baryon
      eu = efac*e_eng
c                          electron free energy per baryon
      fsube = ye*musube-epress/brydns
c                          electron entropy per baryon
      es = (eu-fsube)/t
c
c                          derivative of the electron eta w.r.t. t
      detdt = -(dtfdt*nen+tfac*dnendb*dbtdt)/(tfac*dnende)
c                          derivative of the electron eta w.r.t. nsube
      detdne = 1.0d0/(tfac*dnende)
c
c                    derivatives of chem. potential w.r.t. t,
c                    brydns, ye
      demudt = t*detdt+musube/t
      demudn = ye*t*detdne
      demudy = brydns*t*detdne
c
c
c                    derivatives of pressure w.r.t. brydns,ye,t
      depdn = ye*pfac*dpde*detdne
      depdy = brydns*pfac*dpde*detdne
      depdt = (4.0d0*pfac/t)*e_pr+pfac*(dpde*detdt+dpdb*dbtdt)
c
c                    derivatives of internal energy w.r.t.
c                    t,brydns,ye
      deudt = (4.0d0*efac/t)*e_eng+efac*(dedb*dbtdt+dede*detdt)
      deudn = (-1.0d0*efac/brydns)*e_eng+ye*efac*dede*detdne
      deudy = brydns*efac*dede*detdne
c
c                    derivatives of entropy w.r.t. t,brydns,ye
      desdt = -es/t+(deudt-ye*demudt+depdt/brydns)/t
      desdn = (deudn-ye*demudn+depdn/brydns-epress/(brydns**2))/t
      desdy = (deudy-musube-ye*demudy+depdy/brydns)/t
c
c
c
c                    photon pressure
      ppress = (pi**2)*(t**4)/(45.0*((hbar*c)**3))
c                    photon entropy per baryon
      ps = 4.0*ppress/(t*brydns)
c
c                    photon internal energy per baryon
      pu = 3.0*ppress/brydns
c
c                    photon free energy per baryon
      pf = pu-t*ps
c
c
c                    derivatives of photon pressure
      dppdn = 0.0
      dppdt = brydns*ps
      dppdy = 0.0
c
c                    derivatives of photon entropy
      dpsdn = -ps/brydns
      dpsdt = 3.0*ps/t
      dpsdy = 0.0
c
c                    derivatives of internal energy
      dpudn = -0.75*t*ps/brydns
      dpudt = 3.0*ps
      dpudy = 0.0
c
 999  return
c
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       imtrule
c
c    type:         subroutine
c
c    author:       f. douglas swesty
c                  laboratory for computational astrophysics
c                  dept. of astronomy & ncsa
c                  university of illinois at urbana-champaign
c
c    email:        dswesty@ncsa.uiuc.edu
c
c    date:         8/25/95
c
c
c    purpose:      do the generalized fermi integrals via the imt rule
c
c
c    call line:          subroutine imtrule(eta,beta,x_lo,x_up,
c                       * n,dndb,dnde,e,dedb,dedn,p,dpdb,dpde)
c
c    inputs:       
c
c    outputs:      none
c
c
c 
c    include files:  none
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine imtrule(eta,beta,x_lo,x_up,
     *                   n,dndb,dnde,e,dedb,dede,p,dpdb,dpde)
c
      implicit none
c
      double precision eta, beta, x_lo, x_up
      double precision n, dndb, dnde
      double precision e, dedb, dede
      double precision p, dpdb, dpde
c
      integer i
      double precision fe, dfedb, dfede
      double precision xroot, xn1, dxn1db
      double precision xn2, dxn2db, xn3, dxn3db
c
      integer n_imt
      parameter(n_imt=64)
      double precision w_imt(1:n_imt-1), x_imt(1:n_imt-1)
c
c                           integration weights & limits
      double precision wght, x
c
c-----------------------------------------------------------------------
c         63 pt imt quadrature abscissae & weights
c-----------------------------------------------------------------------
c
      data x_imt(  1),w_imt(  1 ) /
     *    1.9570031734585200d-30 ,  8.2607320509999998d-27 /
      data x_imt(  2),w_imt(  2 ) /
     *    5.9130760189713709d-16 ,  6.4169110630000000d-13 /
      data x_imt(  3),w_imt(  3 ) /
     *    5.4675241449266301d-11 ,  2.7067693377000001d-08 /
      data x_imt(  4),w_imt(  4 ) /
     *    1.9306978596659198d-08 ,  5.5092726333000005d-06 /
      data x_imt(  5),w_imt(  5 ) /
     *    7.1039233799069106d-07 ,  1.3273484743000000d-04 /
      data x_imt(  6),w_imt(  6 ) /
     *    8.2971031730510406d-06 ,  1.0999107621000000d-03 /
      data x_imt(  7),w_imt(  7 ) /
     *    4.9827396171457499d-05 ,  4.9514468819000000d-03 /
      data x_imt(  8),w_imt(  8 ) /
     *    1.9629223097988500d-04 ,  1.5218120421000000d-02 /
      data x_imt(  9),w_imt(  9 ) /
     *    5.8152639478483402d-04 ,  3.6255696685999997d-02 /
      data x_imt( 10),w_imt( 10 ) /
     *    1.4073702294191501d-03 ,  7.2251913311999996d-02 /
      data x_imt( 11),w_imt( 11 ) /
     *    2.9342961722595702d-03 ,  1.2642067111999999d-01 /
      data x_imt( 12),w_imt( 12 ) /
     *    5.4624206134241907d-03 ,  2.0058566844999998d-01 /
      data x_imt( 13),w_imt( 13 ) /
     *    9.3088894298794592d-03 ,  2.9511499967999999d-01 /
      data x_imt( 14),w_imt( 14 ) /
     *    1.4786180477344599d-02 ,  4.0908191706999997d-01 /
      data x_imt( 15),w_imt( 15 ) /
     *    2.2183887527440597d-02 ,  5.4053256120000004d-01 /
      data x_imt( 16),w_imt( 16 ) /
     *    3.1754957727637700d-02 ,  6.8677770086000001d-01 /
      data x_imt( 17),w_imt( 17 ) /
     *    4.3706349739020198d-02 ,  8.4466152386999993d-01 /
      data x_imt( 18),w_imt( 18 ) /
     *    5.8193563828888392d-02 ,  1.0107865713999999d+00 /
      data x_imt( 19),w_imt( 19 ) /
     *    7.5318304334592898d-02 ,  1.1816897903000001d+00 /
      data x_imt( 20),w_imt( 20 ) /
     *    9.5128529563062905d-02 ,  1.3539730362000000d+00 /
      data x_imt( 21),w_imt( 21 ) /
     *    1.1762022895325901d-01 ,  1.5243949451000001d+00 /
      data x_imt( 22),w_imt( 22 ) /
     *    1.4274038487381499d-01 ,  1.6899319837000000d+00 /
      data x_imt( 23),w_imt( 23 ) /
     *    1.7039069590573200d-01 ,  1.8478160150000000d+00 /
      data x_imt( 24),w_imt( 24 ) /
     *    2.0043174541744399d-01 ,  1.9955546823999999d+00 /
      data x_imt( 25),w_imt( 25 ) /
     *    2.3268738852436499d-01 ,  2.1309397372999999d+00 /
      data x_imt( 26),w_imt( 26 ) /
     *    2.6694920180823400d-01 ,  2.2520473296999999d+00 /
      data x_imt( 27),w_imt( 27 ) /
     *    3.0298089532632300d-01 ,  2.3572333323999999d+00 /
      data x_imt( 28),w_imt( 28 ) /
     *    3.4052262810720502d-01 ,  2.4451259898000002d+00 /
      data x_imt( 29),w_imt( 29 ) /
     *    3.7929519918613303d-01 ,  2.5146175750999999d+00 /
      data x_imt( 30),w_imt( 30 ) /
     *    4.1900410866587701d-01 ,  2.5648562654000000d+00 /
      data x_imt( 31),w_imt( 31 ) /
     *    4.5934349927577800d-01 ,  2.5952390874000000d+00 /
      data x_imt( 32),w_imt( 32 ) /
     *    5.0000000000000000d-01 ,  2.6054065145000003d+00 /
      data x_imt( 33),w_imt( 33 ) /
     *    5.4065650072422200d-01 ,  2.5952390874000000d+00 /
      data x_imt( 34),w_imt( 34 ) /
     *    5.8099589133412299d-01 ,  2.5648562654000000d+00 /
      data x_imt( 35),w_imt( 35 ) /
     *    6.2070480081386692d-01 ,  2.5146175750999999d+00 /
      data x_imt( 36),w_imt( 36 ) /
     *    6.5947737189279498d-01 ,  2.4451259898000002d+00 /
      data x_imt( 37),w_imt( 37 ) /
     *    6.9701910467367700d-01 ,  2.3572333323999999d+00 /
      data x_imt( 38),w_imt( 38 ) /
     *    7.3305079819176600d-01 ,  2.2520473296999999d+00 /
      data x_imt( 39),w_imt( 39 ) /
     *    7.6731261147563501d-01 ,  2.1309397372999999d+00 /
      data x_imt( 40),w_imt( 40 ) /
     *    7.9956825458255598d-01 ,  1.9955546823999999d+00 /
      data x_imt( 41),w_imt( 41 ) /
     *    8.2960930409426803d-01 ,  1.8478160150000000d+00 /
      data x_imt( 42),w_imt( 42 ) /
     *    8.5725961512618498d-01 ,  1.6899319837000000d+00 /
      data x_imt( 43),w_imt( 43 ) /
     *    8.8237977104674103d-01 ,  1.5243949451000001d+00 /
      data x_imt( 44),w_imt( 44 ) /
     *    9.0487147043693705d-01 ,  1.3539730362000000d+00 /
      data x_imt( 45),w_imt( 45 ) /
     *    9.2468169566540714d-01 ,  1.1816897903000001d+00 /
      data x_imt( 46),w_imt( 46 ) /
     *    9.4180643617111159d-01 ,  1.0107865713999999d+00 /
      data x_imt( 47),w_imt( 47 ) /
     *    9.5629365026097979d-01 ,  8.4466152386999993d-01 /
      data x_imt( 48),w_imt( 48 ) /
     *    9.6824504227236230d-01 ,  6.8677770086000001d-01 /
      data x_imt( 49),w_imt( 49 ) /
     *    9.7781611247255940d-01 ,  5.4053256120000004d-01 /
      data x_imt( 50),w_imt( 50 ) /
     *    9.8521381952265541d-01 ,  4.0908191706999997d-01 /
      data x_imt( 51),w_imt( 51 ) /
     *    9.9069111057012049d-01 ,  2.9511499967999999d-01 /
      data x_imt( 52),w_imt( 52 ) /
     *    9.9453757938657583d-01 ,  2.0058566844999998d-01 /
      data x_imt( 53),w_imt( 53 ) /
     *    9.9706570382774040d-01 ,  1.2642067111999999d-01 /
      data x_imt( 54),w_imt( 54 ) /
     *    9.9859262977058083d-01 ,  7.2251913311999996d-02 /
      data x_imt( 55),w_imt( 55 ) /
     *    9.9941847360521519d-01 ,  3.6255696685999997d-02 /
      data x_imt( 56),w_imt( 56 ) /
     *    9.9980370776902017d-01 ,  1.5218120421000000d-02 /
      data x_imt( 57),w_imt( 57 ) /
     *    9.9995017260382857d-01 ,  4.9514468819000000d-03 /
      data x_imt( 58),w_imt( 58 ) /
     *    9.9999170289682693d-01 ,  1.0999107621000000d-03 /
      data x_imt( 59),w_imt( 59 ) /
     *    9.9999928960766205d-01 ,  1.3273484743000000d-04 /
      data x_imt( 60),w_imt( 60 ) /
     *    9.9999998069302143d-01 ,  5.5092726333000005d-06 /
      data x_imt( 61),w_imt( 61 ) /
     *    9.9999999994532474d-01 ,  2.7067693377000001d-08 /
      data x_imt( 62),w_imt( 62 ) /
     *    9.9999999999999944d-01 ,  6.4169110630000000d-13 /
      data x_imt( 63),w_imt( 63 ) /
     *    1.0000000000000000d+00 ,  8.2607320509999998d-27 /
c
c-----------------------------------------------------------------------
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
      do i=1,n_imt-1,1                                         !   |
c JCM:
c         write(*,*) 'loopimt',i,n_imt
        x = (x_up-x_lo)*x_imt(i)+x_lo                          !   |
        wght = w_imt(i)*(x_up-x_lo)/(dble(n_imt))        !   |
c       |                    electron occupation numbers       !   |
        fe = 1.0d0/(dexp(x+beta-eta)+1.0d0)                         !   |
c       |                    derivatives w.r.t. beta           !   |
        dfedb = -fe*(1.0d0-fe)                                 !   |
c       |                    derivatives w.r.t. eta            !   |
        dfede = fe*(1.0d0-fe)                                  !   |
c       |                    number density & derivatives      !   |
        xroot = dsqrt(x*(x+2.0d0*beta))                        !   |
        xn1 = xroot*(x+beta)                                   !   |
        dxn1db = xroot+x*(x+beta)/xroot                        !   |
c       |                    sum for electron # integral       !   |
        n = n+wght*xn1*fe                                      !   |
        dndb = dndb+wght*(dxn1db*fe+xn1*dfedb)                 !   |
        dnde = dnde+wght*xn1*dfede                             !   |
c       |                    energy integral & derivatives     !   |
        xn2 = xroot*(x+beta)**2                                !   |
        dxn2db = 2.0d0*xroot*(x+beta)+x*((x+beta)**2)/xroot    !   |
        e = e+wght*xn2*fe                                      !   |
        dedb = dedb+wght*(dxn2db*fe+xn2*dfedb)                 !   |
        dede = dede+wght*xn2*dfede                             !   |
c       |                    pressure integral & derivatives   !   |
        xn3 = xroot**3                                         !   |
        dxn3db = 3.0d0*xroot*x                                 !   |
        p = p+wght*xn3*fe                                      !   |
        dpdb = dpdb+wght*(dxn3db*fe+xn3*dfedb)                 !   |
        dpde = dpde+wght*xn3*dfede                             !   |
      enddo                                                    !   |
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
c
 999  return
c
      end
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       gl16
c
c    type:         subroutine
c
c    author:       f. douglas swesty
c                  laboratory for computational astrophysics
c                  dept. of astronomy & ncsa
c                  university of illinois at urbana-champaign
c
c    email:        dswesty@ncsa.uiuc.edu
c
c    date:         8/25/95
c
c    purpose:      do the generalized fermi integrals via 16 point
c                  gauss--laguerre quadrature
c
c
c    call line:          subroutine gl16(eta,beta,x_lo,x_up,
c                       * n,dndb,dnde,e,dedb,dedn,p,dpdb,dpde)
c
c    inputs:       
c
c    outputs:      none
c
c
c 
c    include files:  none
c
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine gl16(eta,beta,x_lo,x_up,
     *                   n,dndb,dnde,e,dedb,dede,p,dpdb,dpde)
c
      implicit none
c
      double precision eta, beta, x_lo, x_up
      double precision n, dndb, dnde
      double precision e, dedb, dede
      double precision p, dpdb, dpde
c
      integer i
      double precision fe, dfedb, dfede
      double precision xroot, xn1, dxn1db
      double precision xn2, dxn2db, xn3, dxn3db
c
      integer n_imt
      parameter(n_imt=64)
      double precision w_imt(1:n_imt-1), x_imt(1:n_imt-1)
c
c                           integration weights & limits
      double precision wght, x
c

      integer n_gl16
      parameter(n_gl16=16)
c
      double precision w_gl16(n_gl16), x_gl16(n_gl16)
c
      data x_gl16  /.087649410479d00
     *,.46269632892d00,.11410577748d01,.21292836451d01,.34370866339d01
     *,.50780186145d01,.70703385350d01,.94383143364d01,.12214223369d02
     *,.15441527369d02,.19180156857d02,.23515905694d02,.28578729743d02
     *,.34583398702d02,.41940452648d02,.51701160340d02/
      data w_gl16 /.22503631486d00
     *,.52583605276d00,.83196139169d00,.11460992410d01,.14717513170d01
     *,.18131346874d01,.21755175197d01,.25657627502d01,.29932150864d01
     *,.34712344831d01,.40200440864d01,.46725166077d01,.54874206580d01
     *,.65853612333d01,.82763579844d01,.11824277552d02/
c
c
c-----------------------------------------------------------------------
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
      do i=1,n_gl16,1                                          !   |
c JCM:
c         write(*,*) 'loopg',i,n_gl16
        x = x_gl16(i)                                          !   |
        wght = w_gl16(i)                                       !   |
c       |                    electron occupation numbers       !   |
        fe = 1.0d0/(dexp(x+beta-eta)+1.0d0)                    !   |
c       |                    derivatives w.r.t. beta           !   |
        dfedb = -fe*(1.0d0-fe)                                 !   |
c       |                    derivatives w.r.t. eta            !   |
        dfede = fe*(1.0d0-fe)                                  !   |
c       |                    number density & derivatives      !   |
        xroot = dsqrt(x*(x+2.0d0*beta))                        !   |
        xn1 = xroot*(x+beta)                                   !   |
        dxn1db = xroot+x*(x+beta)/xroot                        !   |
c       |                    sum for electron # integral       !   |
        n = n+wght*xn1*fe                                      !   |
        dndb = dndb+wght*(dxn1db*fe+xn1*dfedb)                 !   |
        dnde = dnde+wght*xn1*dfede                             !   |
c       |                    energy integral & derivatives     !   |
        xn2 = xroot*(x+beta)**2                                !   |
        dxn2db = 2.0d0*xroot*(x+beta)+x*((x+beta)**2)/xroot    !   |
        e = e+wght*xn2*fe                                      !   |
        dedb = dedb+wght*(dxn2db*fe+xn2*dfedb)                 !   |
        dede = dede+wght*xn2*dfede                             !   |
c       |                    pressure integral & derivatives   !   |
        xn3 = xroot**3                                         !   |
        dxn3db = 3.0d0*xroot*x                                 !   |
        p = p+wght*xn3*fe                                      !   |
        dpdb = dpdb+wght*(dxn3db*fe+xn3*dfedb)                 !   |
        dpde = dpde+wght*xn3*dfede                             !   |
      enddo                                                    !   |
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
c
 999  return
c
      end
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       f_2
c    type:         double precision function
c    author:       f. douglas swesty
c                  laboratory for computational astrophysics
c                  dept. of astronomy & ncsa
c                  university of illinois at urbana-champaign
c
c    email:        dswesty@ncsa.uiuc.edu
c    date:         12/16/91
c
c                  version 2: 1/24/93
c
c    call line:    f_2(y)      (2nd fermi integral)
c
c    inputs:       y (double precision)   (argument)
c
c    return:       2nd fermi integral (double precision)
c
c***********************************************************************
      double precision function f_2(y)
      implicit none
      double precision y, yexp
      if(y.gt.3.0d0) then
        write(*,*) ' f_2(y) fails for y .gt. 3; y =',y
        stop
      endif
c
c                       note: this approximation is based on the
c                       bludman & van riper approximation (see
c                       ap. j. vol. 212 page 866-867 (1977))
c                       equation (3.6)
c
      if(y.lt.-1.0d0) then
        yexp = exp(y)
        f_2 = 2.0*yexp*(1.0-0.125*yexp+0.037037*(yexp**2))
      else
        f_2 = 1.803d0+1.645d0*y+0.6931d0*(y**2)+0.1666667d0*(y**3)+
     1        2.0833333d-2*(y**4)-3.4722d-4*(y**6)
      endif
 999  return
      end   
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       matlud
c    type:         subroutine
c    author:       f. douglas swesty
c    date:         6/16/94
c
c    purpose:      lu decomposes an nxn matrix
c
c    note:         this subroutine employs crout's algorithm with
c                  implicit pivoting as described in "numerical
c                  recipes", by press, flannery, teukolsky, and
c                  vetterling, (pub. cambridge univ. press) first
c                  edition, whose authors desrve the credit for
c                  this algorithm.
c
c    call line:    call matlud(a,lu,n,ipvt)
c
c    inputs:       a = nxn array to be lu decomposed  (d)
c                  n = size of a (i)
c
c    outputs:      lu = array containing lu decomposition of a (d)
c                  ipvt = vector of pivot indices (i)
c
c    calls :       none
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine matlud(a,lu,n,ipvt)
c
      implicit none
c
c                 parameters
      integer n
      integer ipvt(n)
      double precision a(n,n), lu(n,n)
c
c                 local variables
c
c                 loop variables
      integer i, j, l
c                 a small floating point value to prevent division
c                 by zero
      double precision smallv
      parameter(smallv=1.0d-20)
c
c                 value & original row index or largest l value
      double precision e_max, tstval
      integer rowpvt
c
c                 a scratch variable to use when swapping rows
      double precision scrvar
c
c                 an array for the largest value of each row of a
      integer nmax
      parameter(nmax=100)
      double precision maxval(nmax)
c-----------------------------------------------------------------------
c
c                 find the maximum absolute value in each row of a
      do 20 i=1,n,1
        maxval(i) = -1.0d0
        do 10 j=1,n,1
c JCM:
c           write(*,*) 'loopmaxval',j,i,n
          maxval(i)  = dmax1(smallv,dabs(a(i,j)),maxval(i))
          lu(i,j) = a(i,j) 
10     continue
 20   continue
c
c                 now employ crout's algorithm with implicit pivoting
      do 90 j=1,n,1
c
c                 calculate column j, u matrix elements
        do 40 i=1,j-1,1
          do 30 l=1,i-1,1
c JCM:
c             write(*,*) 'loopcol',j,i,l,n
            lu(i,j) = lu(i,j)-lu(i,l)*lu(l,j)
 30       continue
 40     continue
c
c                 calculate column j, l matrix elements.  the element is
c                 scaled by the largest element of the original matrix.
c                 also, the column from the diagonal down is searched
c                 for the largest element.
        e_max = -1.0d0
        do 60 i=j,n,1
          do 50 l=1,j-1,1
c JCM:
c             write(*,*) 'looplu',i,l
            lu(i,j) = lu(i,j)-lu(i,l)*lu(l,j)
 50       continue
          tstval = dabs( lu(i,j) )/maxval(i)
          if(tstval.gt.e_max) then
            e_max=tstval
            rowpvt = i
          endif
 60     continue
c                    keep track of which row was pivoted into row j
       ipvt(j) = rowpvt
c
c
c
c                 if the original diagonal element wasn't the
c                 largest, then swap row j with row rowpvt
        if(rowpvt.ne.j) then
          do 70 l=1,n,1 
c JCM:
c            write(*,*) 'loopdiag',l
             scrvar = lu(rowpvt,l)
             lu(rowpvt,l) = lu(j,l)
             lu(j,l) = scrvar
 70       continue
c                    set the maxval for row rowpvt to that of the one
c                    that was just swapped in. 
          maxval(rowpvt) = maxval(j)
c
        endif
c
c                  now divide the rest of the l column by the pivot
c                  element
        if(dabs(lu(j,j)).gt.smallv )  then
          scrvar = lu(j,j)
        else
          scrvar = sign(smallv,lu(j,j))
        endif
        do 80 i=j+1,n,1
c JCM:
c           write(*,*) 'loopsign',i
          lu(i,j) = lu(i,j)/scrvar
 80     continue
c
c
 90   continue
c
c
 999  return
c
      end
c
c***********************************************************************
c
c    module:       mluslv
c    type:         subroutine
c    author:       f. douglas swesty
c    date:         6/16/94
c
c    purpose:      lu decomposes an nxn matrix
c
c    note:         this subroutine employs a forward & back substitution
c                  algorithm with implicit pivoting as described in 
c                  "numerical recipes", by press, flannery, teukolsky, &
c                  vetterling, (pub. cambridge univ. press) first
c                  edition, whose authors desrve the credit for
c                  this algorithm.
c
c    call line:    call mluslv(lu,x,b,ipvt,n)
c
c    inputs:       a = nxn lu decomposed array (d)
c                  b  = rhs vector (d)
c                  ipvt = vector of pivots (i)
c                  n = size of a (i)
c
c    outputs:      x = solution vector(d)
c
c    calls :       none
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine mluslv(lu,x,b,ipvt,n)
c
      implicit none
c
c                 parameters
      integer n
      integer ipvt(n)
      double precision lu(n,n), x(n), b(n)
c
c                 local variables
c
c                 loop variables
      integer i, j, l
c
c                 a scratch variable
      double precision scrv
c
c-----------------------------------------------------------------------
c
c                 copy the rhs into x so we don't destroy b by
c                 unscrambling the pivots
      do 10 i=1,n,1
c JCM:
c         write(*,*) 'unsc',i
        x(i) = b(i)
 10   continue
c
c                 do the forward substitution
      do 30 i=1,n,1
        l = ipvt(i)
        scrv = x(l)
        x(l) = x(i)
        x(i) = scrv
        do 20 j=1,i-1,1
c JCM:
c         write(*,*) 'for',i,j
          x(i) = x(i)-lu(i,j)*x(j)
 20     continue
 30   continue
c
c
c                 do the backward substitution
      do 50 i=n,1,-1
        do 40 j=i+1,n,1
c JCM:
c         write(*,*) 'back',i,j
          x(i) = x(i)-lu(i,j)*x(j)
 40     continue
        x(i) = x(i)/lu(i,i)
 50   continue
c
c
 999  return
c
      end


c JCM: presumed already included elsewhere
c      include 'jon_eosfxt.f'
