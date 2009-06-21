c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       maxwel.inc
c    type:         include file
c    author:       f. douglas swesty
c    date:         6/13/90
c
c
c    call line:    include 'maxwel.inc/list'
c
c*************************************************************************
c
c                            number of t points in each boundary
      integer numtmp
      parameter (numtmp=201)
cc      parameter (numtmp=21)
c
c                            number of ye points in each boundary
      integer numye
cc      parameter(numye=49)
cc      parameter(numye=34)
      parameter(numye=49)
c
c                            number of points in alpha-nuclei
c                            boundary
      integer nbpnts, numlow, numhi
      parameter(nbpnts=101,numlow=51,numhi=51)
c
c
c                            baryon densities
      double precision brylow(numtmp,numye), bryhi(numtmp,numye)
      double precision lowdns, hidns, dns_1, dns_2
      double precision dnl_dt, dnh_dt, dnl_dy, dnh_dy
c
c
      double precision lnlow, lnhi, lncut, lnmins, lnplus
      double precision logbry, logbch, dltln1, dltln2, lnfrac
      parameter(lnlow=-8.92d0,lnhi=-0.92d0,lncut=-2.92d0)
c
c
      double precision ylow,yhi,y_cut
cccc 180 mev parms      parameter(ylow=0.03,yhi=0.51,y_cut=0.175)
cccc 220 mev parms      parameter(ylow=0.03,yhi=0.51,y_cut=0.135)
cccc 375 mev parms      parameter(ylow=0.03,yhi=0.51,y_cut=0.155)
cccc pearson et al. parms parameter(ylow=0.03,yhi=0.51,y_cut=0.125)
      parameter(ylow=0.03d0,yhi=0.51d0,y_cut=0.155d0)
c
c
      double precision tchk_b, tchk_n, t_mxwl, d_mxwl
c
c
c                            highest temperatures where coexistence
c                            of bulk and nuclear phases occurs
      double precision t_h(numye), d_h(numye)
c
c
c                            arrays containing number of boundary
c                            points
      double precision lbound(nbpnts,numye), ubound(nbpnts,numye)
c
c                            minimum density code will work at
c
      double precision mindns
      parameter(mindns=1.0d-10)
c
c
c
      double precision newdns
      double precision t_low, t_hi, tfrac
      double precision y_low, y_low2, y_hi, y_hi2, yfrac, ych, ycut
      double precision lnl, lnh, lnc
c
      double precision yminus, yplus, yintrp, tminus, tplus, tintrp
      double precision delt_y, delt_t
      integer i_mxwl, j_mxwl
      integer i_bd, i_bndy, j_bd, j_bndy
c
c                            total chemical potential
      double precision mutild, mutlow, muthi
c
c                            total pressure
      double precision prtild, prlow, prhi
c
c
c                            muhat
      double precision muhlow, muhhi
c
c                            electron chemical potential
      double precision muelow, muehi
c
c                            total entropy per baryon
      double precision s_low, s_hi
c
c                            total free energy density
      double precision f_low, f_hi
c
c
c                            phase fraction
      double precision phasef
c
c
      common /maxcom/ t_low,t_hi,y_low,y_hi,
     1       brylow,bryhi,lbound,ubound,t_h,d_h,
     2       lowdns,hidns,tchk_n,tchk_b,t_mxwl,d_mxwl,
     3       dnl_dt,dnh_dt,dnl_dy,dnh_dy
