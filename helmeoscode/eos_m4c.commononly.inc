c Only those things shared out of nuclear EOS


      double precision mun,mun_i,muprot,mualfa,muhat,etapls,etanls,
     1    bu, bs, bpress, utot, stot, ptot,
     2    xprot, xnut, xalfa, xh, x, a,
     3    buout, bunuc, bprout, bprnuc, bpralf, bsout,bsnuc,bfout,
     4    bfnuc, bfalfa, ftot, bftot,u_nuc, nsubi, dni, rsubn,
     5    vsubn, fbary, pbary, mubary, mu_mat,
     6    xalfa2, xnut2, xprot2, it_num


      common /eosvar/ mun,mun_i,muprot,mualfa,muhat,etapls,etanls,
     1    bu, bs, bpress, utot, stot, ptot,
     2    xprot, xnut, xalfa, xh, x, a,
     3    buout, bunuc, bprout, bprnuc, bpralf, bsout,bsnuc,bfout,
     4    bfnuc, bfalfa, ftot, bftot,u_nuc, nsubi, dni, rsubn,
     5    vsubn, fbary, pbary, mubary, mu_mat,
     6    xalfa2, xnut2, xprot2, it_num



c                   this common block contains the derivatives of
c                   the thermodynamic potentials (both baryon and
c                   the totals)
c
      double precision dbmudt, dbmudn, dbmudy, dmudt, dmudn, dmudy,
     1                dbpdt, dbpdn, dbpdy, dpdt, dpdn, dpdy,
     2                dbsdt, dbsdn, dbsdy, dsdt, dsdn, dsdy,
     3                dbudt, dbudn, dbudy, dudt, dudn, dudy, gam_s,
     4                du_dt, dx_dt, dni_dt, dep_dt, den_dt,
     5                du_dy, dx_dy, dni_dy, dep_dy, den_dy,
     6                du_dn, dx_dn, dni_dn, dep_dn, den_dn

      common /derivs/ dbmudt, dbmudn, dbmudy, dmudt, dmudn, dmudy,
     1                dbpdt, dbpdn, dbpdy, dpdt, dpdn, dpdy,
     2                dbsdt, dbsdn, dbsdy, dsdt, dsdn, dsdy,
     3                dbudt, dbudn, dbudy, dudt, dudn, dudy, gam_s,
     4                du_dt, dx_dt, dni_dt, dep_dt, den_dt,
     5                du_dy, dx_dy, dni_dy, dep_dy, den_dy,
     6                du_dn, dx_dn, dni_dn, dep_dn, den_dn

