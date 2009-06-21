
c
cJCM:
c All these things used in LSEOS.  There are several parameters.
c Near the bottom are the global variables
c These globals are used to return result of LSEOS
c None of these globals should be same name as other globals that want to be
c  written to file
c
c
cc JCM: Some convergence parameters for LSEOS code itself
c..
c JCM:
      double precision lstol1
c   original is 1E-5
      parameter (lstol1=1E-8)
c      parameter (lstol1=1E-5)




c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       eos_m4c.inc      (version 2.7)
c    type:         include file
c    author:       f. douglas swesty
c    date:         (original) 6/27/90
c                  (v2.7)     9/15/95
c
c
c    call line:    include 'eos_m4c.inc'
c
c*************************************************************************
c
c                   input & output varibles
c
c                            electron fraction
      double precision ye
c                            baryon density
      double precision brydns
c                            previous value of the nuclear proton
c                            fraction (must be supplied)
      double precision xprev
c
c                            previous value of the proton number
c                            density (must be supplied)
      double precision p_prev
c
c                            temperature, density, neutron chem. pot.,
c                            and proton chem. pot.
      double precision inpvar(4)
c
c                            flag to tell code what input is being
c                            supplied (i.e. what inpvar is)
c                            1 --> temperature
c                            2 --> internal energy (not implemented)
c                            3 --> entropy (not implemented)
      integer iflag
c
c                            flag returned by code to tell user what
c                            scheme was used
c                            1 --> "no nuclei" scheme
c                            2 --> full scheme
c                            3 --> bulk scheme (above saturation dnsty)
      integer eosflg
c
      integer forflg, ngflag, newflg, ssflag, fflag, adflag
      integer nf_flg, rsflag, dbflag
cc      data nf_flg/1/
cc      data rsflag/0/
c
c*************************************************************************
c                        local variables
c*************************************************************************
c
c                            flag used internally to this routine to
c                            record whether or not another scheme
c                            has already been tried
      integer swtflg
c
c
c                            accuracy in x
      double precision xacc
c     JCM:
      parameter (xacc=1.0d-10)
c      parameter (xacc=1.0d-8)
c
c                            maximum number of x iterations
      integer maxit, maxit2
c
cJCM:
      parameter (maxit=500, maxit2=200)
c      parameter (maxit=200, maxit2=60)
c
c                            update flag
      integer upflag
      parameter(upflag=3)
c
c                            external interaction flag
c                            1 = on, 0 = off
      double precision eiflag
      parameter(eiflag=1.0)
c
c                            alpha particle enable flag
c                            1= alphas, 0 = no alphas
      integer alflag
      parameter(alflag=1)
c
c                            translational energy scaling factor
c                            1 = on, 0 = off
      double precision trscal
      parameter(trscal=1.0)
c
c                            surface and coulomb scaling factor
c                            1.0 = on, 0 = off
      double precision csscal
      parameter(csscal=1.0)
c
c
c
c                            cutoff fraction of nuclei (if xh is below
c                            this the "no nuclei" scheme is used)
      double precision heavct
cJCM:
      parameter(heavct=1.0d-8)
c      parameter(heavct=1.0d-6)
c
c                            iteration loop variable
c      integer i, j, l, it_num
      integer i, j, l
c
c                            combination of coulomb and surface
c                            coefficients
      double precision alpha, beta
      parameter (alpha=224.4, beta=1.3927)
c
c                            surface tension (mev/fm**2)
      double precision sig_0
      parameter(sig_0=46.5)
c                            binding energy of alpha particles (mev)
      double precision balpha
      parameter (balpha=28.3)
      double precision v_alfa
      parameter (v_alfa=24.0)
c
c                            nuclear symmetry energy (mev)
      double precision sv
      parameter (sv=29.3)
c
c                            nuclear compression modulus (mev)
      double precision comprs, k
      parameter (comprs=370.0)
c
c                            nuclear level density (per mev)
      double precision levdns, av
      parameter (levdns=0.0666667)
c
c                            approximate size of nuclei
      double precision azero
      parameter (azero=60.0)
c
c
c                           plancks constant & speed of light
      double precision hbar, c
c     JCM: hbar in "MeV s" and c in "fm/s"
      parameter (hbar=6.58217317d-22,c=2.997924581d23)
c
c                           pi and square root of pi
c      double precision pi, rootpi
c JCM:
c      parameter(pi=3.1415927,rootpi=1.7724539)
c      parameter(pi=3.1415926535897932384d0,rootpi=1.7724538509055160273d0)
c
c                           square root of two
      double precision root2
c JCM:
      parameter (root2=1.4142135623730950488)
c      parameter (root2=1.4142136)
c
c                           parameters containing powers of 1/3
c
      double precision ovr3, ovr23, ovr43, ovr53, ovr29
      double precision movr3, m2ovr3, m4ovr3, m5ovr3
      double precision ovr49
      parameter(ovr49=0.444444444444444)
      parameter (ovr3=0.3333333333333, ovr23=0.6666666666666, ovr43=1.333333333333)
      parameter (ovr53=1.666666666666, ovr29=0.222222222222222)
      parameter (movr3=-0.3333333333333333, m2ovr3=-0.6666666666666666)
      parameter (m4ovr3=-1.333333333333333, m5ovr3=-1.66666666666666666)
c
c                          ratio of baryon density to saturation density
      double precision y
c
      double precision tfinal
c
c                          quantum concentration & fermi integral
c                          coefficent
      double precision nq, mq, lq, kq, uq
c
c                          estimate of x (used in "no nuclei" scheme)
      double precision xguess
c
c                          cutoff used to determine how to calc the
c                          number density of heavy nuclei (used in
c                          the "no nuclei" scheme)
      double precision mhcut, mhovt        
      parameter(mhcut=30.0)
c
c                          chem. pot., mass, # density, & quantum
c                          concentration of heavy nuclei (used in
c                          the "no nuclei" scheme)
      double precision muheav, nucq, masshv, nheavy
c
c
c
c
c                          bulk nuclear pressures
      double precision psubi, povrni
c
c                          neutron & proton rest masses
      double precision massn, massp
      parameter(massn=939.5731,massp=938.2796)
c
c                          saturation density, nuclear density, the
c                          difference of the two and their ratio
c      double precision nsubi, n_i, dni, nratio, n_iold
      double precision n_i, nratio, n_iold
      double precision niold
c
c
c                          surface energy coefficants
      double precision n_s
      double precision q, r_0, zeta_0
c
c
c
c                          coulomb liquid correction functions
c                          & derivatives
      double precision w, dwdx, dwdt, dwdtdx, dwdxdx, dwdtdt
      double precision tzero, dtzdx, dtzdxx
c
c
c                          coulomb functions and derivatives
      double precision du, dmu, dup, dmup, dupp, dmupp
      double precision dux, dmux, dupx, dmupx
      double precision dut, dmut, dupt, dmupt, duppt
      double precision duxt, dmuxt
      double precision duxx, dmuxx, dutt, dmutt
c
      double precision scrdu, scrdup, scrdpp, scrd, scrd2
      double precision scrdut, scrdux, scrdxx, scrdt, scrdx, scrdxt
      double precision scrdtt
      double precision scrd2t, scrd2x, scrdpt, scrdpx
      double precision utilde, dutil, dutilp, omeg_r
c
c
c                          finite size energy coefficent
      double precision zeta, dzdt, dzdx, dzdni
c
c                          equilibrium equation finite size terms
      double precision bn, bp, bsubp
c
c                          baryon densities
      double precision nnout, alfdns, nucdns, nout
c
c                          baryon fractions
c      double precision xalfa, xnut, xprot, xh, xhchk
      double precision xhchk
c      double precision xalfa2, xnut2, xprot2
c
c                          fraction of space occupied by nuclei
c      double precision u_nuc, u_n, compu, ry
      double precision u_n, compu, ry
c                          fraction of space available to remaining
c                          baryons  (1-u_nuc)
      double precision exalfa, exclu
c
c
c                          surface tension variables
      double precision compx, sigma, ahbn, ovrx4, dsigdx, dsgrdx, sigsgp
      double precision sigsg2
c
c
c                          critical temperature variables
      double precision tsubc, tsc1, dtcdx, dtcdxx, tsc_12
c
c                          surface & translational temperature variables
      double precision dhdx, dhdxx, hx, dhdtdx, chova
      double precision h, hprim, hpprim
      double precision caph, caphp
c
c                          neutron-proton mass difference
      double precision deltam
      parameter(deltam=-1.2935)
ccccccc      parameter(deltam=0.0)
c
c                          chemical potentials
      double precision munovt
      double precision mup_i, mun_o, mup_o, mun_z, mup_z
      double precision eta_no, eta_po, eta_ni, eta_pi, etamax
      double precision etp_hi, etp_lo, etn_hi, etn_lo
      double precision deltmn, deltmp, deltap
c
      double precision dmndpi, dmpdpi, dpidni
c
      double precision omega
c
      double precision deta_p, deta_n, neta_n, neta_p
      double precision chng_n, chng_p
c
c
c
c                          n-n interaction variables
      double precision vnout, vpout, vni, vpi
      double precision zno, zpo, zni, zpi
c
c                          nucleon kinetic energy densities
      double precision tau_po, tau_no, tau_ni, tau_pi
      double precision f32_po, f32_no, f32_ni, f32_pi
      double precision frat_n, frat_p
c
c                          newton-raphson equation variables
      double precision g, dgdx, dgdprt, gold
c
c                          atomic weight & number
c      double precision a, z
      double precision z
c
c                          nuclear radius & volume
c      double precision rsubn, vsubn
c
c                          surface, coulomb, translational, and 
c                          bulk free energies (per baryon)
      double precision fsubs, fsubc, fsubsc, fsubi, ftrans, musubt
      double precision dmutdt
c
      double precision ftr, dftrdt, dftrdu, dftrdx, dftrdn
      double precision f_sc, dfscdt, dfscdu, dfscdx, dfscdn
c
      double precision e_tr, detrdt, detrdu, detrdx, detrdn
      double precision e_sc, descdt, descdu, descdx, descdn
c
      double precision s_tr, dstrdt, dstrdu, dstrdx, dstrdn
      double precision s_sc, dsscdt, dsscdu, dsscdx, dsscdn
c
      double precision s_nuc, s_out, s_alfa
      double precision e_nuc, e_out, e_alfa
c
      double precision s_dens, e_dens
c
      double precision detpdx, detndx, detpdn, detndn
      double precision dtpidx, dtnidx, dtpdni, dtndni
c
      double precision dsidt, dsidx, dsidn, deidt, deidx, deidn
c
      double precision dsodt, dsodep, dsoden, deodt, deodep, deoden
c
      double precision dsadt, dsadep, dsaden, deadt, deadep, deaden
c
      double precision dnpdep, dnpden, dnndep, dnnden
      double precision dtpdep, dtpden, dtndep, dtnden
c
c
c
c
c
c                            baryon thermodynamic variables
c
c                          pressures
c      double precision bprout, bprnuc, bpralf
c
c                          entropies (per baryon)
c      double precision bsout, bsnuc, bsalfa
      double precision bsalfa
      double precision dsfsdu, dsfsdx, dsfsdn, dsfsdt
c
c
c
      double precision dbsodt, dbsndt, dbsadt, sfs_dt
c
c                          internal energies (per baryon)
c      double precision buout, bunuc, bualfa
      double precision bualfa
c
c                          helmholtz free energies (per baryon)
c      double precision bfout, bfalfa, bfnuc, bftot

c                          adiabatic index
c      double precision gam_s
c
c                            total energy, entropy, and pressure
c
c      double precision utot, stot, ptot, ftot
c      double precision fbary, pbary, mubary, mu_mat
c
c                            x variables
c
c                          internal proton fraction (in nuclei)
c      double precision x
c                          outside proton fraction (nucleon vapor)
      double precision xout
c
      double precision xnew, xold, dx
c                          limits on x
      double precision x_min
cc      double precision xl, xhi
cc      parameter(xl=0.01,xhi=0.52)
c
c
c                            proton variables
      double precision npout, prtnew, prtold, dprt
c
c
c                            temperature variables
c
      double precision t
c
c                            diagonal line in the temp-density plane
      double precision tcheck
c
c                            temporary storage variables
      double precision tmp1, tmp1p, tmp1pp
      double precision tmp2, tmp2p, tmp2pp, tmp2pt, tmp2px
      double precision tmp2t, tmp2x, tmp2xx, tmp2xt, tmp2tt
      double precision tmp3, tmp3p, tmp3pp, tmp3pt, tmp3px
      double precision tmp3t, tmp3x, tmp3xx, tmp3xt, tmp3tt
      double precision tmp4, tmp5
c
c
c
      double precision nsubin, netap, netan
c
      double precision nsiacc, prtacc, nutacc
c..      parameter(nsiacc=1.0d-8,prtacc=1.0d-9,nutacc=1.0d-9)
      parameter(nsiacc=1.0d-10,prtacc=1.0d-11,nutacc=1.0d-11)

c
      double precision determ
c
      double precision dnsubi, detap, detan
c
      double precision a1, a2, a3
c
      double precision a1h, a1l, da1dn, da1etp, da1etn
      double precision a2h, a2l, da2dn, da2etp, da2etn
      double precision a3h, a3l, da3dn, da3etp, da3etn
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      double precision gpi,gni,gpo,gno
      double precision dvpidp, dvpidn, dvnidp, dvnidn
      double precision dvpodp, dvpodn, dvnodp, dvnodn
      double precision msscon
      double precision dtpidp, dtpidn, dtnidp, dtnidn
      double precision dtpodp, dtpodn, dtnodp, dtnodn
      double precision dmpidp, dmpidn, dmnidp, dmnidn
      double precision dmpodp, dmpodn, dmnodp, dmnodn
      double precision dpidp, dpidn, dpodp, dpodn, dpadp, dpadn
      double precision n1, n2, dudpo, dudno, dudni
      double precision dxdpo, dxdno, dxdni
      double precision db1dni, db1dx, db1du
      double precision db2dni, db2dx, db2du
      double precision db3dni, db3dx, db3du
      double precision da1id1, da1id2, da1id3, da1od1, da1od2, da1od3
      double precision da2id1, da2id2, da2id3, da2od1, da2od2, da2od3
      double precision da3id1, da3id2, da3id3, da3od1, da3od2, da3od3
      double precision da1d1, da1d2, da1d3, db1d1, db1d2, db1d3
      double precision da2d1, da2d2, da2d3, db2d1, db2d2, db2d3
      double precision da3d1, da3d2, da3d3, db3d1, db3d2, db3d3
      double precision dndetn, dpdetp
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c
c/////////////////////////////////////////////////////////////////////
c                   derivative variables
c/////////////////////////////////////////////////////////////////////
c
c
c                   temperature derivatives
c
      double precision dmpodt, dmnodt, dmpidt, dmnidt
      double precision dnadpo, dnadno
c
      double precision dtpidt, dtnidt, dtpodt, dtnodt
      double precision dpodt, dpidt
c
      double precision dnpodt, dnnodt
c
      double precision dv_dt
c
      double precision dnadt, dpadt, dmuadt
c
c
c      double precision du_dt, dx_dt, dni_dt, dep_dt, den_dt
c
c      double precision du_dn, dx_dn, dni_dn, dep_dn, den_dn
c
c      double precision du_dy, dx_dy, dni_dy, dep_dy, den_dy
c
      double precision dna_dt, dna_dy, dna_dn
c
      double precision db1dt, db2dt, db3dt
c
c
c
c                   exterior particle density derivatives
c
c
      double precision dv_dpo, dv_dno
c
c
c
c
c                   eta derivatives
c
c
      double precision dmpdep, dmpden, dmndep, dmnden
c
      double precision dv_dep, dv_den
c
      double precision dnadep, dnaden, dpadep, dpaden
c
      double precision dpodep, dpoden
c
      double precision dnpodn, dnpidn, dnnidn
c
      double precision dmadep, dmaden
c
c
c
      integer ipvt(5)
      double precision dfdom(5,5), dfdomi(5,5), dfdmlu(5,5)
      double precision dfdl_1(5), dfdl_2(5), dfdl_3(5)
c
      double precision dg1do1, dg1do2, dg2do1, dg2do2, det_gt
      double precision dg1dl1, dg1dl2, dg1dl3
      double precision dg2dl1, dg2dl2, dg2dl3
c
      double precision a_1x, a_1u, a_1n
      double precision a_2x, a_2u, a_2n
      double precision a_3x, a_3u, a_3n
      double precision det_al
c
      double precision at_11, at_12, at_13
      double precision at_21, at_22, at_23
      double precision at_31, at_32, at_33
c
      double precision b_1p, b_2p, b_3p
      double precision b_1n, b_2n, b_3n
c
c
      double precision dnidpo, dnidno
c
c
      double precision c_p, c_n, cp_p, cp_n
      double precision bb_p, bb_n, bp_p, bp_n
c
      double precision dmp_dn, dmn_dn
c
      double precision dmp_dy, dmn_dy
c
      double precision dfdndn, dfdydy, dfdndy, dfdndt, dfdydt
c
      double precision dfdtdt
c
c      double precision dbpdt, dbpdn, dbpdy
c
c      double precision dbsdt, dbsdn, dbsdy
c
c      double precision dbudt, dbudn, dbudy
c
c      double precision dbmudt, dbmudn, dbmudy
c
c      double precision dmudt, dmudn, dmudy, dsdt, dsdn, dsdy
c      double precision dpdt, dpdn, dpdy, dudt, dudn, dudy
c
c
      double precision npi, nni, dnpidt, dnnidt
c
      double precision dxhdt, dxodt, dxadt
c
      double precision vp_in, vn_in
c
c
      double precision duodt, duidt, duadt
c
      double precision dbfdt, dbfdy, dbfdn


c
c
c
c
c
c
c/////////////////////////////////////////////////////////////////////
c
c
c
c
c                   this common block contains all of the thermodynamic
c                   and compositional variables that the eos code
c                   calculates & it should be included anywhere
c                   these variables need to be referenced
c
c
c      common /eosvar/ mun,mun_i,muprot,mualfa,muhat,
c     1    bu, bs, bpress, utot, stot, ptot,
c     2    xprot, xnut, xalfa, xh, x, a,
c     3    buout, bunuc, bprout, bprnuc, bpralf, bsout,bsnuc,bfout,
c     4    bfnuc, bfalfa, ftot, bftot,u_nuc, nsubi, dni, rsubn,
c     5    vsubn, fbary, pbary, mubary, mu_mat,
c     6    xalfa2, xnut2, xprot2, it_num
c
c
c                    this common contains the variables that are used 
c                    in solving the chemical and pressure equilibrium 
c                    equations
c
      common /eqleqn/ nnout, npout, nout, vnout, vpout, f32_no, f32_po,
     1 mun_o, mup_o, alfdns, compx, sigma, sigsgp, du, dmu, dup, dmup,
     2 scrdu, scrdup, scrd, zeta, fsubs, fsubc, musubt, ftrans, nratio,
     3 vni, vpi, eta_ni, eta_pi, mup_i, f32_ni, f32_pi, bn, bp, bsubp,
     4 h, hprim, psubi, tsubc, dhdx, dhdxx, scrd2, scrdpp, tfinal
c
c
c
c
c                   this common block contains control flags
      common /flags/ adflag, nf_flg, rsflag, dbflag
c
c
c                   this common block contains the quality of solution
c                   variables that determine how "zeroed" the equilbrium
c                   equations are
      common /qcheck/ a1,a2,a3
c
      double precision ncomp
      common /dtest/ dfdom, dfdl_1,dfdl_2,dfdl_3,dna_dt,dna_dn,dna_dy,
     1 ncomp
c





c..this controls the addition of baryons, e+e-, and photons to the totals
      double precision radadd,eleadd,ionadd
      parameter        (radadd = 1.0d0,
     1                  eleadd = 1.0d0,
     2                  ionadd = 1.0d0)


      include 'eos_m4c.commononly.inc'










