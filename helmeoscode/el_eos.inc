c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c    module:       el_eos.inc
c    type:         include file
c    author:       f. douglas swesty
c    date:         2/12/91
c
c
c    call line:    include 'el_eos.inc'
c
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c
c                              electron density and concentration
      double precision nsube, qsube
c
c                              positron density (not set in our
c                              simple ultrarelativistic eos)
      double precision neplus
c
c                              temporary coefficients
      double precision acoef, bcoef, ccoef, dbdt
c
c                              electron chemical potential
      double precision musube
c
c                              electron pressure, internal energy,
c                              helmholtz free energy, and entropy
      double precision epress, eu, fsube, es
c
c                              photon pressure, internal energy,
c                              helmholtz free energy, and entropy
      double precision ppress, ps, pu, pf
c
c
      double precision demudt, demudn, demudy
      double precision depdt, depdn, depdy
      double precision desdt, desdn, desdy
      double precision deudt, deudn, deudy
c
c
      double precision dppdn, dppdt, dppdy
      double precision dpsdn, dpsdt, dpsdy
      double precision dpudn, dpudt, dpudy
c
c
c                              this common block contains all the
c                              electron and phton variables needed
c                              by the outside world
      common /electr/ nsube, neplus, musube, 
     1  epress, eu, fsube, es, ppress, ps, pu, pf,
     2  demudt, demudn, demudy, depdt, depdn, depdy,
     3  desdt, desdn, desdy, deudt, deudn, deudy,
     4  dppdn, dppdt, dppdy, dpsdn, dpsdt, dpsdy,
     5  dpudn, dpudt, dpudy
c

c
c



