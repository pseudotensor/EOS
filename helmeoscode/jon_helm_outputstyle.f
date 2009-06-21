

      subroutine pretty_eos_out(whose)
      implicit none
      save

      include 'eosparms.f'
      include 'vector_eos.dek'
      include 'vector_eos.single.dek'
      include 'extra_vector_sneos.dek'
      include 'kazeos.dek'
c      include 'kaz_state.dek'

c..bring in the lattimer-swesty commons
      include 'eos_m4c.commononly.inc'
      include 'vector_sneos.dek'


c..   
c..   writes a pretty output for the eos tester
c..   
c..   declare
      integer      jj
      character*7  whose


      real*8 utot_helm, uion_helm,uele_helm,urad_helm,ucou_helm

      include 'const.dek'







c..   popular formats
 01   format(1x,t2,a,t11,'total',t24
     1     ,'ion',t34,'ion+cou',t46,'e- & e+'
     1     ,t58,'radiation',t71,'coulomb')
 02   format(1x,t2,a,1p6e12.5)
 03   format(1x,t2,a6,1pe12.5,t22,a6,1pe12.5,
     1     t42,a6,1pe12.5,t62,a6,1pe12.5)


      open(18,file='eoscoulomb.dat',status='old',access='append')
      open(19,file='eosazbar.dat',status='old',access='append')
      open(20,file='eosdetails.dat',status='old',access='append')




ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     loop over vector of states
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do jj=jlo_eos,jhi_eos




c     Modify the variables so in consistent physical format and Kaz-like variables/physics format
c     Converts _row (both  HELM/TIMMES format and new pure Kaz _row's) into single non-row quantities
         call compute_kazlike_eos(jj)









ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        NOW OUTPUT
c
c No more physics, just dumping
c


cccccccccccccccccc
c
c Output Kaz-like variables (including primary data used by to interpolate to use in HARM)
c
c
         call kaz_output()




cccccccccccccccccc
c
c Rest of quantities are only generated with HELM/TIMMES/LSEOS
c
c
c         if( (whicheleeos.eq.0).OR.(whicheleeos.eq.1).OR.(whicheleeos.eq.2)
c     1 ) then
            

c     GODMARK: Coulomb terms depend upon abar and zbar
c     Output ion as separate since "N" above includes ion+Coulomb
            write(18,100) pcou_row(jj), ecou_row(jj)*rhob, scou_row(jj)*rhob,
     |           sion_row(jj)*rhob, plasg_row(jj) !eoscoulomb.dat

c     output abar and zbar and other things of interest
c     abar, zbar, n_i, n_{e in ion}, cv, cp, gam1, gam2, gam3, c_s
            write(19,100) abar_row(jj), zbar_row(jj)
     |           ,xni_row(jj), xnem_row(jj),xne_row(jj)
     |           ,xneut_row(jj), xprot_row(jj), xalfa_row(jj),xheav_row(jj),aheav_row(jj),zheav_row(jj),xcheck_row(jj),muhat_row(jj)
     |           ,cv_row(jj), cp_row(jj), gam1_row(jj), gam2_row(jj), gam3_row(jj)
     |           ,cs_row(jj),didconverge_row(jj)
     |           ,dse_row(jj),dpe_row(jj),dsp_row(jj)
     |           ,dse_ls_row(jj),dpe_ls_row(jj),dsp_ls_row(jj)

c         end if


ccccccccccccccccccccccccccccccccccccccccccccccc











ccccccccccccccc
c
c Below quantities are outputted only if HELM/TIMMES/LSEOS
c
         if((OUTPUTDETAILS.eq.1)
     1        .AND.(
     1        (whicheleeos.eq.0).OR.(whicheleeos.eq.1)
     1        .OR.(whicheleeos.eq.2))
     1   ) then

         
c..   the input 
            write(20,*) 'jj=',jj
            write(20,03) 'temp =',temp_row(jj),'den  =',den_row(jj),
     1           'abar =',abar_row(jj),'zbar =',zbar_row(jj)
            write(20,*) ' ' 

            utot_helm=etot_row(jj)*den_row(jj)
            uion_helm=eion_row(jj)*den_row(jj)
            uele_helm=eele_row(jj)*den_row(jj)
            urad_helm=erad_row(jj)*den_row(jj)
            ucou_helm=ecou_row(jj)*den_row(jj)

c..   and the output
c..   first the totals from each of the components
c            write(20,01)  whose
c            write(20,02) 'Tpres=',
c     1           p_tot,p_N,p_eleposi,
c     2           p_photon,0.0
c            write(20,02) 'u    =',
c     1           u_tot,u_N,rho_eleposi,
c     2           rho_photon,0.0
c            write(20,02) 'entr =',
c     1           s_tot,s_N,s_eleposi,
c     2           s_photon,0.0



c..   and the output
c..   first the totals from each of the components
            write(20,*)  ' '
            write(20,01)  whose
            write(20,02) 'pres =',
     1           ptot_row(jj),pion_row(jj),pion_row(jj)+pcou_row(jj),pele_row(jj)+ppos_row(jj),
     2           prad_row(jj),pcou_row(jj)
            write(20,02) 'u    =',
     1           utot_helm,uion_helm,uion_helm+ucou_helm,uele_helm,
     2           urad_helm,ucou_helm
            write(20,02) 'sner =',
     1           etot_row(jj),eion_row(jj),eion_row(jj)+ecou_row(jj),eele_row(jj)+epos_row(jj),
     2           erad_row(jj),ecou_row(jj)
            write(20,02) 'entr =',
     1           stot_row(jj),sion_row(jj),sion_row(jj)+scou_row(jj),sele_row(jj)+spos_row(jj),
     2           srad_row(jj),scou_row(jj)




c..   derivatives of the totals with respect to the input variables
            write(20,*)  ' '
            write(20,03) 'dp/dd=',dpd_row(jj),'dp/dt=',dpt_row(jj),
     1           'dp/da=',dpa_row(jj),'dp/dz=',dpz_row(jj)
            write(20,03) 'de/dd=',ded_row(jj),'de/dt=',det_row(jj),
     1           'de/da=',dea_row(jj),'de/dz=',dez_row(jj)
            write(20,03) 'ds/dd=',dsd_row(jj),'ds/dt=',dst_row(jj),
     1           'ds/da=',dsa_row(jj),'ds/dz=',dsz_row(jj)


c..   derivatives of the electron-positron compoenets with
c..   respect to the input variables
            write(20,*) ' ' 
            write(20,03) 'dpepd=',dpepd_row(jj),'dpept=',dpept_row(jj),
     1           'dpepa=',dpepa_row(jj),'dpepz=',dpepz_row(jj)
            write(20,03) 'deepd=',deepd_row(jj),'deept=',deept_row(jj),
     1           'deepa=',deepa_row(jj),'deepz=',deepz_row(jj)
            write(20,03) 'dsepd=',dsepd_row(jj),'dsept=',dsept_row(jj),
     1           'dsepa=',dsepa_row(jj),'dsepz=',dsepz_row(jj)


c..   the thermodynamic consistency relations, these should all be
c..   at the floating poiint limit of zero
            write(20,*) ' ' 
            write(20,03) 'maxw1=',dse_row(jj),'maxw2=',dpe_row(jj),
     1           'maxw3=',dsp_row(jj)


c..   number density of electrons, poistrons, matter electrons, and ions
            write(20,03) 'xne  =',xne_row(jj),'xnp  =',xnp_row(jj),
     1           'xnem =',xnem_row(jj),'xni  =',xni_row(jj)


c..   derivatibves of the electron number density with 
c..   respect to the input variables
            write(20,03) 'dxned=',dxned_row(jj),'dxnet=',dxnet_row(jj),
     1           'dxnea=',dxnea_row(jj),'dxnez=',dxnez_row(jj)


c..   electron chemical potential, positron chemical potential
c..   and derivatives of electron chemical potential with respect
c..   to the input variables
            write(20,03) 'eta  =',etaele_row(jj),'etap =',etapos_row(jj)
            write(20,03) 'detad=',detad_row(jj),'detat=',detat_row(jj),
     1           'detaa=',detaa_row(jj),'detaz=',detaz_row(jj)


c..   specific heats, and ratio of electostatic to thermal energy
            write(20,03) 'cp   =',cp_row(jj),'cv   =',cv_row(jj),
     1           'plasg=',plasg_row(jj)

c..   the 3 gammas and the sound speed
            write(20,03) 'gam1 =',gam1_row(jj),'gam2 =',gam2_row(jj),
     1           'gam3 =',gam3_row(jj),'csond=',cs_row(jj)
            write(20,*) ' '

         end if

        


cccccccccccc
c
c Below output only for LSEOS
c

         if((OUTPUTDETAILS.eq.1)
     1        .AND.(whichnucleareos.eq.1 .OR. whichnucleareos.eq.3)) then


c..   popular format statements
 11         format(1x,t11
     &           ,'total',t23,'nuclear',t41,'e+ e-',t59,'radiation')
 13         format(1x,t2,a,1pe12.4,3(1pe12.4,' (',i3,'%)'),a)
 16         format(1x,t2,a6,1pe12.5,a,t31,a6,1pe12.5,a)
 17         format(1x,t2,a6,1pe12.5,t22,a6,1pe12.5,t42,a6,1pe12.5)
 18         format(1x,t2,a6,1pe12.5,t22,a6,1pe12.5,t42,a6,1pe12.5,a6,1pe12.5)


c..   reflect the input
            write(20,*)
            write(20,*) 'LSEOS'
c            write(20,13)  whose
            write(20,*) 'input variables:'

            write(20,17) 'ye   =',ye_inp
            write(20,16) 'temp =',temp_cgs,' kelvin',
     &           'temp =',temp_nuc,' mev'
            write(20,16) 'den  =',den_cgs, ' g/cm**3',
     &           'den  =',den_nuc, ' 1/fm**3'


c..   give the output
            write(20,*)
            write(20,*)  'output composition:'

            write(20,18) 'xneut=',xnut,'xprot=',xprot,'xalfa=',xalfa
            write(20,17) 'xheav=',xh,  'aheav=',a,    'zheav=',x*a
            write(20,17) 'abar =',abar,'zbar =',zbar,    'ye_chk=',zbar/abar
            write(20,17) '1-sum=',1.0d0-(xnut+xprot+xalfa+xh),'muhat=',muhat



            write(20,*)
            write(20,*)  'output thermodynamics and derivatives in cgs:'

            write(20,11)  
            write(20,13) 'pres =',ptot_cgs,
     &           pbulk_cgs,nint(100.0d0*pbulk_cgs/ptot_cgs),
     &           pele_cgs,nint(100.0d0*pele_cgs/ptot_cgs),
     &           prad_cgs,nint(100.0d0*prad_cgs/ptot_cgs),
     &           ' erg/cm**3'

            write(20,13) 'ener =',etot_cgs,
     &           ebulk_cgs,nint(100.0d0*ebulk_cgs/etot_cgs),
     &           eele_cgs,nint(100.0d0*eele_cgs/etot_cgs),
     &           erad_cgs,nint(100.0d0*erad_cgs/etot_cgs),
     &           ' erg/gr'

            write(20,13) 'entr =',stot_cgs,
     &           sbulk_cgs,nint(100.0d0*sbulk_cgs/stot_cgs),
     &           sele_cgs,nint(100.0d0*sele_cgs/stot_cgs),
     &           srad_cgs,nint(100.0d0*srad_cgs/stot_cgs),
     &           ' erg/(gr k)'

            write(20,13) 'entr =',stot_nuc,
     &           sbulk_nuc,nint(100.0d0*sbulk_nuc/stot_nuc),
     &           sele_nuc,nint(100.0d0*sele_nuc/stot_nuc),
     &           srad_nuc,nint(100.0d0*srad_nuc/stot_nuc),
     &           ' dimensionless'

            write(20,13) 'free =',ftot_cgs,
     &           fbulk_cgs,nint(100.0d0*fbulk_cgs/ftot_cgs),
     &           fele_cgs,nint(100.0d0*fele_cgs/ftot_cgs),
     &           frad_cgs,nint(100.0d0*frad_cgs/ftot_cgs),
     &           ' erg/gr'


 111        format(1x,a,1pe24.16)
            write(20,17) 'dp/dd=',dpdd_cgs,'dp/dt=',dpdt_cgs,'dp/dy=',dpdy_cgs
            write(20,17) 'de/dd=',dedd_cgs,'de/dt=',dedt_cgs,'de/dy=',dedy_cgs
            write(20,17) 'ds/dd=',dsdd_cgs,'ds/dt=',dsdt_cgs,'ds/dy=',dsdy_cgs
            write(20,17) 'df/dd=',dfdd_cgs,'df/dt=',dfdt_cgs,'df/dy=',dfdy_cgs


            write(20,*)
            write(20,*)  'output specific heats,
     1                    sound speed(cgs), gammas:'

            write(20,17) 'cp   =',cp,'cv   =',cv,'csond=',sound
            write(20,17) 'gam1 =',gam1,'gam2 =',gam2,'gam3 =',gam3

            write(20,*)
            write(20,*)  'output thermodynamic consistency relations',
     1           ' (should be near zero):'
            write(20,17) 'dse  =',dse,'dpe =',dpe,'dsp =',dsp


c     write(20,*)
c     write(20,*)  'output thermodynamics and derivatives in nuc units:'

c     write(20,11)  
c     write(20,13) 'pres =',ptot_nuc,
c     &                     pbulk_nuc,nint(100.0d0*pbulk_nuc/ptot_nuc),
c     &                     pele_nuc,nint(100.0d0*pele_nuc/ptot_nuc),
c     &                     prad_nuc,nint(100.0d0*prad_nuc/ptot_nuc),
c     &                     ' mev/fm**3'

c     write(20,13) 'ener =',etot_nuc,
c     &                     ebulk_nuc,nint(100.0d0*ebulk_nuc/etot_nuc),
c     &                     eele_nuc,nint(100.0d0*eele_nuc/etot_nuc),
c     &                     erad_nuc,nint(100.0d0*erad_nuc/etot_nuc),
c     &                     ' mev/baryon'

c     write(20,13) 'entr =',stot_nuc,
c     &                     sbulk_nuc,nint(100.0d0*sbulk_nuc/stot_nuc),
c     &                     sele_nuc,nint(100.0d0*sele_nuc/stot_nuc),
c     &                     srad_nuc,nint(100.0d0*srad_nuc/stot_nuc),
c     &                     ' dimensionless'

c     write(20,13) 'free =',ftot_nuc,
c     &                     fbulk_nuc,nint(100.0d0*fbulk_nuc/ftot_nuc),
c     &                     fele_nuc,nint(100.0d0*fele_nuc/ftot_nuc),
c     &                     frad_nuc,nint(100.0d0*frad_nuc/ftot_nuc),
c     &                     ' mev/baryon'

c     write(20,*)
c     write(20,17) 'dp/dd=',dpdd_nuc,'dp/dt=',dpdt_nuc,'dp/dy=',dpdy_nuc
c     write(20,17) 'de/dd=',dedd_nuc,'de/dt=',dedt_nuc,'de/dy=',dedy_nuc
c     write(20,17) 'ds/dd=',dsdd_nuc,'ds/dt=',dsdt_nuc,'ds/dy=',dsdy_nuc
c     write(20,17) 'df/dd=',dfdd_nuc,'df/dt=',dfdt_nuc,'df/dy=',dfdy_nuc

            write(20,*)



         end if






      enddo



      close(18)
      close(19)
      close(20)



c 100  format(40(1pe27.15,' '))
 100  format(40(E27.16E5,' '))


      return
      end






