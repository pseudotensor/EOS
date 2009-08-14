
C======================================================================
      subroutine outputkazheader()
      implicit none

c     Bring-in Kaz single global quantities

c     only  need below for lsoffset,fakelsoffset,fakeentropylsoffset:
      include 'kazeos.dek'
c     Parms:
      include 'kazeos.parms.dek'
C======================================================================
c     Loop quantities
      include 'kazeos.loopvars1.dek'
      include 'kazeos.loopparms.dek'
C======================================================================

c     Note that if ntdynorynu==1 or ntdynorye==1 then chosen value ignored, while for hcm to be ignored one must choose whichhcmmethod==0
c     Some checks
      if( (whichynumethod.eq.0 .OR. whichynumethod.eq.3).AND.(ntdynorynu.gt.1)) then
         write(*,*) 'whichynumethod==0/3 and ntdynorynu>1 makes no sense'
         stop
      end if

c     actually this is the full table now
c      if( (whichynumethod.eq.1).AND.(whichhcmmethod.eq.1)) then
c         write(*,*) 'whichynumethod==1 and whichhcmmethod==1
c     1        not likely desired'
c         stop 
c      end if

c     actually this is now the simple table
c      if( (whichynumethod.eq.1).AND.(ntdynorynu.le.1)) then
c         write(*,*) 'whichynumethod==1 and ntdynorynu<=1
c     1        not likely desired'
c         stop
c      end if

      if( (whichynumethod.eq.2).AND.(ntdynorynu.ne.1)) then
         write(*,*) 'whichynumethod==2 and ntdynorynu!=1
     1        makes no sense'
         stop
      end if


      if( (whichrnpmethod.eq.0).AND.(ntdynorye.gt.1)) then
         write(*,*) 'whichrnpmethod==0 and ntdynorye>1 makes no sense'
         stop
      end if

      if( (whichhcmmethod.eq.0).AND.(nhcm.gt.1)) then
         write(*,*) 'whichhcmmethod==0 and nhcm>1 makes no sense'
         stop
      end if


c     Output header file
      open(50,file='eos.head',status='unknown')
      write(50,*) whichrnpmethod,whichynumethod,whichhcmmethod
c     First row indicates:
c     1) which data output type (column stuff same, then same number)
c     2) number independent variables
c     3) number of primary data columns
c     4) number of auxillary data columns
c     5) out of primary data, how many "extra" variables
c     8 original values (5 indeps and 3 vars)
      if(whichynumethod.eq.0 .OR. whichynumethod.eq.3) then
         if(whichrnpmethod.eq.0) then
            write(50,*) '1 4 9 80 1'
         else if(whichrnpmethod.eq.1) then
            if(whichhcmmethod.eq.0) then
               write(50,*) '2 3 24 80 16'
            else 
               write(50,*) '3 4 21 80 13'
            end if
         end if
      else if(whichynumethod.eq.1 .OR. whichynumethod.eq.2) then
         if(whichrnpmethod.eq.0) then
            write(50,*) '1 5 9 78 1'
         else if(whichrnpmethod.eq.1) then
            if(whichhcmmethod.eq.0) then ! present default method and most tested
c     GODMARK: 23 not all really extras since need derivatives of them
               write(50,*) '4 4 32 80 24'
            else
               write(50,*) '3 5 21 80 13'
            end if
         end if
      end if

c     Ranges of quantities outputted
      write(50,102) nrhob,rhobmin,rhobmax
      write(50,102) ntk,tkmin,tkmax
      write(50,102) ntdynorye,tdynoryemin,tdynoryemax
      write(50,102) ntdynorynu,tdynorynumin,tdynorynumax
      write(50,102) nhcm,hcmmin,hcmmax

c     Output other quantities
      write(50,103) lsoffset,fakelsoffset,fakeentropylsoffset
      close(50)


 102  format(I4,1x,1pe30.20E5,1x,1pe30.20E5,1x,1pe30.20E5,1x)
 103  format(1pe30.20E5,1x,1pe30.20E5,1x,1pe30.20E5)

      return
      end


c 18         format(1x,t2,a6,1pe12.5,t22,a6,1pe12.5,t42,a6,1pe12.5,a6,1pe12.5)




C======================================================================
      subroutine setupkazloops
     1     (intrhoblocal,inttklocal,inthcmlocal,inttdynoryelocal,inttdynorynulocal)
      implicit none

c     To be returned:
      real*8 intrhoblocal,inttklocal,inthcmlocal,inttdynoryelocal,inttdynorynulocal

c     Bring-in Kaz single global quantities
c      include 'kazeos.dek'
      include 'kazeos.parms.dek'
C======================================================================
c     Loop quantities
      include 'kazeos.loopvars1.dek'
      include 'kazeos.loopparms.dek'
C======================================================================


      if(nhcm.ne.1) then
         inthcmlocal=log10(hcmmax/hcmmin)/dble(nhcm-1)
      else
         inthcmlocal=0.d0
      end if

      if (ntdynorynu.ne.1) then
         inttdynorynulocal=log10(tdynorynumax/tdynorynumin)/dble(ntdynorynu-1)
      else
         inttdynorynulocal=0.d0
      end if


      if (ntdynorye.ne.1) then
         inttdynoryelocal=log10(tdynoryemax/tdynoryemin)/dble(ntdynorye-1)
      else
         inttdynoryelocal=0.d0
      end if
      


      if (ntk.ne.1) then
         inttklocal=log10(tkmax/tkmin)/dble(ntk-1)
      else
         inttklocal=0.d0
      end if
      
      if (nrhob.ne.1) then
         intrhoblocal=log10(rhobmax/rhobmin)/dble(nrhob-1)
      else
         intrhoblocal=0.d0
      end if


      return
      end





C======================================================================
      subroutine setkazindeps(
     1     irhob,itk,ihcm,itdynorye,itdynorynu
     1     ,intrhob,inttk,inthcm,inttdynorye,inttdynorynu
     1     ,rhob,tk,hcm,tdynorye,tdynorynu)
      implicit none
      save

      integer irhob,itk,ihcm,itdynorye,itdynorynu
      real*8 intrhob,inttk,inthcm,inttdynorye,inttdynorynu
      real*8 rhob,tk,hcm,tdynorye,tdynorynu

c     Bring-in Kaz single global quantities
c      include 'kazeos.dek'
      include 'kazeos.parms.dek'
C======================================================================
c     Loop quantities
      include 'kazeos.loopvars1.dek'
      include 'kazeos.loopparms.dek'
C======================================================================
c      integer firsttime
c      data firsttime/1/


c      if(firsttime.eq.1) then
c         firsttime=0
c         open(15,file='eoshelm.debug',status='old',access='append')
c      end if


c     These choices can be overrided by a calling code that doesn't use setkazindeps() or uses it but changes values before rest of code called

      if(whichynumethod.eq.0 .OR. whichynumethod.eq.3 .OR. whichynumethod.eq.2) then
c     Choose dummy hcm of 1 since hcm scales out of problem
         tdynorynu = 1.0
      else
         tdynorynu=tdynorynumin*10.d0**(inttdynorynu*dble(itdynorynu-1))
      end if

      if(whichrnpmethod.eq.0) then
c     Choose dummy hcm of 1 since hcm scales out of problem
         tdynorye = 1.0
      else
         tdynorye=tdynoryemin*10.d0**(inttdynorye*dble(itdynorye-1))
      end if

c      write(15,105) tdynoryemin,inttdynorye,itdynorye,tdynorye

      
      if(whichhcmmethod.eq.0) then
c     Choose dummy hcm of 1 since hcm scales out of problem

c     DEBUG:
c         hcm = 36497750.345734d0

         hcm = 1.0d0
      else
         hcm=hcmmin*10.d0**(inthcm*dble(ihcm-1))
      end if
      
      
      tk=tkmin*10.d0**(inttk*dble(itk-1))
      
      rhob=rhobmin*10.d0**(intrhob*dble(irhob-1))


c 105  format(E30.20E5,1x,E30.20E5,1x,I2,1x,E30.20E5)

c 102  format(I4,1x,1pe30.20E5,1x,1pe30.20E5,1x,1pe30.20E5,1x)

      return
      end




C======================================================================
      subroutine kaz_output()
      implicit none

c     Bring-in Kaz single global quantities
      include 'kazeos.dek'
      include 'kazeos.parms.dek'
C======================================================================
c     Loop quantities
      include 'kazeos.loopvars.dek'
      include 'kazeos.loopparms.dek'
C======================================================================


c     Open file for appended writing.  Assumes outside function
c     first opened it
      open(11,file='eos.dat',status='old',access='append')
      open(17,file='eosother.dat',status='old',access='append')



      if(whichynumethod.eq.0 .OR. whichynumethod.eq.3) then
         if(whichrnpmethod.eq.0) then
c     Here tdynorye is TDYN
            
            write(11,100) rhob, tk, tdynorye, tdynorynu, hcm !eos.dat
     |           ,p_tot, u_tot, s_tot, Qm
            
         else if(whichrnpmethod.eq.1) then
c     Here tdynorye is Y_e

            if(whichhcmmethod.eq.0) then
               
c     In this case, hcm scales out of problem and one can
c     reconstruct final Qm from tau/H's
c     In this case all quantities are independent of H
               write(11,100) rhob, tk, tdynorye,tdynorynu, hcm !eos.dat
     |              ,p_tot-p_nu, u_tot-rho_nu, s_tot-s_nu
     |              ,qtautelohcm, qtauaelohcm
     |              ,qtautmuohcm, qtauamuohcm
     |              ,qtauttauohcm, qtauatauohcm
     |              ,ntautelohcm, ntauaelohcm
     |              ,ntautmuohcm, ntauamuohcm
     |              ,ntauttauohcm, ntauatauohcm
     |              ,gammapeglobal+gammaAeglobal,gammapnuglobal+gammapenuglobal
     |              ,gammanglobal + gammaneglobal,gammannuglobal
            else
c     In this case HCM is assumed fully in calculations and HCM will be extra dimension
               write(11,100) rhob, tk, tdynorye,tdynorynu, hcm !eos.dat
     |              ,p_tot, u_tot, s_tot
     |              ,Qphoton, Qm, graddotrhouye, Tthermaltot, Tdifftot
     |              ,lambdatot,lambdaintot
     |              ,Enuglobal,Enueglobal,Enuebarglobal
     |              ,Ynuthermal,Ynu,Ynu0

            end if
         end if
      else if(whichynumethod.eq.1 .OR. whichynumethod.eq.2) then
         if(whichrnpmethod.eq.0) then
c     Here tdynorye is TDYN
            
            write(11,100) rhob, tk, tdynorye, tdynorynu, hcm !eos.dat
     |           ,p_tot, u_tot, s_tot, Qm
            
         else if(whichrnpmethod.eq.1) then
c     Here tdynorye is Y_e

            if(whichhcmmethod.eq.0) then
               
c     In this case, hcm scales out of problem and one can
c     reconstruct final Qm from tau/H's
c     In this case all quantities are independent of H
               write(11,100) rhob, tk, tdynorye,tdynorynu, hcm !eos.dat
     |              ,p_tot-p_nu, u_tot-rho_nu, s_tot-s_nu
     |              ,qtautnueohcm, qtauanueohcm
     |              ,qtautnuebarohcm, qtauanuebarohcm
     |              ,qtautmuohcm, qtauamuohcm
     |              ,ntautnueohcm, ntauanueohcm
     |              ,ntautnuebarohcm, ntauanuebarohcm
     |              ,ntautmuohcm, ntauamuohcm
     |              ,unue0,unuebar0,unumu0
     |              ,nnue0,nnuebar0,nnumu0
     |              ,lambdatot,lambdaintot
     |              ,tauphotonohcm,tauphotonabsohcm
     |              ,nnueth0,nnuebarth0

            else
c     In this case HCM is assumed fully in calculations and HCM will be extra dimension
               write(11,100) rhob, tk, tdynorye,tdynorynu, hcm !eos.dat
     |              ,p_tot, u_tot, s_tot
     |              ,Qphoton, Qm, graddotrhouye, Tthermaltot, Tdifftot
     |              ,lambdatot,lambdaintot
     |              ,Enuglobal,Enueglobal,Enuebarglobal
     |              ,Ynuthermal,Ynu,Ynu0

            end if
         end if
      end if





      
c     npratiofree is n_n/n_p for only free nucleons
      write(17,100) etae, etap, etan, etanu, xnuc, npratiofree   !eosother.dat
     |     ,p_tot,p_photon, p_eleposi, p_N, p_nu
     |     ,u_tot,rho_photon, rho_eleposi, rho_N, rho_nu
     |     ,s_tot,s_photon, s_eleposi, s_N, s_nu
     |     ,qtautelohcm, qtauaelohcm
     |     ,qtautmuohcm, qtauamuohcm
     |     ,qtauttauohcm, qtauatauohcm
     |     ,ntautelohcm, ntauaelohcm
     |     ,ntautmuohcm, ntauamuohcm
     |     ,ntauttauohcm, ntauatauohcm
     |     ,qtausohcm
     |     ,Qmel, Qmmu, Qmtau ! 37
     |     ,qminusel, qminusmu, qminustau
     |     ,ntausohcm
     |     ,Nmel, Nmmu, Nmtau
     |     ,nminusel, nminusmu, nminustau
     |     ,Qphoton, tauphotonohcm,tauphotonabsohcm ! 50
     |     ,yetot,yefree,yebound,dyedt,dyedtthin,graddotrhouyenonthermal,graddotrhouye,thermalye
     |     ,gammapeglobal,gammaAeglobal,gammapnuglobal,gammapenuglobal
     |     ,gammanglobal,gammaneglobal,gammannuglobal
     |     ,gammap2nglobal,gamman2pglobal
     |     ,Qm,Nm,Tdifftot,Tthermaltot,lambdatot
     |     ,Enuglobal,Enueglobal,Enuebarglobal
     |     ,RufNm,RufQm,Rufgraddotrhouye
     |     ,Ynu,Ynu0


c     JCM: without E5 a very small number like 1E-100 will print out as 1-100
 100  format(100(E30.20E5,' '))
c 100  format(100(E24.18,' '))

      close(11)
      close(17)


      return
      end


      
