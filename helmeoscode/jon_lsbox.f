cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Code outline:
c
c     1) any_eos()
c     2) Using nuclear EOS call jonnucbox(), else call full_nonnuclear_eos()
c     3) jonnucbox() 
c     --- sets: lsindex=jhi_eos+1 and loops over index=jlo_eos to jhi_eos
c     --- sets: yelocal,tklocal,rhoblocal used for lsbox() and shenbox()
c     --- calls check_convergence()
c     --- If goodconverge, then call store_row_fromcgsnuc2helmeos(index),
c     ------else set loci=lsindex=jhi_eos+1, save jlo_eos, jhi_eos, set as loci
c     ------udi0: Set den,temp,etc.(loci<-index), i.e. revert A,Z and call  full_nonnuclear_eos(index)
c     ------udi2: Set loci<-index only for den,temp, set abar,zbar,etc.(loci,index<-LSEOS singles)c     ------udi3: (see code, but extend use of A,Z,etc. revert den,temp, but match energy/baryon
c     --- Revert jlo_eos, jhi_eos
c                    
c
c     lsbox():
c     1) Set iunit,ye_inp,temp_cgs,den_cgs,usediffinputs,limitedrange,temp_cgs_lseos,den_cgs_lseos, temp_nuc,temp_nuc_lseos,den_nuc,den_nuc_lseos
c     2) Call lsguess, set ipvar(), call inveos(ye_inp,temp_nuc_lseos,den_nuc_lseos)
c     2.5) Internally calls any_electron()
c     3) Call store_row_fromlseos2cgs(didconverge)
c     --- Sets lsabar,etc. <--abar,etc.
c     --- Contains lsoffset
c     --- convert to ???_cgs single globals, also computes xnuc,npratiofree,etc. single globals
c     
c     If lsbox returns didconverge, then check convergence.  If goodconverge, then:
c     --- Call store_row_fromcgsnuc2helmeos(index)
c     ------ Stores abar_row(loci) <- lsabar,etc.
c     ------ etaele_row(loci) <- musube
c     ------ sets xnuc_row(loci) <- xnuc, etc.
c     
c     
c Note on mb vs. mn or whatever rest-mass definition
c
c LSEOS defines "internal" energy to be w.r.t. m_n, not m_u
c But LSEOS computes n_b = rhob/m_u (den_nuc)
c So for LSEOS, total true energy per baryon is $rhotrue/n_b = m_n + utot/nb$
c So $rhotrue/n_b = m_n - m_b + rhob/n_b + utot/n_b$
c So $rhotrue/n_b = rhob/n_b + [utot/n_b + m_n - m_b]$
c So $rhotrue = rhob + [utot + (m_n - m_b)*n_b]$
c So there is an additional effective internal energy that needs to be added for sound speed and self-gravity to be correct
c
c
c
c
c
cccccccccc List of subroutines:
c
cccccccccccccccccccccccccccccccccccc
c  Primary function called externally to get full EOS:
c      subroutine anyeos()
c
cccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccc
c  Set how extrapolation should be done and limit rho,T,Ye based upon whether within nuclear EOS or also if NSE holds
c      subroutine tableminmaxfixes(tmin, tmax, rhomin, rhomax, ypmin, ypmax)
cccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccc
c  Wrapper to call LS and Shen EOSs
c      subroutine jonnucbox(numconverged,numgoodconverged)
c  Wrapper for calling LS EOS
c      subroutine lsbox(iunitlocal, yelocal, tklocal, rhoblocal, didconverge)
c      subroutine loadguess
c      subroutine lsguess(rho,temp,ye,etanlin,etaplin,yplin)
cccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccc
c  Routines that check convergence and accuracy of general and specific checks for LS EOS
c      subroutine check_convergence(didconverge,goodconverge,numconverged,numgoodconverged)
c      subroutine check_convergence_lseos(didconverge,goodconverge,numconverged,numgoodconverged)
cccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccc
c  Wrappers to call HELM, TIMMES, and old LS electron EOS
c      subroutine full_nonnuclear_eos(whichnonnucleareos, loci, dostore,dostorespecies)
c      subroutine any_electron(tin,yein,dnsin)
c      subroutine any_electron_new(tin,yein,dnsin)
c Some calculations that want that normally done in Kaz EOS so also do for HELM/TIMMES EOSs
c      subroutine helmtimmes_extracalculation(loci,whichnonnucleareos)
c
c  Wrapper to call Kaz EOS for its electron or simple nuclear EOS:
c      subroutine kazeoswrapper()
c Prepare to call Kaz for either EOS or for neutrinos at end when writing to file:
c      subroutine preparecall2kazeos(loci)
c  After calling kazeoswrapper(), stores Kaz variables into normal EOS variables
c      subroutine store_row_kazeos2helminternals()
c  Convert species information from normal EOS to Kaz EOS format
c      subroutine eos2kazeos_species(nbtotal)
cccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccc
c Initialize EOS variables (not arrays) used by HELM/TIMMES:
c      subroutine init_row()
c Initialize EOS variables (not arrays) used by LS EOS
c      subroutine zero_cgsnuc_quantities()
c
c Store EOS results from individual variables to array location of similar name
c      subroutine store_row(loci,whichnonnucleareos,storenewspecies)
c
c Store LSEOS results (_nuc) into _cgs variable names (not arrays) [includes energy/baryon offset!]
c      subroutine store_row_fromlseos2cgs(didconverge)
c Calles store_row() to get the ele EOS and then LSEOS/Shen EOS results in _cgs/_nuc are put into HELM (normal) EOS *arrays* of similar names
c      subroutine store_row_fromcgsnuc2helmeos(loci)
c Store "Species" information from LS EOS (and Shen EOS) variable names into *arrays* of simliar names
c      subroutine store_row_lsspecies(loci)
c
c Used before file writing in order to take arrays and store them back into individual variable names:
c      subroutine storeback_row(loci)
c
c
cccccccccccccccccccccccccccccccccccc








cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c JCM: I added anyeos() to call any eos
c
c     1) Call anyeos() having den_row(), temp_row(), abar_row(), abarnum_row(), zbar_row()
c     set for jlo_eos to jhi_eos
c
c     2) If non-nuclear EOS used as nuclear EOS, then call full_nonnuclear_eos() and done.
c
c     3) If nuclear EOS used, then call jonnucbox() and done.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine anyeos()
      implicit none
      save
      
c..   bring in the data structure for HELM or TIMMES EOS
      include 'eosparms.f'
      include 'vector_eos.dek'


c      variable to indicate how many converged solutions out of each 1D line
      integer numconverged,numgoodconverged



cccccccccccccccccccccccccccccccccccccc
c
c..   call HELM/TIMMES/KAZ EOS
c
ccccccccccccccccccccccccccccccccccccccc

      if( (whichnucleareos.eq.0).OR.(whichnucleareos.eq.2).OR.(whichnucleareos.eq.4)) then
c     1 means do store result
         call full_nonnuclear_eos(whichnucleareos,index,1,1)
      end if

cccccccccccccccccccccccccccccccccccccc
c
c.. call LS or Shen EOS
c
ccccccccccccccccccccccccccccccccccccccc

      if((whichnucleareos.eq.1).OR.(whichnucleareos.eq.3)) then

         
         numconverged=0
         numgoodconverged=0

c     Uses jlo_eos,jhi_eos,index, etc.
c         write(*,*) 'Calling jonnucbox'
         call jonnucbox(numconverged,numgoodconverged)
c..   Back-stores many things in ??_row(?) that sit in vector_eos.dek and extra_vector_sneos.dek

         totalnumconverged=totalnumconverged+numconverged
         totalnumgoodconverged=totalnumgoodconverged+numgoodconverged

c         write(*,*) totalnumconverged,totalnumgoodconverged

      end if
c     Otherwise use HELM EOS "nuclear" EOS (i.e. only call helmeos and not jonnucbox)




      return
      end







cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Further restrict and change rho,T limits (includes NSE fix)
c     Controls behavior of how to handle being outside desired use of nuclear table
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tableminmaxfixes(tmin, tmax, rhomin, rhomax, ypmin, ypmax)
      implicit none
      save
      
c..   does all the stuff needed for calling the lattimer-swesty eos
      
      
c..   bring in the SHEN data structures
      include 'vector_sheneos.dek'
c     contains single globals
      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
      include 'eosparms.f'
c     Constants
      include 'const.dek'


      double precision AA,BB,CC

      double precision tmin,rhomin,tmax,rhomax, ypmin, ypmax
      double precision tminnse,rhominnse,tmaxnse,rhomaxnse
      double precision taunse,taunselimit
      integer resettmin1,resettmin2,resettmin3



c     defaults:
      limitedrange=0

c     Controls which density and temperature to use for LSEOS vs. electron EOS
c     0 = let LS EOS fail outside range of validity
c     1 = fix LS EOS to boundary values and use them
c     2 = fix LS EOS abar,abarnum,zbar to boundary values, but then override nuclear solution by using normal non-boundary rho,T

c     Controls which density and temperature to use for LSEOS vs. electron EOS
c      usediffinputs=2
c      usediffinputs=1
c      usediffinputs=1
c      usediffinputs=3
c     Just use simpler approach of interpolation but correct for conflicting energy/baryon at same species between nuclear and non-nuclear EOSs
      usediffinputs=4



c     Below 0 indicates not using yet -- if ever
      if(usediffinputs.eq.1 .OR. usediffinputs.eq.2 .OR. usediffinputs.eq.3 .OR. usediffinputs.eq.4) then
c     Fix rhob and T near region where LSEOS is bad, but continue using normal rho,T for electron EOS


c     taunsefix now here:
c     Note that as rho\to 0, taunse\to 0 so NSE well-established
c     Note that as T\to \infty, taunse\to 0, so NSE well-established
c     So we use taunse to force maximum on density

         AA=179.7*1.0D9
         BB=39.0d0
         CC=0.2d0
         taunse=(den_cgs)**(CC)*dexp(AA/(temp_cgs)-BB)
         taunselimit=1.0

         resettmin1=0
         resettmin2=0
         resettmin3=0

         if(taunse.gt.taunselimit) then
c     Then NSE only established on second timescales, so don't use NSE table -- just revert to desired A,Z from stellar model and use simplified nuclear EOS.
c     Make tmin a bit larger than temp_cgs so that if(temp_cgs.lt.tmin) is forced to be true
c            tminnse=temp_cgs

c     Need to choose temperature for which NSE was broken at this density so can back-track to that temperature as the minimum allowed
            tminnse = AA/(BB+dlog(taunselimit/(den_cgs)**(CC)))

c     Make rhomax a bit smaller than den_cgs so that if(den_cgs.gt.rhomax) is forced to be true
c     Currently rhomax isn't used because always sufficiently in NSE
            rhomaxnse=den_cgs
            resettmin1=1


            if(tmin.lt.tminnse .AND. tmax.gt.tminnse) then
c     Then within nuclear EOS table but need even more strict limit on tmin
               tmin = tminnse
               resettmin2=1
            end if
            if(tmax.lt.tminnse) then
c     Then entire table is not in NSE
               resettmin3=1
            end if

c     If tmin>tnse, then already outside table when not in NSE, so fine

c     write(*,*) 'taunsestuff',den_cgs,temp_cgs,taunse,taunselimit
c     write(*,*) 'taunsestuff2',tmin,tminnse
c            write(*,*) 'tmins',tmin,tminnse

         end if


c         write(*,*) 'resettmin',resettmin1,resettmin2,resettmin3


ccccccccccccccccccccc
c     Control Temperature
ccccccccccccccccccccc
         if( (temp_cgs.gt.tmax .OR. resettmin3.eq.1)) then
            limitedrange=1
            temp_cgs_lseos = tmax
         else if(temp_cgs.lt.tmin .OR. resettmin1.eq.1) then
            limitedrange=1
            temp_cgs_lseos = tmin
         else
            temp_cgs_lseos = temp_cgs
         end if

ccccccccccccccccccccc
c     Control Density
ccccccccccccccccccccc
         if(den_cgs.gt.rhomax) then
            limitedrange=1
            den_cgs_lseos = rhomax
         else if(den_cgs.lt.rhomin) then
            limitedrange=1
            den_cgs_lseos = rhomin
         else
            den_cgs_lseos = den_cgs
         end if
         

      end if



ccccccccccccccccccccc
c     Control Y_p
ccccccccccccccccccccc

c     Just don't allow ye to go out of range
      if(ye_inp.lt.ypmin .OR. ye_inp.gt.ypmax) then
         write(*,*) 'ye_inp out of range',ye_inp
         stop
      end if



c     Override even inside LSEOS domain
      if(den_cgs.lt.rhotaunsefix .AND. temp_cgs.lt.temptaunsefix) then
c     No longer use this method to fix non-NSE
c     GODMARK:
         taunsefix=0
      else
         taunsefix=0
      end if




      return
      end









cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c JCM: I added jonnucbox() as wrapper to wrapper called lsbox() or shenbox()
c
c
c Thompson et al. (2003) way to merge Helmholtz EOS and LS EOS
c rho<4E7g/cc use Helmholtz
c rho>6E7g/cc use LS between quadratically interpolate
c pressure difference order 1%, entropy order 5-10%.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine jonnucbox(numconverged,numgoodconverged)
       implicit none
       save


c..does all the stuff needed for calling the lattimer-swesty eos


c..bring in the lattimer-swesty data structures
c      include 'eos_m4c.commononly.inc'
c      include 'el_eos.inc'
c..bring in the data structure for LSEOS
c  contains single globals
      include 'eosparms.f'
      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
c      include 'const.dek'
c..bring in the data structure for HELM or TIMMES EOS
      include 'vector_eos.dek'



c     passed variable to indicate how many converged answers (out of jlo_eos to jhi_eos)
      integer numconverged,numgoodconverged


c     Local varaibles
      integer fortindex
      integer didconverge,goodconverge

      integer saved_jlo_eos,saved_jhi_eos
      integer loci

      integer iunitlocal
      double precision yelocal,tklocal,rhoblocal

      double precision localpdiff,localediff,localsdiff

      double precision difffactor
      double precision difffactor1,difffactor2,difffactorfinal
      double precision yekeep

c   Index used by LSEOS so can preserve original HELM EOS if LS EOS fails to converge
      lsindex=jhi_eos+1


ccccccccccccccccc
c     No longer use this method to fix non-NSE
      if(1.eq.0) then
c     
c     taunsefix
c     r= 10**8.186 cm taunse=1sec
c     
c     rhob<10**7.69549 in star
c     temp<10**9.76274=5.791e+09 in star
c     
c     Below these the NSE timescale is 1-10seconds, so don't have to strictly use nuclear EOSs values
c     rhotaunsefix=10**7.69549
c     rhotaunsefix=10**8.5
         rhotaunsefix=5.0D9
         temptaunsefix=1.0D10
      end if



c     initial star has rho=3.698e+09   temp=8.183e+09
cccccccccccccccccccccccccccccccccccccccccccccccc         
c
c     Loop over states to get LS EOS solution for each
c
c  NOTE that LSEOS only accepts 1 state at a time, so here loop over states and store results each time into new row vector for later use
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      do fortindex=jlo_eos,jhi_eos

c   Set index consistent with above (global variable used in any_electron())
         index=1+(fortindex-1)



cccccccccccccccccccccccccccccccccc
c
c     Set input to nuclear EOSs
c
ccccccccccccccccccccccccccccccccccc

         yelocal=zbar_row(index)/abarnum_row(index)
         tklocal=temp_row(index)
         rhoblocal=den_row(index)



cccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Get LS EOS solution
ccccccccccccccccccccccccccccccccccccccccccccccccc
         if(whichnucleareos.eq.1) then


c     iunit: Controls units used by LSEOS, where 1=cgs 2=nuclear
            iunitlocal=1

c     Uses iunit, ye_inp, temp_cgs, den_cgs
c..   This internally uses electron EOS, and default is to use already-computed values from HELM eos
c     write(*,*) "Calling lsbox",index,lsindex
            call lsbox(iunitlocal, yelocal, tklocal, rhoblocal, didconverge)

         end if


cccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Get Shen EOS solution
ccccccccccccccccccccccccccccccccccccccccccccccccc
         if(whichnucleareos.eq.3) then

c     write(*,*) "Calling shenbox",index,lsindex
            call shenbox(yelocal, tklocal, rhoblocal, didconverge)
c      DEBUG:
c            write(*,*) 'Called shenbox',yelocal,tklocal,rhoblocal,didconverge

         end if


ccccccccccccccccccccccccccccccccccc
c
c     See if and how converged
c
cccccccccccccccccccccccccccccccccccc
         call check_convergence(didconverge,goodconverge,numconverged,numgoodconverged)
c            write(*,*) didconverge,goodconverge,numconverged,numgoodconverged


c      DEBUG:
c         write(*,*) 'Called shenbox2',didconverge,goodconverge,numconverged,numgoodconverged



cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Store LSEOS/ShenEOS into HELM form
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
         if(goodconverge.eq.1) then
c..   NOW BACK-STORE RESULTS into row vectors replacing results from HELM EOS if was using it
c     write(*,*) "Calling store_row_fromlseos2helmeos",index,lsindex
            call store_row_fromcgsnuc2helmeos(index)
         end if






         if(goodconverge.eq.0) then
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  If no good convergence (or using limited range), then extrapolate into non-nuclear region (e.g. reduce to HELM/TIMMES/KAZ/etc. non-nuclear EOS)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc



c     If did not converge in good way, then avoid LSEOS and revert to electron EOS using newly computed abar and zbar
c     That is, convergence is not required to compute abar and zbar


c     "Loop" is over LSEOS extra storage for HELM/TIMMES EOS array
         lsindex=jhi_eos+1
         loci=lsindex

c     HELM and TIMMES EOSs use global variables already processed.  Save the ones that need to be changed before running then, so that they can be restored after that run
         saved_jlo_eos=jlo_eos
         saved_jhi_eos=jhi_eos
         jlo_eos = loci
         jhi_eos = loci


ccccccccccccccccccccc
c
c     HELM/TIMMES need den,temp,abar,abarnum,zbar _row versions set
c
ccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccc
c
c usediffinputs==1
c Hard switch to non-nuclear EOS
c
cccccccccccccccccccccccccccccccccccccccccccc

         
         if(usediffinputs.eq.0 .OR. didconverge.eq.0) then
cccccccccc STEP1 ccccccccccc

c     Get true rho,T
            den_row(loci)=den_row(index)
            temp_row(loci)=temp_row(index)

c     Reset A,Z to old-non-nuclear values and  hope that similar and no discontinuities
            zbar_row(loci)=zbar_row(index)
            abar_row(loci)=abar_row(index)
            abarnum_row(loci)=abarnum_row(index)

cccccccccc STEP2 ccccccccccc
c     Get non-nuclear solution and store all results into final solution ???_row(index) that is true solution states below lsindex
            call full_nonnuclear_eos(whicheleeos,index,1,1)


         end if


cccccccccccccccccccccccccccccccccccccccccccc
c
c usediffinputs==2
c Use original A,Z for all non-nuclear regimes
c Not good idea in general, and there will be jumps
c
cccccccccccccccccccccccccccccccccccccccccccc

         if(usediffinputs.eq.2  .AND. didconverge.eq.1) then
cccccccccc STEP1 ccccccccccc
c     Select nuclear EOS version of A,Z only and hope that energy/baryon, etc. are continuous

c     Get true rho,T
            den_row(loci)=den_row(index)
            temp_row(loci)=temp_row(index)
c     At least use consistent zbar,abar
            call store_row_lsspecies(loci)
c     This is the one set of things we store back overwriting original choice by non-nuclear EOS
            call store_row_lsspecies(index)

c            write(*,*) 'TEMPDEN',temp_row(loci),den_row(loci)
c            write(*,*) 'LSBAR',lszbar,lsabar,lsabarnum

cccccccccc STEP2 ccccccccccc
c     Get non-nuclear solution (store in lsindex space)
            call full_nonnuclear_eos(whicheleeos,lsindex,1,1)
c     Now store pure electron quantities back to index space (avoids overwriting xnuc, etc. species info)
            call store_row(index,whicheleeos,0)


         end if


cccccccccccccccccccccccccccccccccccccccccccc
c
c usediffinputs==3
c Force non-nuclear energy per baryon to match nuclear energy/baryon
c Note that species for nuclear not right if usingn old taunse, so jumps can occur that cannot fix
c This is why went to usediffinputs==4
c
cccccccccccccccccccccccccccccccccccccccccccc



         if(usediffinputs.eq.3  .AND. didconverge.eq.1) then
c     Select nuclear EOS version of A,Z and obtain non-nuclear energy/baryon, etc. so can determine correct offset

c     Setup lsindex version of density and temperature from limited nuclear table
            den_row(loci)=den_cgs_lseos
            temp_row(loci)=temp_cgs_lseos

cccccccccc STEP1 ccccccccccc
c     Store nuclei terms from nuclear EOS at nuclear version of temperature and density
c     Get in terms of per baryon
            localpdiff = pbulk_cgs/den_row(loci)
            localediff = ebulk_cgs
            localsdiff = sbulk_cgs

cccccccccc STEP2 ccccccccccc
c     Setup electron EOS inputs for A and Z to have consistent zbar,abar,density,temperature as nuclear EOS

            call store_row_lsspecies(loci)
c     This is the one set of things we store back overwriting original choice by non-nuclear EOS
            call store_row_lsspecies(index)
            

cccccccccc STEP3 ccccccccccc
c     Get non-nuclear EOS solution for nuclei at nuclear version of temperature and density
c     Store results in lsindex so don't overwrite nuclear species information
c            write(*,*) 'atincorrect',den_row(loci),temp_row(loci)
            call full_nonnuclear_eos(whicheleeos,lsindex,0,1)

cccccccccc STEP4 ccccccccccc
c     Store non-nuclear nuclei terms as part of offset (nuclear EOS includes Coulomb term in "ion")

c     Shouldn't matter whether include other non-nuclei terms
         if(1) then
            localpdiff = localpdiff - (pion+pcoul)/den_row(loci)
            localediff = localediff - (eion+ecoul)
            localsdiff = localsdiff - (sion+scoul)
         else
            localpdiff = localpdiff - (pres)/den_row(loci)
            localediff = localediff - (ener)
            localsdiff = localsdiff - (entr)
         end if

c         write(*,*) 'DIFFS',localpdiff,localediff,localsdiff

cccccccccc STEP5 ccccccccccc
c     Get non-nuclear EOS solution for nuclei at normal (correct) temperature and density
c     Again, store results in lsindex space so don't overwrite species information
         den_row(loci)=den_row(index)
         temp_row(loci)=temp_row(index)
c         write(*,*) 'atcorrect',den_row(loci),temp_row(loci)

         call full_nonnuclear_eos(whicheleeos,lsindex,0,1)



c     p and s don't have binding energy to shift and any mismatch is simply error
         localpdiff=0
c     GODMARK:
c         localediff=0
         localsdiff=0

         if(1.eq.0) then
c     DEBUG GODMARK

c     Interpolate pdiff from full value to 0 after 50% lower in density or temperature
c     den_row(loci)=den_cgs_lseos
c     temp_row(loci)=temp_cgs_lseos
c     den_row(loci)=den_row(index)
c     temp_row(loci)=temp_row(index)
            difffactor=(den_cgs_lseos-den_row(index))/den_cgs_lseos
c     write(*,*) 'difffactor1a',difffactor
            if(difffactor.lt.0.0) then
               difffactor=0.0
            else if(difffactor.gt.0.5) then
               difffactor=0.0
            else if(difffactor.ge.0.0 .AND. difffactor.le.0.5) then
               difffactor=(1.0-0.0)/(0.0-0.5)*(difffactor-0.5)
            end if

            localpdiff=localpdiff*difffactor
            localediff=localediff*difffactor
            localsdiff=localsdiff*difffactor
c     write(*,*) 'difffactor1b',difffactor

            difffactor=(temp_cgs_lseos-temp_row(index))/temp_cgs_lseos
c     write(*,*) 'difffactor2a',difffactor
            if(difffactor.lt.0.0) then
               difffactor=0.0
            else if(difffactor.gt.0.5) then
               difffactor=0.0
            else if(difffactor.ge.0.0 .AND. difffactor.le.0.5) then
               difffactor=(1.0-0.0)/(0.0-0.5)*(difffactor-0.5)
            end if

            localpdiff=localpdiff*difffactor
            localediff=localediff*difffactor
            localsdiff=localsdiff*difffactor
c     write(*,*) 'difffactor2b',difffactor
         end if

cccccccccc STEP6 ccccccccccc
c     Offset non-nuclear EOSs values so matches well to nuclear EOS at boundary where goodconverge==0
c     Only add it to ion term instead of Coulomb term
c     If nuclear EOS matched perfectly already with non-nuclear EOS, then correction cancels
c     As nuclear EOS has, say, a larger energy/baryon, then we add that to the "ion" term
         pion = pion + localpdiff*den_row(loci)
         eion = eion + localediff
         sion = sion + localsdiff

c     Have to correct totals and individuals since already totalled into pres,ener,entr inside non-nuclear EOS
         pres = pres + localpdiff*den_row(loci)
         ener = ener + localediff
         entr = entr + localsdiff

c         write(*,*) 'diffs',localpdiff,localediff,localsdiff

cccccccccc STEP7 ccccccccccc
c     Finally, store into ???_row @ index where normally would have been put
c     This only stores back non-nuclear nucleon+electron+total into index space, not overwriting nuclear species information
         call store_row(index,whicheleeos,0)




         end if



cccccccccccccccccccccccccccccccccccccccccccc
c
c usediffinputs==4
c Interpolation between nuclear and non-nuclear since generally jump in species
c But also account for energy/baryon offset between nuclear and non-nuclear EOS at boundary where nuclear EOS is used
c so that in general this is continuous so consistent
c Pressure and entropy not changed since assume any error is small
c Species change can make pressure and entropy jump, but nuclear pressure is not dominant where doing this interpolation
c
cccccccccccccccccccccccccccccccccccccccccccc


         if(usediffinputs.eq.4  .AND. didconverge.eq.1) then
c     Select nuclear EOS version of A,Z and obtain non-nuclear energy/baryon, etc. so can determine correct offset

c     Setup lsindex version of density and temperature from limited nuclear table
            den_row(loci)=den_cgs_lseos
            temp_row(loci)=temp_cgs_lseos

cccccccccc STEP1 ccccccccccc
c     Store nuclei terms from nuclear EOS at nuclear version of temperature and density
c     Get in terms of per baryon
            localpdiff = pbulk_cgs/den_row(loci)
            localediff = ebulk_cgs
            localsdiff = sbulk_cgs

cccccccccc STEP2 ccccccccccc
c     Setup electron EOS inputs for A and Z to have consistent zbar,abar,density,temperature as nuclear EOS
c            write(*,*) 'store_row_lsspecies',loci
            call store_row_lsspecies(loci)
c     This is the one set of things we store back overwriting original choice by non-nuclear EOS
c     Required for getting nuclear version of (e.g.) etapls,etanls that is needed by neutrinos later
c            call store_row_lsspecies(index)
c            etap_row(index) = etap_row(loci)
c            etan_row(index) = etan_row(loci)

           

cccccccccc STEP3 ccccccccccc
c     Get non-nuclear EOS solution for nuclei at nuclear version of temperature, density, and species
c     Store results in lsindex so don't overwrite nuclear species information
c            write(*,*) 'atincorrect',den_row(loci),temp_row(loci)
            call full_nonnuclear_eos(whicheleeos,lsindex,0,1)


cccccccccc STEP4 ccccccccccc
c     Store non-nuclear nuclei terms as part of offset (nuclear EOS includes Coulomb term in "ion")

c     Shouldn't matter whether include other non-nuclei terms
         if(1) then
            localpdiff = localpdiff - (pion+pcoul)/den_row(loci)
            localediff = localediff - (eion+ecoul)
            localsdiff = localsdiff - (sion+scoul)
         else
            localpdiff = localpdiff - (pres)/den_row(loci)
            localediff = localediff - (ener)
            localsdiff = localsdiff - (entr)
         end if

c         write(*,*) 'DIFFS',localpdiff,localediff,localsdiff

cccccccccc STEP5 ccccccccccc
c     Get non-nuclear EOS solution for nuclei at normal (correct) temperature and density and (correct) species
c     Again, store results in lsindex space so don't overwrite species information
         den_row(loci)=den_row(index)
         temp_row(loci)=temp_row(index)

         
         if(1.eq.0) then
c     Reset A,Z to old-non-nuclear values and  hope that similar and no discontinuities
            zbar_row(loci)=zbar_row(index)
            abar_row(loci)=abar_row(index)
            abarnum_row(loci)=abarnum_row(index)
         else
c     Interpolate A for a bit to avoid induced jumps in A (note, real star or whatever can have real compositional jumps)
            difffactor1=(den_cgs_lseos-den_row(index))/den_cgs_lseos
            difffactor2=(temp_cgs_lseos-temp_row(index))/temp_cgs_lseos
            difffactor=dmax1(difffactor1,difffactor2)


c     Allow for 50% change in T or rho before completely changing A
            if(difffactor.lt.0.0) then
               difffactorfinal=1.0
            else if(difffactor.gt.0.2) then
               difffactorfinal=0.0
            else if(difffactor.ge.0.0 .AND. difffactor.le.0.2) then
               difffactorfinal=(1.0-0.0)/(0.0-0.2)*(difffactor-0.2)
            end if

c     yekeep could be computed using loci or index, same ye
            yekeep=zbar_row(loci)/abarnum_row(loci)
            abar_row(loci)=abar_row(index)*(1.0d0-difffactorfinal) + abar_row(loci)*difffactorfinal
            abarnum_row(loci)=abarnum_row(index)*(1.0d0-difffactorfinal) + abarnum_row(loci)*difffactorfinal
            zbar_row(loci)=yekeep*abarnum_row(loci)
            zbar_row(index)=zbar_row(loci)
            abar_row(index)=abar_row(loci)
            abarnum_row(index)=abarnum_row(loci)

c            write(*,*) 'diffstuff',difffactor1,difffactor2,difffactor,difffactorfinal

         end if

c         write(*,*) 'atcorrect',den_row(loci),temp_row(loci)

c     abar,zbar,abarnum still old non-nuclear values since nuclear didn't converge so didn't call store_row_lsspecies()
c     The 0,1 means don't overwrite any _row's or species information (i.e. don't call store_row_lsspecies())
         call full_nonnuclear_eos(whicheleeos,lsindex,0,1)


c     p and s don't have binding energy to shift and any mismatch is simply error
         localpdiff=0
c     GODMARK:
c         localediff=0
         localsdiff=0

         if(1.eq.0) then
c     DEBUG GODMARK

c     Interpolate pdiff from full value to 0 after 50% lower in density or temperature
c     den_row(loci)=den_cgs_lseos
c     temp_row(loci)=temp_cgs_lseos
c     den_row(loci)=den_row(index)
c     temp_row(loci)=temp_row(index)
            difffactor=(den_cgs_lseos-den_row(index))/den_cgs_lseos
c     write(*,*) 'difffactor1a',difffactor
            if(difffactor.lt.0.0) then
               difffactor=0.0
            else if(difffactor.gt.0.5) then
               difffactor=0.0
            else if(difffactor.ge.0.0 .AND. difffactor.le.0.5) then
               difffactor=(1.0-0.0)/(0.0-0.5)*(difffactor-0.5)
            end if

            localpdiff=localpdiff*difffactor
            localediff=localediff*difffactor
            localsdiff=localsdiff*difffactor
c     write(*,*) 'difffactor1b',difffactor

            difffactor=(temp_cgs_lseos-temp_row(index))/temp_cgs_lseos
c     write(*,*) 'difffactor2a',difffactor
            if(difffactor.lt.0.0) then
               difffactor=0.0
            else if(difffactor.gt.0.5) then
               difffactor=0.0
            else if(difffactor.ge.0.0 .AND. difffactor.le.0.5) then
               difffactor=(1.0-0.0)/(0.0-0.5)*(difffactor-0.5)
            end if

            localpdiff=localpdiff*difffactor
            localediff=localediff*difffactor
            localsdiff=localsdiff*difffactor
c     write(*,*) 'difffactor2b',difffactor
         end if

cccccccccc STEP6 ccccccccccc
c     Offset non-nuclear EOSs values so matches well to nuclear EOS at boundary where goodconverge==0
c     Only add it to ion term instead of Coulomb term
c     If nuclear EOS matched perfectly already with non-nuclear EOS, then correction cancels
c     As nuclear EOS has, say, a larger energy/baryon, then we add that to the "ion" term
         pion = pion + localpdiff*den_row(loci)
         eion = eion + localediff
         sion = sion + localsdiff

c     Have to correct totals and individuals since already totalled into pres,ener,entr inside non-nuclear EOS
         pres = pres + localpdiff*den_row(loci)
         ener = ener + localediff
         entr = entr + localsdiff

c         write(*,*) 'diffs',localpdiff,localediff,localsdiff

ccccccccccSTEP7 ccccccccccc
c     Finally, store into ???_row @ index where normally would have been put
c     This only stores back non-nuclear nucleon+electron+total into index space, ALSO overwriting nuclear species information
c     call store_row(index,whicheleeos,0)
         call store_row(index,whicheleeos,1)
         
         call backstore_row_lsspecies(index)

c     Store "LS" results back into normal array
         call store_row_lsspecies(index)
         



         end if





c     write(*,*) 'TEMPDEN',temp_row(loci),den_row(loci)
c     write(*,*) 'LSBAR',lszbar,lsabar,lsabarnum


c     NOW restore saved versions
c     Now any EOS processes called before reaching here are preserved
         jlo_eos=saved_jlo_eos
         jhi_eos=saved_jhi_eos


         end if
c  Otherwise assume calling function keeps decent values







         if(didconverge_row(index).ge.2) then
            write(*,*) "Why here4?"
         end if




cccccccccccccccccccccccccccccc
c
c.. done looping over states for single-input function producing LSEOS
c
cccccccccccccccccccccccccccccc
      end do




      return
      end







c Should roughly be inverse of store_row_lsspecies()
c But should be consistent with how set abar and zbar initially in jon_helm.f or jon_helm_stellarinput.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine backstore_row_lsspecies(loci)
      implicit none
      save
      
      integer loci
      double precision nblocal




c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
c      include 'el_eos.inc'
c..bring in the data structure for LSEOS
c  contains single globals
      include 'vector_eos.dek'
      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
      include 'extra_vector_sneos.dek'
c      include 'const.dek'
c..bring in the data structure for HELM or TIMMES EOS
      include 'kazeos.dek'

c     Constants
      include 'const.dek'



c     We want to use store_row_lsspecies() below, but first assign right-hand-side inside there to be correct
         lszbar=zbar_row(loci)
         lsabar=abar_row(loci)
         lsabarnum=abarnum_row(loci)

         if(lsabar.gt.4.5) then
c     call computemutotfit_xnuc0(den_row(),temp_row(loci),abar_row(loci))
            xnut =1D-49
            xprot =1D-49
            xalfa = 1D-49
            xh = 1.0d0
            a = lsabar
            x = lszbar/lsabar
            muhat = lsabar

c     Note that etapls and etanls are set from "lsindex" unlike others, which means we use nuclear versions of these, where nuclear versions are assigned at the boundary of the table where nuclear EOS exists (goodconverge.eq.1)
            etapls = etap_row(lsindex)
            etanls = etan_row(lsindex)

            xnuc = xnut + xprot
            yetot = x
            yefree = 0.5d0      !doesn't matter
            yebound = x
            yeheav = x
            npratiofree = 1D-49
            npratiobound = 1.0d0-x
            npratiototal = npratiobound
            abarbound = lsabar

            nblocal = (den_row(loci)/mb)
            
            nptotal = nblocal/(1.0d0 + npratiototal)
            nntotal = nblocal - nptotal
            npfree = 1D-49
            nnfree = 1D-49
            npbound = nptotal
            nnbound = nntotal
            npheav = nptotal
            nnheav = nntotal

         else if(lsabar.gt.3.5) then
c     call computemutotfit_xnuc0(den_row(),temp_row(loci),abar_row(loci))
            xnut =1D-49
            xprot =1D-49
            xalfa = 1.0d0
            xh = 1D-49
            a = 1D-49
            x = 0.50d0 ! doesn't matter
            muhat = lsabar

c     Note that etapls and etanls are set from "lsindex" unlike others, which means we use nuclear versions of these, where nuclear versions are assigned at the boundary of the table where nuclear EOS exists (goodconverge.eq.1)
            etapls = etap_row(lsindex)
            etanls = etan_row(lsindex)

            xnuc = xnut + xprot
            yetot = x
            yefree = 0.5d0      !doesn't matter
            yebound = x
            yeheav = 0.5d0 ! doesn't matter
            npratiofree = 1D-49
            npratiobound = 1.0d0-x
            npratiototal = npratiobound
            abarbound = lsabar

            nblocal = (den_row(loci)/mb)
            
            nptotal = nblocal/(1.0d0 + npratiototal)
            nntotal = nblocal - nptotal
            npfree = 1D-49
            nnfree = 1D-49
            npbound = nptotal
            nnbound = nntotal
            npheav = 1D-49
            nnheav = 1D-49

         else
c     This extrapolation of etap and etan works unless in hydrogen region.  Just force etap=etan=0 if in that region where etae>>1
            xnut = 1D-49
            xprot =1
            xalfa = 1D-49
            xh = 1D-49
            a = 1D-49
            x = 0.5d0
            muhat = lsabar

            etapls = 1D-49
c     Assume beta equilibrium as valid at high temperatures
c     see Kohri & Mineshige (2002)
            if(temp_row(loci) .gt. me*light2/kerg) then
               etanls = etapls + etae
            end if

c     Except sometimes assume degenerate nucleons at high density and low temp
c     see Kohri & Mineshige (2002)
            if( (temp_row(loci) .lt. me*light2/kerg) .AND. (den_row(loci).gt. 1.0D10) ) then
               etanls = etae + (mn-mb)*light2
               etapls = etae + (mp-mb)*light2
            end if


            xnuc = xnut + xprot
            yetot = lszbar/lsabar
            yefree = lszbar/lsabar
            yebound = 1D-49
            yeheav = 0.5d0 !doesn't matter
            npratiofree = 1.0d0 - yefree
            npratiobound = 1D-49
            npratiototal = npratiofree
            abarbound = 1D-49

            nblocal = (den_row(loci)/mb)
            
            nptotal = nblocal/(1.0d0 + npratiototal)
            nntotal = nblocal - nptotal
            npfree = nptotal
            nnfree = nntotal
            npbound = 1D-49
            nnbound = 1D-49
            npheav = 1D-49
            nnheav = 1D-49

         end if



      


      return
      end







c Should roughly be inverse of preparecall2kazeos() except for independent variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine store_row_lsspecies(loci)
      implicit none
      save
      
      integer loci

c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
c      include 'el_eos.inc'
c..bring in the data structure for LSEOS
c  contains single globals
      include 'vector_eos.dek'
      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
      include 'extra_vector_sneos.dek'
c      include 'const.dek'
c..bring in the data structure for HELM or TIMMES EOS
      include 'kazeos.dek'

c     This is the one thing we store back overwriting original choice by non-nuclear EOS
      zbar_row(loci)=lszbar
      abar_row(loci)=lsabar
      abarnum_row(loci)=lsabarnum


c     NEW LSEOS quantities
      xneut_row(loci)    = xnut
      xprot_row(loci)    = xprot
      xalfa_row(loci)    = xalfa
      xheav_row(loci)    = xh
      aheav_row(loci)    = a
      zheav_row(loci)    = x*a
      xcheck_row(loci)   = 1.0d0-(xnut+xprot+xalfa+xh)
      muhat_row(loci)    = muhat


      if(aheav_row(loci).lt.0.0) then
         write(*,*) 'aheavneg',xnut,xprot,xalfa,xh,a,x,xcheck_row(loci)
      end if

      if(zheav_row(loci).lt.0.0) then
         write(*,*) 'zheavneg',xnut,xprot,xalfa,xh,a,x,xcheck_row(loci)
      end if

c     GODMARK: treat the below as species information?
c      write(*,*) 'etaprecheck',loci,etapls,etanls
      etap_row(loci) = etapls
      etan_row(loci) = etanls
      etanu_row(loci) = 0.0

c     Also store-over Kaz-like species quantities
      xnuc_row(loci)=xnuc
      yetot_row(loci)=yetot
      yefree_row(loci)=yefree
      yebound_row(loci)=yebound
      yeheav_row(loci)=yeheav
      npratiofree_row(loci)=npratiofree
      npratiobound_row(loci)=npratiobound
      npratioheav_row(loci)=npratioheav
      npratiototal_row(loci)=npratiototal
      abarbound_row(loci)=abarbound

      if(yebound_row(loci).lt.0.0) then
         write(*,*) 'yeboundneg',xnut,xprot,xalfa,xh,a,x,xcheck_row(loci)
      end if


      nptotal_row(loci)=nptotal
      nntotal_row(loci)=nntotal
      npfree_row(loci)=npfree
      nnfree_row(loci)=nnfree
      npbound_row(loci)=npbound
      nnbound_row(loci)=nnbound
      npheav_row(loci)=npheav
      nnheav_row(loci)=nnheav

      


      return
      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Call electron/positron/photon/nuclei EOS (non-nuclear)
c     whichnonnucleareos corresponds to whicheleeos numbers
c     
c     For HELM/TIMMES, uses jlo_eos,jhi_eos,index
c     to determine what memory location to look at
c     for arrays den_row(), temp_row(), abarnum_row(), abar_row(),zbar_row()
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine full_nonnuclear_eos(whichnonnucleareos, loci, dostore,dostorespecies)
      implicit none
      save
      

c     include 'const.dek'
c..bring in the data structure for HELM or TIMMES EOS
c      include 'vector_eos.dek'

c     Passed
      integer whichnonnucleareos,loci,dostore,dostorespecies


cccccccccccccccccccccccccccccccccccccc
c     
c..   call HELM EOS
c     
ccccccccccccccccccccccccccccccccccccccc
      if(whichnonnucleareos.eq.0) then
c..   call the eos with assumed global vectors of temp_row, den_row, abar_row, abarnum_row, zbar_row
c     write(*,*) 'Calling helmeos'
         call helmeos
c..   Stores many things in ??_row(?) that sit in vector_eos.dek

c     ISSUES:
c     HELM does very poorly for \eta_e at low density, high temperatures

      end if

cccccccccccccccccccccccccccccccccccccc
c     
c..   call TIMMES EOS
c     
ccccccccccccccccccccccccccccccccccccccc
      if(whichnonnucleareos.eq.2) then
c..   call the eos with assumed global vectors of temp_row, den_row, abar_row, abarnum_row, zbar_row
c     write(*,*) 'Calling eosfxt'
         call eosfxt
c..   Stores many things in ??_row(?) that sit in vector_eos.dek

c     ISSUES:
c     Sometimes EOSFXT gets negative entropy for electron/positron, while HELM does not.

      end if



cccccccccccccccccccccccccccccccccccccc
c     
c..   call Kaz EOS
c     
ccccccccccccccccccccccccccccccccccccccc
      if(whichnonnucleareos.eq.4) then
c..   call the eos with assumed global vectors of temp_row, den_row, abar_row, abarnum_row, zbar_row
c     write(*,*) 'Calling kazeoswrapper'
         call kazeoswrapper
c..   Stores many things in ??_row(?) that sit in vector_eos.dek

c     ISSUES:
c     Relatively slow
c     No Coulomb corrections
      end if



c     Now store from lsindex to normal index
      if(dostore.eq.1) then
         call store_row(loci,whichnonnucleareos,dostorespecies)
      end if





      return
      end







cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Wrapper for KAZ EOS
c
c     Uses global jlo_eos,jhi_eos to set where to store in HELM ???_row() arrays
c     For indep vars uses globals den_row, temp_row,abarnum_row,zbar_row,abar_row
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine kazeoswrapper()
      implicit none
      save
      
c..bring in the data structure for HELM or TIMMES EOS (only needs jlo_eos jhi_eos)
      include 'vector_eos.dek'
c     Bring in Kaz data structures
      include 'kazeos.dek'
      include 'eosparms.f'

      integer kazindex
      integer whichnonnucleareos
c      integer computespecies


c     First setup globals used by kaz_eos()

      do kazindex=jlo_eos,jhi_eos

         

         call preparecall2kazeos()


c         write(*,*) kazindex,jlo_eos,jhi_eos,rhob,tk,hcm,tdynorye,tdynorynu
c     Call Kaz EOS to get those things set in kazeos.dek

         if(whichnucleareos.eq.4) then
            computespecies=1
         else
c     Then assume nuclear EOS called and set species
c     Needed to have set xnuc, yetot, yefree, yebound, npratiofree, kazabar, kazzbar, etc.
            computespecies=0
         end if


         call kaz_eos()
         
         
         
c     Kaz quantities need to be put into HELM internals so can
c     treat like HELM or TIMMES EOS was called
         call store_row_kazeos2helminternals()

c     Act like HELM/TIMMES and store internals into ???_row() arrays for a given indexs
         whichnonnucleareos=4
         call store_row(kazindex,whichnonnucleareos,0)


      end do


      return
      end






cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Sets up call to kazeos by assigning variables NEEDED BY Kaz EOS, not necessarily generated by Kaz EOS
c     Should roughly be inverse of store_row_lsspecies() except for independent variables
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine preparecall2kazeos(loci)
      implicit none
      save

      integer loci
      
c..bring in the data structure for HELM or TIMMES EOS (only needs jlo_eos jhi_eos)
      include 'vector_eos.dek'
      include 'extra_vector_sneos.dek'
c     Bring in Kaz data structures
      include 'kazeos.dek'
      include 'eosparms.f'


      rhob=den_row(loci)
      tk=temp_row(loci)
      hcm=hcm_row(loci)
c     Note that this "recomputes" tdynorye.  Assume any error is in zbar or abar, so don't recompute this and assume same as set initially
c     Noticed that can  have small errors since abar,zbar are changed
c      tdynorye=zbar_row(loci)/abarnum_row(loci)
      tdynorye=tdynorye_row(loci)
      tdynorynu=tdynorynu_row(loci)


      kazabar=abar_row(loci)
      kazzbar=zbar_row(loci)

      kazxneut=xneut_row(loci)
      kazxprot=xprot_row(loci)
      kazxalfa=xalfa_row(loci)
      kazxheav=xheav_row(loci)
      kazaheav=aheav_row(loci)
      kazzheav=zheav_row(loci)
c      yeheav = (kazzheav/kazaheav) ! no _row, but set in kaz_species_etae()
c      xcheck_row() and muhat_row() are not needed by kaz EOS

c     If not getting Kaz EOS to compute this, then needed
c     If Kaz EOS will compute this, won't hurt to make this assignment
c     Below is ONLY time eta,etap,etan,etanu should be assigned or used since LS/Shen EOS set etapls,etanls that should elsewhere be used.
c     Note that etaele_row() is non-nuclear so doesn't appear in store_row_lsspecies()
      etae=etaele_row(loci)
      etap=etap_row(loci)
      etan=etan_row(loci)
      etanu=etanu_row(loci)


      npratiototal = (1.0-tdynorye)/tdynorye  ! no _row, but set in kaz_species_etae()

c     Rest are for overriding Kaz-set species in case using nuclear EOS
      xnuc=xnuc_row(loci)
      yetot=yetot_row(loci)
      yefree=yefree_row(loci)
      yebound=yebound_row(loci)
      yeheav=yeheav_row(loci)
      npratiofree=npratiofree_row(loci)
      npratiobound = npratiobound_row(loci)
      npratioheav = npratioheav_row(loci)
      abarbound = abarbound_row(loci)

c     DEBUG:
c      write(*,*) 'kazazbar',kazabar,kazzbar,abarbound



c      write(*,*) 'etacheck',loci,etae,etap,etan,etanu

      nptotal = nptotal_row(loci)
      nntotal = nntotal_row(loci)
      npfree = npfree_row(loci)
      nnfree = nnfree_row(loci)
      npbound = npbound_row(loci)
      nnbound = nnbound_row(loci)
      npheav = npheav_row(loci)
      nnheav = nnheav_row(loci)

c DEBUG if things are set or not
c      write(*,*) loci,rhob,tk,hcm,tdynorye,tdynorynu,npratiototal,xnuc
c     1 ,yetot,yefree,yebound,npratiofree,npratiobound,npratioheav
c     1 ,kazabar,kazzbar,abarbound,kazxneut,kazxprot,kazxalfa
c     1 ,kazxheav,kazaheav,kazzheav,yeheav,etae,etap,etan,etanu
c     1 ,nptotal,nntotal,npfree,nnfree,npbound,nnbound,npheav,nnheav


      
      return
      end







cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Check convergence properties of LSEOS and Shen EOS        
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_convergence(didconverge,goodconverge,numconverged,numgoodconverged)
      implicit none
      save
      
      include 'eosparms.f'
      include 'vector_eos.dek'
      include 'vector_sneos.dek'

c     Passed or returned variables
      integer didconverge,goodconverge,numconverged,numgoodconverged



      if(whichnucleareos.eq.1) then
         
c     See if and how converged
         call check_convergence_lseos(didconverge,goodconverge,numconverged,numgoodconverged)
         
      else if(whichnucleareos.eq.3) then
         
c     With Shen EOS, output is predetermined to be good if within lookup table
c     If not within lookup table, treat as if failed just as if lsbox generated goodconverge=0
         if(didconverge.eq.1) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     if usediffinputs==2 or 3, then setting only A,Z to boundary values, and when detect fixing to boundary, enforce that solution was googconverge==0
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            numconverged=numconverged+1
            didconverge_row(index)=didconverge

            if( (usediffinputs.eq.2 .OR. usediffinputs.eq.3 .OR. usediffinputs.eq.4) .AND. limitedrange.eq.1) then
               goodconverge=0
c     Using didconverge single value to indicate success of limited range, so can still use _row version to indicate approach
               didconverge_row(index)=-500
            else
               goodconverge=1
               numgoodconverged=numgoodconverged+1
            end if
         else
            didconverge_row(index)=didconverge
            goodconverge=0
         end if
         
      end if


      



      return
      end







cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Check convergence properties of LSEOS        
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine check_convergence_lseos(didconverge,goodconverge,numconverged,numgoodconverged)
       implicit none
       save


c  contains single global: xnut, and similar things
c..bring in the lattimer-swesty data structures
      include 'eosparms.f'
      include 'eos_m4c.commononly.inc'
c..bring in the data structure for LSEOS
c  contains single globals
      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
      include 'const.dek'
c..bring in the data structure for HELM or TIMMES EOS
      include 'vector_eos.dek'

c     Local varaibles
      integer consistent,positivepressure,positiveenergy,yeconsistent,xconsistent,cvcpcheck,muhatcheck
      double precision xcheck

c     Passed or returned variables
      integer didconverge,goodconverge,numconverged,numgoodconverged

      double precision yelocal

      double precision negative1,mynan,cpinf



         if((abs(dse).lt.LSTOL)
     2      .AND.(abs(dpe).lt.LSTOL)
     3      .AND.(abs(dsp).lt.LSTOL)) then
            consistent=1
         else
            consistent=0
         end if

c     Override, assume user knew what they were doing
         if(usediffinputs.eq.1 .OR. usediffinputs.eq.2 .OR. usediffinputs.eq.3 .OR. usediffinputs.eq.4) then
            consistent=1
         end if

         if(ptot_cgs.gt.0.0) then
            positivepressure=1
         else
            positivepressure=0
         end if

c     LSEOS/SHEN EOS ok with etot<0
         if((1).OR.(etot_cgs.gt.0.0)) then
            positiveenergy=1
         else
            positiveenergy=0
         end if

c     if(abs(zbar_row(lsindex)/abarnum_row(lsindex) - ye_inp)<1E-5) then
         if(abs(lszbar/lsabarnum - ye_inp)<1E-5) then
            yeconsistent=1
         else
            yeconsistent=0
            call enforce_ye_consistency(ye_inp, xnut, xprot, xalfa, xh, a, x, yeconsistent)
         end if

ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        Recompute these since adjusted above
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
         call compute_nuclear_azbar(xnut, xprot, xalfa, xh, a, x,
     1        abarnum,abar,zbar,yelocal)
         lszbar = zbar
         lsabarnum = abarnum
         lsabar = abar

         
c  .428431667757438 0.2146141197858406E+00011 0.1644676177994669E+00015   

         xcheck=1.0d0-(xnut+xprot+xalfa+xh)
         if(abs(xcheck)<1E-5) then
            xconsistent=1
         else
c if(abs(xcheck)<1E-2) then
            cpinf=1.0/0.0
c     Use fix for special case of somewhat bad xcheck but cp=Inf.  Assume rest of variables ok no matter how bad xcheck is
            if((1).AND.(abs(cp).eq.abs(cpinf))) then
               xconsistent=1
               cp = 0.2829658353841933E+00008
            else
               xconsistent=0
            end if
         end if

         cpinf=1.0/0.0
c     Use fix for special case of somewhat bad xcheck but cp=Inf.  Assume rest of variables ok no matter how bad cp is
         if((1).AND.(abs(cp).eq.abs(cpinf))) then
            xconsistent=1
            cp = 0.2829658353841933E+00008
         end if



         negative1=-1.0
         mynan=sqrt(negative1)
         
c         if(muhat.eq.mynan) then
c            muhatcheck=0
c     Still want to be able to output muhat even if didn't converge
c            muhat=-1
c         else
c            muhatcheck=1
c         end if


c         if((cv.ge.0.0).AND.(cp.ge.0.0)) then
c            cvcpcheck=1
c         else
c            cvcpcheck=0
c         end if


c         write(*,*) didconverge,consistent,positivepressure,positiveenergy,yeconsistent,xconsistent

c         if((didconverge.eq.1).AND.(yeconsistent.eq.0)) then
c            write (*,*) 'death'
c         end if


         if(didconverge.eq.1) then
            numconverged=numconverged+1
            didconverge_row(index)=didconverge
            dse_ls_row(index)=dse
            dpe_ls_row(index)=dpe
            dsp_ls_row(index)=dsp

c DEBUG: 1.OR.
            if(
c     1     1.OR.
     1           (
     1           (consistent.eq.1)
     1           .AND.(positivepressure.eq.1)
     1           .AND.(positiveenergy.eq.1)
     1           .AND.(yeconsistent.eq.1)
     1           .AND.(xconsistent.eq.1)
c     1           .AND.(cvcpcheck.eq.1)
c     1           .AND.(muhatcheck.eq.1)
c     1           .AND. (0.AND.(usediffinputs.ne.2 .OR. limitedrange.ne.1))
     1           .AND. (usediffinputs.ne.2 .OR. limitedrange.ne.1)
     1           .AND. (usediffinputs.ne.3 .OR. limitedrange.ne.1)
     1           .AND. (usediffinputs.ne.4 .OR. limitedrange.ne.1)
     1           )
     3           ) then
               goodconverge=1
               numgoodconverged=numgoodconverged+1
            else if( (usediffinputs.eq.2.OR.usediffinputs.eq.3.OR.usediffinputs.eq.4) .AND. limitedrange.eq.1) then
               goodconverge=0
c     Don't report how other things behave if this critical problem occurs
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     if usediffinputs==2 or 3, 4(kinda) then setting only A,Z to boundary values, and when detect fixing to boundary, enforce that solution was googconverge==0
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Using didconverge single value to indicate success of limited range, so can still use _row version to indicate approach
               didconverge_row(index)=-500
            else if(yeconsistent.ne.1) then
               goodconverge=0
c     Don't report how other things behave if this critical problem occurs
               didconverge_row(index)=-100
            else if(xconsistent.ne.1) then
               goodconverge=0
c     Don't report how other things behave if this critical problem occurs
               didconverge_row(index)=-200
c            else if(cvcpcheck.ne.1) then
c            else if(muhatcheck.ne.1) then
c               goodconverge=0
cc     Don't report how other things behave if this critical problem occurs
c               didconverge_row(index)=-300
            else
               goodconverge=0
               if(consistent.eq.1) then
                  didconverge_row(index)=-1
                  if((positivepressure.eq.1).AND.(positiveenergy.eq.1)) then
                     didconverge_row(index)=didconverge_row(index)-1
                  else if((positivepressure.eq.0).AND.(positiveenergy.eq.1)) then
                     didconverge_row(index)=didconverge_row(index)-2
                  else if((positivepressure.eq.1).AND.(positiveenergy.eq.0)) then
                     didconverge_row(index)=didconverge_row(index)-3
                  else if((positivepressure.eq.0).AND.(positiveenergy.eq.0)) then
                     didconverge_row(index)=didconverge_row(index)-4
                  else
                     write(*,*) "Why here1?"
                     write(*,*) consistent,positivepressure,positiveenergy,yeconsistent,xconsistent
                  end if
               else
                  didconverge_row(index)=-10
                  if((positivepressure.eq.1).AND.(positiveenergy.eq.1)) then
                     didconverge_row(index)=didconverge_row(index)-1
                  else if((positivepressure.eq.0).AND.(positiveenergy.eq.1)) then
                     didconverge_row(index)=didconverge_row(index)-2
                  else if((positivepressure.eq.1).AND.(positiveenergy.eq.0)) then
                     didconverge_row(index)=didconverge_row(index)-3
                  else if((positivepressure.eq.0).AND.(positiveenergy.eq.0)) then
                     didconverge_row(index)=didconverge_row(index)-4
                  else
                     write(*,*) "Why here2?"
                     write(*,*) consistent,positivepressure,positiveenergy,yeconsistent,xconsistent
                  end if
               end if
            end if
         else
c  Did not converge at all
            goodconverge=0
            didconverge_row(index)=0
            dse_ls_row(index)=1.0
            dpe_ls_row(index)=1.0
            dsp_ls_row(index)=1.0
         end if

         if(didconverge_row(index).ge.2) then
            write(*,*) "Why here3?"
            write(*,*) didconverge,consistent,positivepressure,positiveenergy,yeconsistent,xconsistent
         end if




         return
         end


















cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Convert from ??_cgs output to HELM or TIMMES EOS output (with perhaps some extra new things)
c     Go from single global to global array
c
c
c     Should NEVER have anything but ???_cgs terms on the right, since HELM/TIMMES modify _row terms directly when having called store_row().  Can't rely on single globals from HELM/TIMMES
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine store_row_fromcgsnuc2helmeos(loci) 
      implicit none

c     Passed quantity
c     Local index
      integer loci

c  contains array globals
      include 'eosparms.f'
      include 'vector_eos.dek'
c  contains array globals
      include 'extra_vector_sneos.dek'
c  contains single globals (only for those HELM/TIMMES-like quantities that LSEOS also computes)
      include 'vector_eos.single.dek'
c  contains single globals
      include 'vector_sneos.dek'
c  contains single global: xnut, and similar things
      include 'eos_m4c.commononly.inc'




c     First copy over pure electron quantities that LSEOS doesn't recompute by changed by the call to any_electron()
c     We need this since electron stored in lsindex (jhi_eos+1) and now storing back to
c     normal location <=jhi_eos
      call store_row(loci,whicheleeos,1)


c     Rest of quantities are new quantities or newly computed quantities from LSEOS directly or from something done in store_row_fromlseos2cgs()





c     LSEOS modified these from original choice that determined Y_e
      abar_row(loci) = lsabar
      abarnum_row(loci) = lsabarnum
      zbar_row(loci) = lszbar
      

c..store this row
        ptot_row(loci)   = ptot_cgs
        dpt_row(loci)    = dpdt_cgs
        dpd_row(loci)    = dpdd_cgs
c       LSEOS doesn't compute the below 2
c        if(1) then
c           dpa_row(loci)    = dpresda
c           dpz_row(loci)    = dpresdz
c        end if

c       NEW LSEOS quantity
        dpy_row(loci)     = dpdy_cgs

        etot_row(loci)   = etot_cgs
        det_row(loci)    = dedt_cgs
        ded_row(loci)    = dedd_cgs
c       LSEOS doesn't compute the below 2
c        if(1) then
c           dea_row(loci)    = denerda   
c           dez_row(loci)    = denerdz
c        end if

        dey_row(loci)     = dedy_cgs

        stot_row(loci)   = stot_cgs
        dst_row(loci)    = dsdt_cgs
        dsd_row(loci)    = dsdd_cgs
c       LSEOS doesn't compute the below 2
c        if(1) then
c           dsa_row(loci)    = dentrda      
c           dsz_row(loci)    = dentrdz
c        end if

        dsy_row(loci)     = dsdy_cgs


c   NEW LSEOS quantities
        ftot_row(loci)   = ftot_cgs
        dft_row(loci)    = dfdt_cgs
        dfd_row(loci)    = dfdd_cgs
c       LSEOS doesn't compute the below 2 dfda or dfdz
c        dfa_row(loci)    = dfda
c        dfz_row(loci)    = dfdz

        dfy_row(loci)     = dfdy_cgs


        prad_row(loci)   = prad_cgs
        erad_row(loci)   = erad_cgs
        srad_row(loci)   = srad_cgs
c       NEW LSEOS quantity
        frad_row(loci)   = frad_cgs
        

c HELM VERSIONS:
c       LSEOS doesn't compute the below 4
c        if(1) then
c           pion_row(loci)   = pion
c           eion_row(loci)   = eion
c           sion_row(loci)   = sion 
c           xni_row(loci)    = xni
c        end if

c LSEOS uses bulk rather than "ion"
        pion_row(loci)   = pbulk_cgs
        eion_row(loci)   = ebulk_cgs
        sion_row(loci)   = sbulk_cgs
c       NEW LSEOS quantity
        fion_row(loci)   = fbulk_cgs



c     Note that for HELM EOS there is no separate "positron" part, just all called "ele" and other part is 0
c     Note that pele_cgs contains both electrons and positrons since comes from "epress" that contained both
c       So should set ppos_row=0 here so no duplication of positron term

        pele_row(loci)   = pele_cgs
c  LSEOS doesn't produce below, so just assume doesn't change from globally set quantities
c        ppos_row(loci)   = ppos !new
c        ppos_row(loci)   = 0
c        dpept_row(loci)  = dpepdt
c        dpepd_row(loci)  = dpepdd
c        dpepa_row(loci)  = dpepda  
c        dpepz_row(loci)  = dpepdz

        eele_row(loci)   = eele_cgs
c  LSEOS doesn't produce below, so just assume doesn't change from globally set quantities
c        epos_row(loci)   = epos !new
c        epos_row(loci)   = 0
c        deept_row(loci)  = deepdt
c        deepd_row(loci)  = deepdd
c        deepa_row(loci)  = deepda   
c        deepz_row(loci)  = deepdz

        sele_row(loci)   = sele_cgs
c  LSEOS doesn't produce below, so just assume doesn't change from globally set quantities
c        spos_row(loci)   = spos !new
c        spos_row(loci)   = 0
c        dsept_row(loci)  = dsepdt 
c        dsepd_row(loci)  = dsepdd 
c        dsepa_row(loci)  = dsepda        
c        dsepz_row(loci)  = dsepdz


c       NEW LSEOS quantity
        fele_row(loci)   = fele_cgs

c  LSEOS doesn't produce below, so just assume doesn't change from globally set quantities
c        if(whicheleeos.eq.2) then
c           xnem_row(loci)   = xne
c        else if(whicheleeos.eq.0) then
c           xnem_row(loci)   = xnem
c        end if

c     Unlike presure, internal energy, and entropy, xnefer and xnpfer are not combined, so can keep separate
c     Same for etaele and etapos
c        xne_row(loci)    = xnefer
c        dxnet_row(loci)  = dxnedt
c        dxned_row(loci)  = dxnedd
c        dxnea_row(loci)  = dxneda
c        dxnez_row(loci)  = dxnedz
c        xnp_row(loci)    = xnpfer !new (OK to keep as separate!)

c  LSEOS doesn't produce below, so just assume doesn't change from globally set quantities
c     Since HELM/TIMMES/KAZ set _row(loci) in special way, don't want to overwrite that
c        etaele_row(loci) = etaele
c        detat_row(loci)  = detadt
c        detad_row(loci)  = detadd
c        detaa_row(loci)  = detada
c        detaz_row(loci)  = detadz
c        etapos_row(loci) = etapos !new (OK to keep as separate!)


c  LSEOS doesn't use below, so just assume doesn't change from globally set quantities
c        pcou_row(loci)   = pcoul
c        ecou_row(loci)   = ecoul
c        scou_row(loci)   = scoul 
c        plasg_row(loci)  = plasg
c     LSEOS absorbed coulomb terms directly into nuclei terms, so must turn off "cou" terms
        pcou_row(loci)   = 0.0
        ecou_row(loci)   = 0.0
        scou_row(loci)   = 0.0



c     Originally computed versions of these is overwritten
c     Ok even if not called ???_cgs, just use same variables since overwritten
        dse_row(loci)    = dse
        dpe_row(loci)    = dpe
        dsp_row(loci)    = dsp

        cv_row(loci)     = cv
        cp_row(loci)     = cp
        gam1_row(loci)   = gam1
        gam2_row(loci)   = gam2
        gam3_row(loci)   = gam3
        cs_row(loci)     = sound


c     Store species quantities
        call store_row_lsspecies(loci)

        




      return
      end
















c     Initialize TIMMES (and HELM) structures.
c     Note that TIMMES sets some things that HELM does not set
c     such as ppos, epos, etc. for positrons
c     But since rely universal variables, need to uniformly initialize
      subroutine init_row()
      implicit none



      include 'const.dek'

      include 'vector_eos.dek'
      include 'vector_eos.single.dek'
c     Caution:
      include 'kazeos.dek' ! includes a very similar name etae for etaele



      prad     = 0.0d0
      dpraddd  = 0.0d0
      dpraddt  = 0.0d0
      dpradda  = 0.0d0
      dpraddz  = 0.0d0

      erad     = 0.0d0
      deraddd  = 0.0d0
      deraddt  = 0.0d0
      deradda  = 0.0d0
      deraddz  = 0.0d0

      srad     = 0.0d0
      dsraddd  = 0.0d0
      dsraddt  = 0.0d0
      dsradda  = 0.0d0
      dsraddz  = 0.0d0

      xni      = 0.0d0
      dxnidd   = 0.0d0
      dxnidt   = 0.0d0
      dxnida   = 0.0d0
      dxnidz   = 0.0d0

      pion     = 0.0d0
      dpiondd  = 0.0d0
      dpiondt  = 0.0d0
      dpionda  = 0.0d0
      dpiondz  = 0.0d0

      eion     = 0.0d0
      deiondd  = 0.0d0
      deiondt  = 0.0d0
      deionda  = 0.0d0
      deiondz  = 0.0d0

      sion     = 0.0d0
      dsiondd  = 0.0d0
      dsiondt  = 0.0d0
      dsionda  = 0.0d0
      dsiondz  = 0.0d0

      xne      = 0.0d0
      dxnedd   = 0.0d0
      dxnedt   = 0.0d0
      dxneda   = 0.0d0
      dxnedz   = 0.0d0

      etaele   = 0.0d0
      detadd   = 0.0d0
      detadt   = 0.0d0
      detada   = 0.0d0
      detadz   = 0.0d0
      etapos   = 0.0d0

      xnefer    = 0.0d0
      dxneferdd = 0.0d0
      dxneferdt = 0.0d0
      dxneferda = 0.0d0
      dxneferdz = 0.0d0

      xnpfer    = 0.0d0
      dxnpferdd = 0.0d0
      dxnpferdt = 0.0d0
      dxnpferda = 0.0d0
      dxnpferdz = 0.0d0

      pele     = 0.0d0
      dpeledd  = 0.0d0
      dpeledt  = 0.0d0
      dpeleda  = 0.0d0
      dpeledz  = 0.0d0

      eele     = 0.0d0
      deeledd  = 0.0d0
      deeledt  = 0.0d0
      deeleda  = 0.0d0
      deeledz  = 0.0d0

      sele     = 0.0d0
      dseledd  = 0.0d0
      dseledt  = 0.0d0
      dseleda  = 0.0d0
      dseledz  = 0.0d0

      ppos     = 0.0d0
      dpposdd  = 0.0d0
      dpposdt  = 0.0d0
      dpposda  = 0.0d0
      dpeledz  = 0.0d0

      epos     = 0.0d0
      deposdd  = 0.0d0
      deposdt  = 0.0d0
      deposda  = 0.0d0
      deeledz  = 0.0d0

      spos     = 0.0d0
      dsposdd  = 0.0d0
      dsposdt  = 0.0d0
      dsposda  = 0.0d0
      dseledz  = 0.0d0

      pep      = 0.0d0
      dpepdd   = 0.0d0
      dpepdt   = 0.0d0
      dpepda   = 0.0d0
      dpepdz   = 0.0d0

      eep      = 0.0d0
      deepdd   = 0.0d0
      deepdt   = 0.0d0
      deepda   = 0.0d0
      deepdz   = 0.0d0

      sep      = 0.0d0
      dsepdd   = 0.0d0
      dsepdt   = 0.0d0
      dsepda   = 0.0d0
      dsepdz   = 0.0d0

c     JCM: For HELM, add-in eosfxt.f TIMMES-like code for eip
c     assumes potmult .eq. 0
      eip      = 0.0d0
      deipdd   = 0.0d0
      deipdt   = 0.0d0
      deipda   = 0.0d0
      deipdz   = 0.0d0

      sip      = 0.0d0
      dsipdd   = 0.0d0
      dsipdt   = 0.0d0
      dsipda   = 0.0d0
      dsipdz   = 0.0d0

      pcoul    = 0.0d0
      dpcouldd = 0.0d0
      dpcouldt = 0.0d0
      dpcoulda = 0.0d0
      dpcouldz = 0.0d0

      ecoul    = 0.0d0
      decouldd = 0.0d0
      decouldt = 0.0d0
      decoulda = 0.0d0
      decouldz = 0.0d0

      scoul    = 0.0d0
      dscouldd = 0.0d0
      dscouldt = 0.0d0
      dscoulda = 0.0d0
      dscouldz = 0.0d0



c     Initialize Kaz-like quantities
          
c     Needed for normal output of eos.dat
          dyedtthin = 0.0d0
          dyedt = 0.0d0
          graddotrhouyenonthermal=0.0d0
          graddotrhouye=0.0d0
          Qphoton=0.0d0
          Qm=0.0d0
          Nm=0.0d0
          Ynu = 0.0d0
          Ynuthermal = 0.0d0
          Rufgraddotrhouye=0.0d0
          RufQm=0.0d0
          RufNm=0.0d0
          lambdatot=0.0d0
          Enuglobal=0.0d0
          Enueglobal=0.0d0
          Enuebarglobal=0.0d0
          Tdifftot=0.0d0
          Tthermaltot=0.0d0

          lambdaintot=0.0d0

          tauphotonohcm = 0.0d0
          tauphotonabsohcm = 0.0d0

          unue0=0.0d0
          unuebar0=0.0d0
          unumu0 = 0.0d0

          nnue0=0.0d0
          nnuebar0=0.0d0
          nnumu0 = 0.0d0

          qtautelohcm = 0.0d0
          qtautnueohcm = 0.0d0
          qtautnuebarohcm = 0.0d0

          qtauaelohcm = 0.0d0
          qtauanueohcm = 0.0d0
          qtauanuebarohcm = 0.0d0

          qtautmuohcm = 0.0d0
          qtauamuohcm = 0.0d0
          qtauttauohcm = 0.0d0
          qtauatauohcm = 0.0d0

          ntautelohcm = 0.0d0
          ntautnueohcm = 0.0d0
          ntautnuebarohcm = 0.0d0

          ntauaelohcm = 0.0d0
          ntauanueohcm = 0.0d0
          ntauanuebarohcm = 0.0d0

          ntautmuohcm = 0.0d0
          ntauamuohcm = 0.0d0
          ntauttauohcm = 0.0d0
          ntauatauohcm = 0.0d0

c     Needed for normal output of eosother.dat
          xnuc = 0.0d0
          npratiofree = 0.0d0
          npratiobound = 0.0d0
          npratioheav = 0.0d0
          npratiototal = 0.0d0
          yetot = 0.0d0
          yeheav = 0.0d0

c     Don't zero out ls quantities that are already set as desired
c          etapls=0.0d0
c          etanls=0.0d0
          etap=0.0d0
          etan=0.0d0
          etanu=0.0d0

          nptotal =0.0d0
          nntotal =0.0d0
          npfree =0.0d0
          nnfree =0.0d0
          npbound =0.0d0
          nnbound =0.0d0
          npheav =0.0d0
          nnheav =0.0d0


          p_nu = 0.0d0
          rho_nu = 0.0d0
          s_nu = 0.0d0

          qtausohcm = 0.0d0
          Qmel = 0.0d0
          Qmmu = 0.0d0
          Qmtau = 0.0d0
          qminusel = 0.0d0
          qminusmu = 0.0d0
          qminustau = 0.0d0

          ntausohcm = 0.0d0
          Nmel = 0.0d0
          Nmmu = 0.0d0
          Nmtau = 0.0d0
          nminusel = 0.0d0
          nminusmu = 0.0d0
          nminustau = 0.0d0

          thermalye = 0.0d0
          yefree = 0.0d0
          abarbound=0.0d0
          yebound = 0.0d0
          gammapeglobal = 0.0d0
          gammapnuglobal = 0.0d0
          gammapenuglobal = 0.0d0
          gammanglobal = 0.0d0
          gammaneglobal = 0.0d0
          gammannuglobal = 0.0d0
          gammap2nglobal = 0.0d0
          gamman2pglobal = 0.0d0



      return
      end









c Used by jon_helmstandard.f (HELM EOS) and jon_eosfxt.f (TIMMES EOS) for storing single globals used to store later here into row vectors
c "store this row" in those 2 files has been combined into 1 subroutine
c Note that some additional things needed to interface with LSEOS that were in eosfxt.f and added here to use with HELM EOS also
c     Go from single global to global array
      subroutine store_row(loci,whichnonnucleareos,storenewspecies)
      implicit none

c     Passed quantity
c     Local index
      integer loci
      integer whichnonnucleareos,storenewspecies
      

c..bring in the lattimer-swesty data structures
c      include 'eos_m4c.commononly.inc'

      include 'const.dek'

      include 'eosparms.f'
      include 'vector_eos.dek'
      include 'vector_eos.single.dek'
c      include 'vector_eos.extra4kaz.dek'
c     Caution:
      include 'kazeos.dek' ! includes a very similar name etae for etaele
      

c..store this row
        ptot_row(loci)   = pres
        dpt_row(loci)    = dpresdt
        dpd_row(loci)    = dpresdd
        dpa_row(loci)    = dpresda   
        dpz_row(loci)    = dpresdz

        etot_row(loci)   = ener
        det_row(loci)    = denerdt
        ded_row(loci)    = denerdd
        dea_row(loci)    = denerda   
        dez_row(loci)    = denerdz

        stot_row(loci)   = entr
        dst_row(loci)    = dentrdt
        dsd_row(loci)    = dentrdd
        dsa_row(loci)    = dentrda   
        dsz_row(loci)    = dentrdz

        prad_row(loci)   = prad
        if(1) then
           dpradt_row(loci) = dpraddt
           dpradd_row(loci) = dpraddd
           dprada_row(loci) = dpradda
           dpradz_row(loci) = dpraddz
        end if


        erad_row(loci)   = erad
        if(1) then
           deradt_row(loci) = deraddt
           deradd_row(loci) = deraddd
           derada_row(loci) = deradda
           deradz_row(loci) = deraddz
        end if

        srad_row(loci)   = srad 
        if(1) then
           dsradt_row(loci) = dsraddt
           dsradd_row(loci) = dsraddd
           dsrada_row(loci) = dsradda
           dsradz_row(loci) = dsraddz
        end if

        pion_row(loci)   = pion
        eion_row(loci)   = eion
        sion_row(loci)   = sion 

        xni_row(loci)    = xni

        pele_row(loci)   = pele+ppos
c        ppos_row(loci)   = ppos !new
        ppos_row(loci)   = 0.0
        dpept_row(loci)  = dpepdt
        dpepd_row(loci)  = dpepdd
        dpepa_row(loci)  = dpepda  
        dpepz_row(loci)  = dpepdz

        eele_row(loci)   = eele+epos
c        epos_row(loci)   = epos !new
        epos_row(loci)   = 0.0
        deept_row(loci)  = deepdt
        deepd_row(loci)  = deepdd
        deepa_row(loci)  = deepda   
        deepz_row(loci)  = deepdz

        sele_row(loci)   = sele+spos
c        spos_row(loci)   = spos !new
        spos_row(loci)   = 0.0
        dsept_row(loci)  = dsepdt 
        dsepd_row(loci)  = dsepdd 
        dsepa_row(loci)  = dsepda        
        dsepz_row(loci)  = dsepdz

        if(1) then

           if(whicheleeos.eq.2) then
              xnem_row(loci)   = xne
           else if(whicheleeos.eq.0) then
              xnem_row(loci)   = xnem
           end if
        else
           xnem_row(loci)   = xnem
        end if

        xne_row(loci)    = xnefer

        if(whicheleeos.eq.2) then
           dxnet_row(loci)  = dxnedt
           dxned_row(loci)  = dxnedd
           dxnea_row(loci)  = dxneda
           dxnez_row(loci)  = dxnedz
        else if(whicheleeos.eq.0) then
           dxnet_row(loci)  = dxneferdt + dxnpferdt
           dxned_row(loci)  = dxneferdd + dxnpferdd
           dxnea_row(loci)  = dxneferda + dxnpferda
           dxnez_row(loci)  = dxneferdz + dxnpferdz
        end if

        xnp_row(loci)    = xnpfer !new

        if(1) then
c     Below was only set in TIMMES EOS
           zeff_row(loci)   = zeff
        end if




c        Store species information
c       \eta's
c     GODMARK: Careful to choose correct etae vs. etaele below
c        etaele_row(loci) = etae
        etaele_row(loci) = etaele

c     Don't set etap_row and etan_row here!
c        etap_row(loci) = etapls
c        etan_row(loci) = etanls
c        etanu_row(loci) = etanu
c        write(*,*) 'storerow',loci,etaele,etapls,etanls,etanu

        detat_row(loci)  = detadt
        detad_row(loci)  = detadd
        detaa_row(loci)  = detada
        detaz_row(loci)  = detadz
        etapos_row(loci) = etapos !new

        if(1) then
           if(whicheleeos.eq.2) then
              eip_row(loci)    = eip
              sip_row(loci)    = sip
           else if(whicheleeos.eq.0) then
c     Not set since not assigned!           
           end if
        end if

        pcou_row(loci)   = pcoul
        ecou_row(loci)   = ecoul
        scou_row(loci)   = scoul 
        plasg_row(loci)  = plasg

        dse_row(loci)    = dse
        dpe_row(loci)    = dpe
        dsp_row(loci)    = dsp

        cv_row(loci)     = cv
        cp_row(loci)     = cp
        gam1_row(loci)   = gam1
        gam2_row(loci)   = gam2
        gam3_row(loci)   = gam3
        cs_row(loci)     = sound

c..for debugging
c        crap1_row(loci)   = etaele
c        dcrap1d_row(loci) = detadd
c        dcrap1t_row(loci) = detadt
c        dcrap1a_row(loci) = detada
c        dcrap1z_row(loci) = detadz









ccccccccccccccccccccccccccccccccccccccc
c
c     Store Kaz-like quantities
c
ccccccccccccccccccccccccccccccccccccccc

          
c     Needed for normal output of eos.dat
          dyedtthin_row(loci)=          dyedtthin
          dyedt_row(loci)=          dyedt
          graddotrhouyenonthermal_row(loci) = graddotrhouyenonthermal
          graddotrhouye_row(loci) = graddotrhouye
          Qphoton_row(loci)=          Qphoton
          Qm_row(loci)=          Qm
          Nm_row(loci)=          Nm
          Ynu_row(loci) = Ynu
          Ynuthermal_row(loci) = Ynuthermal
          Rufgraddotrhouye_row(loci) = Rufgraddotrhouye
          RufQm_row(loci)=          RufQm
          RufNm_row(loci)=          RufNm

          lambdatot_row(loci)=          lambdatot

c          write(*,*) 'lambdainside1',lambdatot,lambdatot_row(loci)

          Enuglobal_row(loci)=          Enuglobal
          Enueglobal_row(loci)=          Enueglobal
          Enuebarglobal_row(loci)=          Enuebarglobal
          Tdifftot_row(loci)=          Tdifftot
          Tthermaltot_row(loci)=          Tthermaltot


          lambdaintot_row(loci) =           lambdaintot

          tauphotonohcm_row(loci) =           tauphotonohcm
          tauphotonabsohcm_row(loci) =           tauphotonabsohcm

          unue0_row(loci) =           unue0
          unuebar0_row(loci) =           unuebar0
          unumu0_row(loci) =           unumu0

          nnue0_row(loci) =           nnue0
          nnuebar0_row(loci) =           nnuebar0
          nnumu0_row(loci) =           nnumu0

          qtautelohcm_row(loci) =           qtautelohcm
          qtautnueohcm_row(loci) =           qtautnueohcm
          qtautnuebarohcm_row(loci) =           qtautnuebarohcm

          qtauaelohcm_row(loci) =           qtauaelohcm
          qtauanueohcm_row(loci) =           qtauanueohcm
          qtauanuebarohcm_row(loci) =           qtauanuebarohcm

          qtautmuohcm_row(loci) =           qtautmuohcm
          qtauamuohcm_row(loci) =           qtauamuohcm
          qtauttauohcm_row(loci) =           qtauttauohcm
          qtauatauohcm_row(loci) =           qtauatauohcm

          ntautelohcm_row(loci) =           ntautelohcm
          ntautnueohcm_row(loci) =           ntautnueohcm
          ntautnuebarohcm_row(loci) =           ntautnuebarohcm

          ntauaelohcm_row(loci) =           ntauaelohcm
          ntauanueohcm_row(loci) =           ntauanueohcm
          ntauanuebarohcm_row(loci) =           ntauanuebarohcm

          ntautmuohcm_row(loci) =           ntautmuohcm
          ntauamuohcm_row(loci) =           ntauamuohcm
          ntauttauohcm_row(loci) =           ntauttauohcm
          ntauatauohcm_row(loci) =           ntauatauohcm


c     Needed for normal output of eosother.dat
c     If storenewspecies.eq.0 then assume already computed and don't want to replace
          if(storenewspecies.eq.1) then
             xnuc_row(loci)=          xnuc
             npratiofree_row(loci)=          npratiofree
             npratiobound_row(loci)=npratiobound
             npratioheav_row(loci)=npratioheav
             npratiototal_row(loci)=npratiototal
             yetot_row(loci)=          yetot
             yefree_row(loci)=          yefree
             yebound_row(loci)=          yebound
             yeheav_row(loci)=          yeheav

             abarbound_row(loci) = abarbound

             nptotal_row(loci)=nptotal
             nntotal_row(loci)=nntotal
             npfree_row(loci)=npfree
             nnfree_row(loci)=nnfree
             npbound_row(loci)=npbound
             nnbound_row(loci)=nnbound
             npheav_row(loci)=npheav
             nnheav_row(loci)=nnheav


          end if

          p_nu_row(loci)=          p_nu
c     Put \rho and s into HELM form
          rho_nu_row(loci)=          rho_nu/den_row(loci)
          s_nu_row(loci)=          (kerg*s_nu)/den_row(loci)

          qtausohcm_row(loci)=          qtausohcm
          Qmel_row(loci)=          Qmel
          Qmmu_row(loci)=          Qmmu
          Qmtau_row(loci)=          Qmtau
          qminusel_row(loci)=          qminusel
          qminusmu_row(loci)=          qminusmu
          qminustau_row(loci)=          qminustau

          ntausohcm_row(loci)=          ntausohcm
          Nmel_row(loci)=          Nmel
          Nmmu_row(loci)=          Nmmu
          Nmtau_row(loci)=          Nmtau
          nminusel_row(loci)=          nminusel
          nminusmu_row(loci)=          nminusmu
          nminustau_row(loci)=          nminustau

          thermalye_row(loci)=          thermalye
          gammapeglobal_row(loci)=          gammapeglobal
          gammapnuglobal_row(loci)=          gammapnuglobal
          gammapenuglobal_row(loci)=          gammapenuglobal
          gammanglobal_row(loci)=          gammanglobal
          gammaneglobal_row(loci)=          gammaneglobal
          gammannuglobal_row(loci)=          gammannuglobal
          gammap2nglobal_row(loci)=          gammap2nglobal
          gamman2pglobal_row(loci)=          gamman2pglobal





      return
      end




c inverse of store_row() above.  Stores from array to singles -- used before file writing
c     Used to store what Kaz generated and put into arrays to be placed back into singles
      subroutine storeback_row(loci)
      implicit none

c     Passed quantity
c     Local index
      integer loci

c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
      
      include 'const.dek'

      include 'eosparms.f'
      include 'vector_eos.dek'
      include 'vector_eos.single.dek'
c      include 'vector_eos.extra4kaz.dek'
c     Caution:
      include 'kazeos.dek' ! includes a very similar name etae for etaele







c..   store this row back into individual variables
      pres =       ptot_row(loci)
      dpresdt =       dpt_row(loci)
      dpresdd =       dpd_row(loci)
      dpresda    =       dpa_row(loci)
      dpresdz =       dpz_row(loci)

      ener =       etot_row(loci)
      denerdt =       det_row(loci)
      denerdd =       ded_row(loci)
      denerda    =       dea_row(loci)
      denerdz =       dez_row(loci)

      entr =       stot_row(loci)
      dentrdt =       dst_row(loci)
      dentrdd =       dsd_row(loci)
      dentrda    =       dsa_row(loci)
      dentrdz =       dsz_row(loci)

      prad =       prad_row(loci)
      if(1) then
         dpraddt =          dpradt_row(loci)
         dpraddd =          dpradd_row(loci)
         dpradda =          dprada_row(loci)
         dpraddz =          dpradz_row(loci)
      end if


      erad =       erad_row(loci)
      if(1) then
         deraddt =          deradt_row(loci)
         deraddd =          deradd_row(loci)
         deradda =          derada_row(loci)
         deraddz =          deradz_row(loci)
      end if

      srad  =       srad_row(loci)
      if(1) then
         dsraddt =          dsradt_row(loci)
         dsraddd =          dsradd_row(loci)
         dsradda =          dsrada_row(loci)
         dsraddz =          dsradz_row(loci)
      end if

      pion =       pion_row(loci)
      eion =       eion_row(loci)
      sion  =       sion_row(loci)

      xni =       xni_row(loci)

      pele =       pele_row(loci)
      ppos =       ppos_row(loci)
      dpepdt =       dpept_row(loci)
      dpepdd =       dpepd_row(loci)
      dpepda   =       dpepa_row(loci)
      dpepdz =       dpepz_row(loci)

      eele =       eele_row(loci)
      epos =       eele_row(loci)
      deepdt =       deept_row(loci)
      deepdd =       deepd_row(loci)
      deepda    =       deepa_row(loci)
      deepdz =       deepz_row(loci)

      sele =       sele_row(loci)
      spos =       sele_row(loci)
      dsepdt  =       dsept_row(loci)
      dsepdd  =       dsepd_row(loci)
      dsepda         =       dsepa_row(loci)
      dsepdz =       dsepz_row(loci)

      if(1) then

         if(whicheleeos.eq.2) then
            xne =             xnem_row(loci)
         else if(whicheleeos.eq.0) then
            xnem =             xnem_row(loci)
         end if
      else
         xnem =          xnem_row(loci)
      end if

      xnefer =       xne_row(loci)

      if(whicheleeos.eq.2) then
         dxnedt =          dxnet_row(loci)
         dxnedd =          dxned_row(loci)
         dxneda =          dxnea_row(loci)
         dxnedz =          dxnez_row(loci)
      else if(whicheleeos.eq.0) then
         dxneferdt  =          dxnet_row(loci)
         dxnpferdt =          dxnet_row(loci)

         dxneferdd =          dxned_row(loci)
         dxnpferdd =          dxned_row(loci)

         dxneferda =          dxnea_row(loci)
         dxnpferda =          dxnea_row(loci)

         dxneferdz =          dxnez_row(loci)
         dxnpferdz =          dxnez_row(loci)
      end if

      xnpfer =       xnp_row(loci) !new

      if(1) then
c     Below was only set in TIMMES EOS
         zeff =          zeff_row(loci)
      end if


c     \eta's
c     GODMARK: WHICH ONE OF BELOW:
      etaele =       etaele_row(loci)

c     write(*,*) 'storerow',loci,etaele,etapls,etanls,etanu

      detadt =       detat_row(loci)
      detadd =       detad_row(loci)
      detada =       detaa_row(loci)
      detadz =       detaz_row(loci)
      etapos =       etapos_row(loci)!new

      if(1) then
         if(whicheleeos.eq.2) then
            eip =             eip_row(loci)
            sip =             sip_row(loci)
         else if(whicheleeos.eq.0) then
c     Not set since not assigned!           
         end if
      end if

      pcoul =       pcou_row(loci)
      ecoul =       ecou_row(loci)
      scoul  =       scou_row(loci)
      plasg =       plasg_row(loci)

      dse =       dse_row(loci)
      dpe =       dpe_row(loci)
      dsp =       dsp_row(loci)

      cv =       cv_row(loci)
      cp =       cp_row(loci)
      gam1 =       gam1_row(loci)
      gam2 =       gam2_row(loci)
      gam3 =       gam3_row(loci)
      sound =       cs_row(loci)






      etae=etaele_row(loci)
      etaele=etae
      etap=etap_row(loci)
      etapls=etap
      etan=etan_row(loci)
      etanls=etan
      etanu=etanu_row(loci)
c      write(*,*) 'storebackrow',loci,etae,etap,etan,etanu


cccccccccccccccc
c
c     Translate stored true Kaz variables in _row form into normal single form
c     That is, for output we use single name vars but need that value for
c     each ???_row(loci)
c     Other Kaz quantities have HELM/TIMMES analogue so came from a different
c     ???_row name already translated
         
c     Needed for normal output of eos.dat
      dyedtthin  =           dyedtthin_row(loci)
      dyedt  =           dyedt_row(loci)
      graddotrhouyenonthermal = graddotrhouyenonthermal_row(loci)
      graddotrhouye = graddotrhouye_row(loci)
      Qphoton=Qphoton_row(loci)
      Qm  =           Qm_row(loci)
      Nm  =           Nm_row(loci)
      Ynu = Ynu_row(loci)
      Ynuthermal = Ynuthermal_row(loci)
      Rufgraddotrhouye = Rufgraddotrhouye_row(loci)
      RufQm  =           RufQm_row(loci)
      RufNm  =           RufNm_row(loci)

      lambdatot  =              lambdatot_row(loci)

c      write(*,*) 'lambdainside2',lambdatot,lambdatot_row(loci)

      Enuglobal  =              Enuglobal_row(loci)
      Enueglobal  =              Enueglobal_row(loci)
      Enuebarglobal  =              Enuebarglobal_row(loci)
      Tdifftot  =              Tdifftot_row(loci)
      Tthermaltot  =           Tthermaltot_row(loci)


      lambdaintot=lambdaintot_row(loci)
      
      tauphotonohcm=tauphotonohcm_row(loci)
      tauphotonabsohcm=tauphotonabsohcm_row(loci)
      
      unue0=unue0_row(loci)
      unuebar0=unuebar0_row(loci)
      unumu0=unumu0_row(loci)
      
      nnue0=nnue0_row(loci)
      nnuebar0=nnuebar0_row(loci)
      nnumu0=nnumu0_row(loci)



      qtautelohcm  =           qtautelohcm_row(loci)
      qtautnueohcm  =           qtautnueohcm_row(loci)
      qtautnuebarohcm  =           qtautnuebarohcm_row(loci)

      qtauaelohcm  =           qtauaelohcm_row(loci)
      qtauanueohcm  =           qtauanueohcm_row(loci)
      qtauanuebarohcm  =           qtauanuebarohcm_row(loci)

      qtautmuohcm  =           qtautmuohcm_row(loci)
      qtauamuohcm  =           qtauamuohcm_row(loci)
      qtauttauohcm  =           qtauttauohcm_row(loci)
      qtauatauohcm  =           qtauatauohcm_row(loci)

      ntautelohcm  =           ntautelohcm_row(loci)
      ntautnueohcm  =           ntautnueohcm_row(loci)
      ntautnuebarohcm  =           ntautnuebarohcm_row(loci)

      ntauaelohcm  =           ntauaelohcm_row(loci)
      ntauanueohcm  =           ntauanueohcm_row(loci)
      ntauanuebarohcm  =           ntauanuebarohcm_row(loci)

      ntautmuohcm  =           ntautmuohcm_row(loci)
      ntauamuohcm  =           ntauamuohcm_row(loci)
      ntauttauohcm  =           ntauttauohcm_row(loci)
      ntauatauohcm  =           ntauatauohcm_row(loci)

c     Needed for normal output of eosother.dat
      xnuc  =           xnuc_row(loci)
      npratiofree  =           npratiofree_row(loci)
      npratiobound = npratiobound_row(loci)
      npratioheav = npratioheav_row(loci)
      npratiototal = npratiototal_row(loci)
      yetot  =           yetot_row(loci)
      yefree  =           yefree_row(loci)
      yebound  =           yebound_row(loci)
      yeheav  =           yeheav_row(loci)

      abarbound = abarbound_row(loci)

      nptotal = nptotal_row(loci)
      nntotal = nntotal_row(loci)
      npfree = npfree_row(loci)
      nnfree = nnfree_row(loci)
      npbound = npbound_row(loci)
      nnbound = nnbound_row(loci)
      npheav = npheav_row(loci)
      nnheav = nnheav_row(loci)


      p_nu  =           p_nu_row(loci)
c     Convert back to Kaz form for \rho and s
c     Notice the change of energy to per unit energy
      rho_nu  =           rho_nu_row(loci)*den_row(loci)
      s_nu  =           (s_nu_row(loci)/kerg)*den_row(loci)

      qtausohcm  =           qtausohcm_row(loci)
      Qmel  =           Qmel_row(loci)
      Qmmu  =           Qmmu_row(loci)
      Qmtau  =           Qmtau_row(loci)
      qminusel  =           qminusel_row(loci)
      qminusmu  =           qminusmu_row(loci)
      qminustau  =           qminustau_row(loci)

      ntausohcm  =           ntausohcm_row(loci)
      Nmel  =           Nmel_row(loci)
      Nmmu  =           Nmmu_row(loci)
      Nmtau  =           Nmtau_row(loci)
      nminusel  =           nminusel_row(loci)
      nminusmu  =           nminusmu_row(loci)
      nminustau  =           nminustau_row(loci)

      thermalye  =           thermalye_row(loci)
      gammapeglobal  =           gammapeglobal_row(loci)
      gammapnuglobal  =           gammapnuglobal_row(loci)
      gammapenuglobal  =           gammapenuglobal_row(loci)
      gammanglobal  =           gammanglobal_row(loci)
      gammaneglobal  =           gammaneglobal_row(loci)
      gammannuglobal  =           gammannuglobal_row(loci)
      gammap2nglobal  =           gammap2nglobal_row(loci)
      gamman2pglobal  =           gamman2pglobal_row(loci)


      return
      end







c     Some extra Kaz-like calculations to be done by HELM/TIMMES
      subroutine helmtimmes_extracalculation(loci,whichnonnucleareos)
      implicit none

c     Passed quantity
c     Local index
      integer loci
      integer whichnonnucleareos
      
c     Need to add-in rest-mass so consistent with LSEOS and Kaz EOS
      real*8 rhoblocal,nb,neleposi
      real*8 yelocal,yesumlocal,eele_rest
      real*8 betatemp
      real*8 xnuccalc !function
      real*8 nbtotal ! tempvar


      include 'const.dek'

      include 'eosparms.f'
      include 'vector_eos.dek'
      include 'vector_eos.single.dek'
c      include 'vector_eos.extra4kaz.dek'
c     Caution:
      include 'kazeos.dek' ! includes a very similar name etae for etaele
      





c     This wasn't computed by HELM/TIMMES, so compute it now
      xnuc=xnuccalc(den_row(loci),temp_row(loci))
      yetot = zbar_row(loci)/abarnum_row(loci)
      
      kazxheav = xnuc
      abarbound = abarbound_row(loci)

      call computeyefreebound(den_row(loci),temp_row(loci),yetot,xnuc,yefree,yebound)
      yeheav = yebound
      
      rhoblocal = den_row(loci)
      nbtotal = (rhoblocal/mb)
c     Below computes some generally-true things needed for Kaz-EOS
      call eos2kazeos_species(nbtotal)
   



cccccccccccccccc
c     Translate some quantities into right physics format
c     None of these rest-mass corrections affect the thermodynamic consistency calculations
c     But the sound speed is affected through the total specific energy

c     JCM: etaele has rest-mass subtracted out
c     JCM: etaepos has double rest-mass added in
c     JCM: etaepos = -etae-2/betatemp
c     JCM: When outputting, make (LSEOS/Kaz)-like (no subtraction)
      betatemp = kerg*temp_row(loci)/(mecc)



c     electron rest-mass NOT negligible!
      yelocal = zbar_row(loci)/abarnum_row(loci)
      rhoblocal = den_row(loci)
      nb = (rhoblocal/mb)
c     yesumlocal=(xne_row(loci)+xnp_row(loci))/nb
      eele_rest=mecc/mb*yelocal
c     eele_rest=mecc/mb*yesumlocal


c     Add rest-mass of electrons back into chemical potential
c     GODMARK: NOTE THAT the RHS should only have non-loci things unless inputs (den,temp,abar,abarnum,zbar)
      etaele = etaele+1.0/betatemp
      etapos = etapos+1.0/betatemp

      eele = eele + epos + eele_rest
      epos = 0.0
      ener = ener + eele_rest

      deepdz = deepdz + mecc*avo/(abar_row(loci)*yelocal)





      return
      end











c     Some extra Kaz-like calculations to be done by all EOSs
c     Requires already globally set: yetot,yefree,yebound,yeheav,xnuc,kazxheav
      subroutine eos2kazeos_species(nbtotal)
      implicit none

c     Passed quantity
      real*8 nbtotal
      
      include 'const.dek'

c      include 'eosparms.f'
c      include 'vector_eos.dek'
c      include 'vector_eos.single.dek'
c      include 'vector_eos.extra4kaz.dek'
c     Caution:
      include 'kazeos.dek' ! includes a very similar name etae for etaele
      

      npratiototal = (1.0-yetot)/yetot
      npratiofree=(1.0-yefree)/yefree
      npratiobound = (1.0-yebound)/yebound
      
      etanu=0.0
      
      nptotal = nbtotal*yetot
      nntotal = nptotal*npratiototal
      
      npfree  = nbtotal*xnuc*yefree
      nnfree  = npfree*npratiofree
      
      npbound = nptotal-npfree
      nnbound = nntotal-nnfree
      if(npbound<0.0) then
         npbound=0.0
      end if
      if(nnbound<0.0) then
         nnbound=0.0
      end if
      
      npratioheav = (1.0-yeheav)/yeheav
      npheav = nbtotal*kazxheav*yeheav
      nnheav = npheav*npratioheav
      

      return
      end












c     Go from Kaz EOS single globals to HELM/TIMMES internals
c     Go from single global to another single global
      subroutine store_row_kazeos2helminternals()
      implicit none


      include 'const.dek'
c     Bring in HELM single globals
      include 'vector_eos.single.dek'
c     Bring in Kaz single globals
      include 'kazeos.dek'



      pres = p_tot
      dpresdt=0.0
      dpresdd=0.0
      dpresda=0.0
      dpresdz=0.0

      ener=u_tot/rhob
      denerdt=0.0
      denerdd=0.0
      denerda=0.0
      denerdz=0.0

      entr = s_tot/rhob
      dentrdt=0.0
      dentrdd=0.0
      dentrda=0.0    
      dentrdz=0.0

      prad=p_photon
      dpraddt=0.0
      dpraddd=0.0
      dpradda=0.0
      dpraddz=0.0


      erad=rho_photon/rhob
      deraddt=0.0
      deraddd=0.0
      deradda=0.0
      deraddz=0.0

      srad=s_photon/rhob
      dsraddt=0.0
      dsraddd=0.0
      dsradda=0.0
      dsraddz=0.0

      pion=p_N
      eion=u_N/rhob
      sion=s_N/rhob

      xni=0.0

      pele = p_eleposi
      ppos=0.0
      dpepdt=0.0
      dpepdd=0.0
      dpepda=0.0
      dpepdz=0.0

      eele=rho_eleposi/rhob
      epos=0.0
      deepdt=0.0
      deepdd=0.0
      deepda=0.0
      deepdz=0.0

      sele=s_eleposi/rhob
      spos=0.0
      dsepdt=0.0
      dsepdd=0.0
      dsepda=0.0    
      dsepdz=0.0

      xne=0.0
      xnem=0.0
      xnefer=0.0

      dxnedt=0.0
      dxnedd=0.0
      dxneda=0.0
      dxnedz=0.0
      dxneferdt=0.0
      dxnpferdt=0.0
      dxneferdd=0.0
      dxnpferdd=0.0
      dxneferda=0.0
      dxnpferda=0.0
      dxneferdz=0.0
      dxnpferdz=0.0

      xnpfer=0.0

      zeff=0.0


      etaele=etae
      detadt=0.0
      detadd=0.0
      detada=0.0
      detadz=0.0
      etapos=-etae

      eip=0.0
      sip=0.0

      pcoul=0.0
      ecoul=0.0
      scoul =0.0
      plasg=0.0

      dse=0.0
      dpe=0.0
      dsp=0.0

      cv=0.0
      cp=0.0
      gam1=0.0
      gam2=0.0
      gam3=0.0
      sound=0.0

      detadd=0.0
      detadt=0.0
      detada=0.0
      detadz=0.0


c     Kaz returns EOS with electron rest-mass in electron internal energy
c     So no need to translate to another physics form
c     Actually, when HELM store_row() is called next, adds in rest-mass related terms, so need to avoid for Kaz EOS




      return
      end













ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Obtain electron/radiation EOS
c     Nuclear EOSs call this for only electron/radiation terms.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc

c     tin = temperature in MeV
c     yein = Y_e = Y_p = n_p / n_b
c     dnsin = density in number of baryons/fm^3
      subroutine any_electron(tin,yein,dnsin) 
      implicit none

c..declare the pass
      double precision tin,yein,dnsin
      double precision tintrue,dnsintrue

      

      include 'const.dek'
      include 'vector_sneos.dek'

c      Apparently only actually called once with known tin,yein,dnsin for each LSEOS
c      write(*,*) 'Called any_electron()',tin,yein,dnsin


      if(usediffinputs.eq.1 .OR. usediffinputs.eq.2 .OR. usediffinputs.eq.3 .OR. usediffinputs.eq.4) then
         tintrue = temp_cgs*k2mev
         dnsintrue = den_cgs*(avo*fm3)
      else
c     Otherwise no change, but below would be correct, but already set correctly
c         tin = temp_cgs_lseos*k2mev
c         dnsin = den_cgs_lseos*(avo*fm3)
         tintrue=tin
         dnsintrue=dnsin
      end if




c Just a check -- only setup for lsindex=1 and so only 1 state
c      Old one is in jon_ls_2p7.f
c      call any_electron_old(tintrue,yeintrue,dnsintrue)
      call any_electron_new(tintrue,yein,dnsintrue)
      


      return
      end











c JCM: previously in jon_ls_2p7.f
      subroutine any_electron_new(tin,yein,dnsin) 
      implicit none

      include 'eos_m4c.commononly.inc'
c JCM: below used to store electron EOS properties needed by LSEOS to generate final solution
      include 'el_eos.inc'
c  JCM:
      include 'const.dek'

      include 'eosparms.f'
      include 'vector_eos.dek'
      include 'vector_sneos.dek'

c..declare the pass
      double precision tin,yein,dnsin

c     Local variables
      double precision abarnum,abar,zbar,yelocal


c..local variables
c      double precision ytot1,zbarxx,abar,abarnum,zbar,yelocal

c JCM:
c..for unit conversions
c      double precision fm,fm3,avo,kerg,kev,k2mev,mev2erg,
c     1                 clight,light2,me,mecc
c      parameter        (fm      = 1.0d-13,
c     1                  fm3     = fm*fm*fm,
c JCM:
c     2                  avo     = 6.0221367d23,
c     3                  kerg    = 1.380658d-16,
c     4                  kev     = 8.617385d-5,
c     5                  k2mev   = kev * 1.0d-6,
c     6                  mev2erg = ev2erg * 1.0d6,
c     7                  clight  = 2.99792458d10,
c     8                  light2  = clight*clight,
c     9                  me      = 9.1093897d-28,
c     &                  mecc    = me*light2)

      integer loci
      integer saved_jlo_eos,saved_jhi_eos
      integer dostore,dostorespecies



c..debug formats
 111  format(1x,1p3e24.16)
 112  format(1x,1p6e14.6)





c     If whicheleeos==1, then input to any_electron() function is all that's needed
      if(whicheleeos.eq.1) then
c..   the original e+e- eos in ls
         call el_eos(tin,yein,dnsin)
c     Store quantities in memory pointed to by el_eos.inc used by LSOES
c     totally done since already filled el_eos.inc
         return
      else


c     HELM and TIMMES EOSs use global variables already processed.  Save the ones that need to be changed before running then, so that they can be restored after that run
         saved_jlo_eos=jlo_eos
         saved_jhi_eos=jhi_eos



c     Note that "index" specifies what component of ???_row(?) to use while computing HELM or TIMMES EOS
c     Since that is a global quantity, must choose same one as set that got here
c         loci=index
         loci=lsindex
         

         if(loci.ne.(jhi_eos+1)) then
            write(*,*) 'loci PROBLEM',loci
         end if
         
c..   set the input vector in cgs units
c     This should have changed nothing about temp_row and den_row
         temp_row(loci) = tin / k2mev
         den_row(loci)  = dnsin / (avo*fm3)

         call compute_nuclear_azbar(xnut, xprot, xalfa, xh, a, x,
     1        abarnum,abar,zbar, yelocal)

         if(abs(yelocal-yein)>100.0*yetolerance) then
            write (*,*) 'input ye to any_electron
     1           is not consistent with X_i',yein,yelocal
         end if
         
         
c         write(*,*) 'abar,abarnum,zbar',abar,zbar
c
c        Keep same for now
         abar_row(loci) = abar
         abarnum_row(loci) = abarnum
         zbar_row(loci) = zbar


c     "Loop" is over LSEOS extra storage for HELM/TIMMES EOS array
         jlo_eos = loci
         jhi_eos = loci


cccccccccccccccccccccccccccccccc
c     
c     
c     NOW RUN CHOSEN EOS


c         write(*,*) 'eleEOSgets',temp_row(loci),den_row(loci)
c     1        ,abar_row(loci),abarnum_row(loci),zbar_row(loci)



c     Whether to perform extra store in _row() beyond normal store into jlo_eos and jhi_eos
c     Automatically stores once into loci already, so don't need to do again
         dostore=0
         dostorespecies=0
         call full_nonnuclear_eos(whicheleeos, loci, dostore,dostorespecies)



c     NOW restore saved versions
c     Now any EOS processes called before reaching here are preserved
         jlo_eos=saved_jlo_eos
         jhi_eos=saved_jhi_eos





cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Convert units and construct those things needed by LSEOS routines
c     If using prior internal electron EOS within LSEOS, then already in right format for LSEOS
c
c JCM: These ??_row values (HELM/TIMMES output format) are stored into native LSEOS electron variables
c      These native LSEOS electron variables are later stored as ??_cgs and ??_nuc values
c      These ??_cgs values are then later stored as ??_row values
c      GODMARK: These _row's should be same as final _row's since electron/photon don't change
cccccccccccccccccccccccccccccccccccccccccccccccccccc



c..   set the output vector
         nsube  = xnem_row(loci) * fm3
         neplus = xnp_row(loci) * fm3
c     Already added rest-mass into HELM/TIMMES myself
c         musube = (etaele_row(loci) * kerg*temp_row(loci) + mecc)/mev2erg
         musube = (etaele_row(loci) * kerg*temp_row(loci))/mev2erg

c..   electron thermodynamics
         epress = (pele_row(loci) + ppos_row(loci)) * fm3/mev2erg
c     Already added rest-mass into HELM/TIMMES myself
c         eu     = (eele_row(loci) + epos_row(loci) + mecc*avo*yelocal)
c     1        / (avo * mev2erg)
         eu     = (eele_row(loci) + epos_row(loci))
     1        / (avo * mev2erg)
         es     = (sele_row(loci) + spos_row(loci)) / (kerg*avo)
         fsube  = eu - tin * es

c..   photon thermodynamics
         ppress = prad_row(loci) * fm3/mev2erg
         pu     = erad_row(loci) / (avo * mev2erg)
         ps     = srad_row(loci) / (kerg*avo)
         pf     = pu - tin * ps


c..   chemical potential derivaties
         demudt = (detat_row(loci)*temp_row(loci) + etaele_row(loci)) * kerg
     1        / (mev2erg*k2mev)
         demudn = detad_row(loci)*kerg*temp_row(loci) / mev2erg / (avo*fm3)
         demudy = (detaz_row(loci) - detaa_row(loci)/yelocal)
     1        * abar_row(loci)*kerg*temp_row(loci) * yelocal / mev2erg

c..   electron derivatives
         depdt  = dpept_row(loci) * fm3/mev2erg / k2mev
         depdn  = dpepd_row(loci) * fm3/mev2erg / (avo*fm3)
         depdy  = (dpepz_row(loci) - dpepa_row(loci)/yelocal)
     1        * abar_row(loci) * yelocal * fm3/mev2erg

         deudt  = deept_row(loci) / (avo * mev2erg) / k2mev
         deudn  = deepd_row(loci) /(avo * mev2erg) / (avo*fm3)
c     Already added rest-mass into HELM/TIMMES myself
c         deudy  = ((deepz_row(loci) - deepa_row(loci)/yelocal)
c     1        *abar_row(loci)*yelocal  + mecc*avo)
c     2        /(avo * mev2erg)
         deudy  = ((deepz_row(loci) - deepa_row(loci)/yelocal)
     1        *abar_row(loci)*yelocal)
     2        /(avo * mev2erg)

         desdt  = dsept_row(loci) / (kerg*avo) / k2mev
         desdn  = dsepd_row(loci) / (kerg*avo) / (avo*fm3)
         desdy  = (dsepz_row(loci) -dsepa_row(loci)/yelocal)*abar_row(loci)*yelocal 
     1        / (kerg*avo) 

c..   photon derivatives
         dppdt  = dpradt_row(loci) * fm3/mev2erg / k2mev
         dppdn  = dpradd_row(loci) * fm3/mev2erg / (avo*fm3)
         dppdy  = (dpradz_row(loci) - dprada_row(loci)/yelocal)
     1        *abar_row(loci)*yelocal * fm3/mev2erg 

         dpudt  = deradt_row(loci) / (avo * mev2erg) / k2mev
         dpudn  = deradd_row(loci) / (avo * mev2erg) / (avo*fm3)
         dpudy  = (deradz_row(loci) - derada_row(loci)/yelocal)
     1        *abar_row(loci)*yelocal/ (avo * mev2erg)

         dpsdt  = dsradt_row(loci) / (kerg*avo) / k2mev
         dpsdn  = dsradd_row(loci) / (kerg*avo) / (avo*fm3)
         dpsdy  = (dsradz_row(loci) - dsrada_row(loci)/yelocal)
     1        *abar_row(loci)*yelocal / (kerg*avo) 



      end if
c     if here then whicheleeos.ne.1




      return
      end


      





cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     lsbox() has been changed quite a bit from original wrapper
c
c     
c     lsbox() inputs cgs values of Y_e, T[K], and \rho_b[g/cc]
c     lsbox() does:
c     1) Reads LS guess file
c     2) Gets LS solution into LSEOS-type quantities into global singular values stored in eos_m4c.commononly.inc, el_eos.inc and vector_sneos.dek
c     3) Finally converts LSEOS-type quantities to cgs
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       subroutine lsbox(iunitlocal, yelocal, tklocal, rhoblocal, didconverge)
       implicit none
       save

c..does all the stuff needed for calling the lattimer-swesty eos


c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
      include 'el_eos.inc'
c  JCM:
c  contains single globals
      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
      include 'const.dek'

c     JCM: passed parameters
      integer iunitlocal
      double precision yelocal,tklocal,rhoblocal

c     JCM: passed parameter to return to tell if converged or not
      integer didconverge

c  JCM:
c..for unit conversions
c      double precision fm,fm3,avo,kerg,kev,k2mev,mev2erg,clight,light2
c      parameter        (fm      = 1.0d-13,
c     1                  fm3     = fm*fm*fm,
c     2                  avo     = 6.0221367d23,
c     3                  kerg    = 1.380658d-16,
c     4                  kev     = 8.617385d-5,
c     5                  k2mev   = kev * 1.0d-6,
c     6                  mev2erg = 1.602177d-6,
c     7                  clight  = 2.99792458d10,
c     8                  light2  = clight*clight)



c..local variables
      integer          iflg,eflg,fflg,sf
      double precision ipvar(4),told,xpr,ppr,yp,etanguess,etapguess

      integer          ifirst
      double precision temp_old
      data             ifirst/0/, temp_old/1.0d0/

      double precision tmin,rhomin,tmax,rhomax, ypmin, ypmax

ccccccccccccccccccccccccccccccccccccc
c JCM:
c     Set up input state (only allows cgs input)
c     
c     Note that above is true Ye, while rest-mass version would be using abar_row(index)
c
cccccccccccccccccccccccccccccccccccccc

      iunit    = iunitlocal
      ye_inp   = yelocal
      temp_cgs = tklocal
      den_cgs  = rhoblocal

      

ccccccccccccccccccccccccccccccccccccc

c     Fix rhob and T near region where LSEOS is bad, but continue using normal rho,T for electron EOS
c     6.412,   8.87604
      tmin = 10**9
      tmax = 10**14
      rhomin = 10**6.4
      rhomax = 10**15
      ypmin=0.04
      ypmax=0.9


c     Further restrict and change rho,T limits
      call tableminmaxfixes(tmin,tmax,rhomin,rhomax, ypmin, ypmax)

c
cccccccccccccccccccccccccccccccccccccc
      



c..let the user know loading the initial guess file and the 
c..fermi integral tables might take a few seconds

      if (ifirst .eq. 0) then
       write(6,*)
       write(6,*) 'loading initial guess and fermi integral tables'
       call loadguess
       call loadmx
       write(6,*) 'done loading initial guess and fermi integral tables'
       ifirst = 1
      end if

ccccccccccccccccccccccccccccccccccc
c
c..convert the input units
c
ccccccccccccccccccccccccccccccccccc

      if (iunit .eq. 1) then
c..cgs units being input
       temp_nuc = temp_cgs * k2mev
       temp_nuc_lseos = temp_cgs_lseos * k2mev
c  JCM:
c  baryon number density in (number/fm^3):
       den_nuc  = den_cgs * avo * fm3
       den_nuc_lseos  = den_cgs_lseos * avo * fm3
      else if (iunit .eq. 2) then
c..nuclear units being input
       temp_cgs = temp_nuc / k2mev
       temp_cgs_lseos = temp_nuc_lseos / k2mev
       den_cgs  = den_nuc/(avo*fm3)
       den_cgs_lseos  = den_nuc_lseos/(avo*fm3)
      end if



c..make a good guess at the proton and neutron chemical potentials,
c..and exterior proton fraction


c      write(*,*) "Calling lsguess"
      call lsguess(den_cgs_lseos,temp_cgs_lseos,ye_inp,etanguess,etapguess,yp)


c..set the lattimer-swesty eos input vector
      ipvar(1) = temp_nuc_lseos
      ipvar(2) = 0.155
      ipvar(3) = etapguess
      ipvar(4) = etanguess
      ppr      = den_nuc_lseos*yp
      fflg     = 0
      iflg     = 1
      told     = 1.0d0


cccccccccccccccccccccccccccccccccccc
c
c..call the lattimer-swesty eos
c
cccccccccccccccccccccccccccccccccccc

c      write(*,*) "Begin Calling inveos"
      call inveos(ipvar,told,ye_inp,den_nuc_lseos,iflg,
     1            eflg,fflg,sf,xpr,ppr)
c      write(*,*) "End Calling inveos"



c..failure
      if (sf .ne. 1) then
         if(0) then
            write(6,*) 
            write(6,*) ' eos failed at ye temp den (cgs):'
            write(6,111) ye_inp,temp_cgs_lseos,den_cgs_lseos
            write(6,111) ye_inp,temp_nuc_lseos,den_nuc_lseos
            write(6,111) etapguess,etanguess,yp
 111        format(1x,1p3e24.16)
            stop 'ls did not converge'
         else
c     Don't report, just go back and tell calling function if we converged
            didconverge=0
c            return
         end if
      else
c     otherwise assume converged
         didconverge=1
      end if

c      write(*,*) didconverge


c     Convert from LSEOS-like values to cgs values
      call store_row_fromlseos2cgs(didconverge)





      return
      end












c     Convert from LSEOS output to _cgs and _nuc output (with some extra new things)
      subroutine store_row_fromlseos2cgs(didconverge)
      implicit none


c     JCM: passed parameter to return to tell if converged or not
      integer didconverge

c     local variables
      double precision ytot1,zbarxx

c JCM:
c Used to check if cp is Infinity (i.e. must be same type?)
      double precision cpinf

c JCM:
c Used to offset total specific energy
      double precision yelocal
      double precision lsoffset

c      double precision npratiobound

      include 'kazeos.dek'
      include 'eosparms.f'
c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
      include 'el_eos.inc'
c  JCM:
c  contains single globals like sound, dse, etc.
      include 'vector_eos.single.dek'
      include 'vector_eos.dek'
c Below contains _cgs and _nuc globals
      include 'vector_sneos.dek'
      include 'const.dek'

      double precision bpressnew,bsnew,bunew
      double precision nbtotal ! tempvar
      real*8 rhoblocal
      real*8 mutotxnuc0,ataunse
      real*8 aorig,pnogas





      if(didconverge.eq.1) then


c JCM: note that the below means that:
c xprot  = 1 n_p / n_b
c xneut  = 1 n_n / n_b
c xalpha = 4 n_\alpha / n_b
c xheav  = a n_h / n_b  where A_h=a and Z_h=x*a




ccccccccccccccccc
c     No longer use this method to fix non-NSE
      if(1.eq.0) then
c     If taunsefix activated, assume only need to correct a, keeping x same
         if(taunsefix.eq.1) then
            aorig=a
c     Override use of nuclear EOS if NSE timescale is too long
c     Purpose is to fit presupernovae model better
c     Note use of true density and temperature, not _lseos versions
            call computemutotfit_xnuc0(den_cgs,temp_cgs,mutotxnuc0)
c     ataunse = (xh+xalfa)/((xh+xalfa)/mutotxnuc0-xalfa/4.0)
c     Ensure that abar~ataunse for heavy components
            ataunse = xh/((xh+xalfa)/mutotxnuc0-xalfa/4.0)
            if(ataunse.lt.0.5) then
               ataunse = 1.0
               write(*,*) 'TYPE1'
            end if
            if(den_cgs.lt.rhotaunsefix/1.5 .OR. temp_cgs.lt.temptaunsefix/1.5) then
               a=ataunse
               write(*,*) 'TYPE2'
            else
c     average in between
               write(*,*) 'TYPE3'
               a=0.5*(a+ataunse)
            end if
c     Keep x (i.e. Y_e for heavy nuclei) the same


c     Must modify pion
c     pnogas=ptot-bpress
c     Doesn't account for change in Coulomb correction
c     bpress=bpress*aorig/a
c     ptot = pnogas + bpress

         end if
      end if








         call compute_nuclear_azbar(xnut, xprot, xalfa, xh, a, x,
     1        abarnum,abar,zbar,yelocal)
         lszbar = zbar
         lsabarnum = abarnum
         lsabar = abar


c     DEBUG:
c      write(*,*) 'INLSEOS',xnut,xprot,xalfa,xh,a,x,'LSBAR',abarnum,abar,zbar,yelocal


c     abarnum,abar,zbar here are stored for later storing back to ??_row versions
c      lsabar=abar
c      lsabarnum=abar
c      lszbar=zbar



c     LS EOS OFFSET
      if(whichnucleareos.eq.1) then
c     Modify the total for energy/gram
c         lsoffset=QEMEV*yelocal+16.0*yelocal*2.0
c         lsoffset=QEMEV*yelocal+16.0
c         lsoffset=16.0
c         lsoffset=9.14
c     Incorrect to enforce offset, must leave so that utot< possible
c         lsoffset=0.0

c     See top of jon_lsbox.f.  One must ensure total mass-energy density is actually correct -- can't be arbitrarily offset
c     For LSEOS this means adding 8.07131747535936MeV to the nuclear term
         lsoffset=8.07131747535936
         bu = bu + lsoffset
         utot  = utot + lsoffset

c     Most negative specific energy gets is about 10MeV/baryon
c     GODMARK: This still leaves poor connection to ideal gas
c         lsoffset=20.0*yelocal

c     Was using below:
c         bu = bu + lsoffset
c         utot  = utot + lsoffset
c         if(utot.lt.0.0) then
c            write(*,*) bu,utot,lsoffset
c         end if
      end if


c     SHEN EOS OFFSET
      if(whichnucleareos.eq.3) then
c     Modify the total for energy/gram
c         lsoffset=QEMEV*yelocal*8.0
c     Incorrect to enforce offset, must leave so that utot< possible
c         lsoffset=9.14
c         lsoffset = 0.0
c Absolute minimum Shen ebulk is -9.132329941 according to table (happens at large Y_e and low temperatures when nuclei are very bound)
c         lsoffset=0         

c     See top of jon_lsbox.f.  One must ensure total mass-energy density is actually correct -- can't be arbitrarily offset
c     For LSEOS this means adding 8.07131747535936MeV to the nuclear term
         lsoffset=8.07131747535936
         bu = bu + lsoffset
         utot  = utot + lsoffset

c
c     Was using below
c         bu = bu + lsoffset
c         utot  = utot + lsoffset
c         if(utot.lt.0.0) then
c            write(*,*) bu,utot,lsoffset
c         end if
      end if



cccccccccccccccccccccccccccccccccccccccccccccc
c
c     Try for enforce good connection between LSEOS and normal EOS by avoiding very negative pressures in nuclear terms and very negative internal energies after already above accounted for 0-point energy
c
cccccccccccccccccccccccccccccccccccccccccccccc
c      if(usediffinputs.eq.1 .OR. usediffinputs.eq.2 .OR. usediffinputs.eq.3 .OR. usediffinputs.eq.3) then
c
c         if(
c     1        ((bpress.lt.0.0).AND.(abs(bpress)/abs(ptot).gt.0.3)).AND.(den_cgs.lt.1E6)
c     1        ) then
c     Then don't trust negative correction by nuclei
c            bpressnew = bpress*0.5
c            bpressnew = 0
c            ptot = ptot - bpress + bpressnew
c            bpress = bpressnew
c         end if
c
c     Below not needed as long as set lsoffset above
c         if(
c     1        ((bu.lt.0.0).AND.(abs(bu)/abs(utot).gt.0.3)).AND.(den_cgs.lt.1E6)
c     1        ) then
c     Then don't trust negative correction by nuclei
c            bunew=bu*0.5
c            bunew=0
c            utot = utot - bu + bunew
c            bu = bunew
c         end if
c      end if

c reduced form that was using:
c      if(usediffinputs.eq.1 .OR. usediffinputs.eq.2 .OR. usediffinputs.eq.3 .OR. usediffinputs.eq.3) then
c         if(
c     1        ((bpress.lt.0.0).AND.(abs(bpress)/abs(ptot).gt.0.3)).AND.(den_cgs.lt.1E6)
c     1        ) then
c            bpressnew = 0
c            ptot = ptot - bpress + bpressnew
c            bpress = bpressnew
c         end if
c      end if


c..photon, e+e-, bulk, and total pressures in mev/fm**3 and erg/cm**3

      prad_nuc  = ppress
      prad_cgs  = prad_nuc * mev2erg/fm3

      pele_nuc  = epress
      pele_cgs  = pele_nuc * mev2erg/fm3

      pbulk_nuc = bpress
      pbulk_cgs = pbulk_nuc * mev2erg/fm3

      ptot_nuc  = ptot
      ptot_cgs  = ptot_nuc * mev2erg/fm3

      dpdt_nuc  = dpdt
      dpdt_cgs  = dpdt_nuc * mev2erg/fm3 * k2mev
      dpdd_nuc  = dpdn
      dpdd_cgs  = dpdd_nuc * mev2erg * avo
      dpdy_nuc  = dpdy
      dpdy_cgs  = dpdy_nuc * mev2erg/fm3



c..photon, e+e-, bulk, and total entropies dimensionless and erg/g/k

      srad_nuc  = ps
      srad_cgs  = srad_nuc * (kerg*avo)

      sele_nuc  = es
      sele_cgs  = sele_nuc * (kerg*avo)

      sbulk_nuc = bs
      sbulk_cgs = sbulk_nuc * (kerg*avo)

      stot_nuc  = stot
      stot_cgs  = stot_nuc * (kerg*avo)
      dsdt_nuc  = dsdt
      dsdt_cgs  = dsdt_nuc * (kerg*avo) * k2mev
      dsdd_nuc  = dsdn
      dsdd_cgs  = dsdd_nuc * (kerg*avo) * avo * fm3
      dsdy_nuc  = dsdy
      dsdy_cgs  = dsdy_nuc * (kerg*avo)



c..photon, e+e-, bulk, and total internal energies in mev/baryon and erg/g

      erad_nuc  = pu
      erad_cgs  = erad_nuc * (avo * mev2erg)

      eele_nuc  = eu
      eele_cgs  = eele_nuc * (avo * mev2erg)

      ebulk_nuc = bu
      ebulk_cgs = ebulk_nuc * (avo * mev2erg)

      etot_nuc  = utot
      etot_cgs  = etot_nuc * (avo * mev2erg)
      dedt_nuc  = dudt
      dedt_cgs  = dedt_nuc * (avo * mev2erg) * k2mev
      dedd_nuc  = dudn
      dedd_cgs  = dedd_nuc * (avo * mev2erg) * avo * fm3
      dedy_nuc  = dudy
      dedy_cgs  = dedy_nuc * (avo * mev2erg)



c..photon, e+e-, bulk, and total helmholtz free energy in mev/baryon and erg/g

      frad_nuc  = pf
      frad_cgs  = frad_nuc * (avo * mev2erg)

      fele_nuc  = fsube
      fele_cgs  = fele_nuc * (avo * mev2erg)

      fbulk_nuc = bftot
      fbulk_cgs = fbulk_nuc * (avo * mev2erg)

      ftot_nuc  = ftot
      ftot_cgs  = ftot_nuc * (avo * mev2erg)
      dfdt_nuc  = -stot_nuc
      dfdt_cgs  = dfdt_nuc * (avo * mev2erg) * k2mev
      dfdd_nuc  = ptot_nuc/(den_nuc*den_nuc)
      dfdd_cgs  = dfdd_nuc * (avo * mev2erg) * avo * fm3
      dfdy_nuc  = ptot_nuc/(den_nuc*den_nuc)
      dfdy_cgs  = dfdy_nuc * (avo * mev2erg)




c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)
      zz    = ptot_cgs/den_cgs
      zzi   = den_cgs/ptot_cgs
      chit  = temp_cgs/ptot_cgs * dpdt_cgs
      chid  = dpdd_cgs * zzi
      cv    = dedt_cgs
      xx    = zz * chit/(temp_cgs * cv)
      gam3  = xx + 1.0d0
      gam1  = chit*xx + chid
      nabad = xx/gam1
      gam2  = 1.0d0/(1.0d0 - nabad)
      cp    = cv * gam1/chid
      zz    = 1.0d0 + (etot_cgs + light2)*zzi
      if(gam1>0.0 .AND. zz>0.0) then
         sound = clight * sqrt(gam1/zz)
      else
         sound = 0.1*clight
      end if

c..maxwell relations; each is zero if the consistency is perfect
      xx  = den_cgs * den_cgs
      dse = temp_cgs*dsdt_cgs/dedt_cgs - 1.0d0
      if(ptot_cgs>0.0) then
         dpe = (dedd_cgs*xx + temp_cgs*dpdt_cgs)/ptot_cgs - 1.0d0
      else
         dpe = 1E30
      end if
      dsp = -dsdd_cgs*xx/dpdt_cgs - 1.0d0

c      write(*,*) 'rho,tk,ye,zbar,abar',den_cgs,temp_cgs,zbar,abar
c      write(*,*) 'dse,dpe,dsp',dse,dpe,dsp



      if(1) then
c     JCM: For LSEOS, at \rho=1.644676...E14g/cc T=5E8-1E11, cp=Inf and cv~-Inf
c     And sound=Nan
c     If this happens, assume cp just large or phase transition
c     In any case reset cp to be large
         if(cp>1D49) then
            cp=1D49
            sound=clight*0.9999
         end if
         if(cp<-1D49) then
            cp=-1D49
            sound=clight*0.9999
         end if

         if(cv>1D49) then
            cv=1D49
            sound=clight*0.9999
         end if
         if(cv<-1D49) then
            cv=-1D49
            sound=clight*0.9999
         end if

         if(sound.ne.sound) then
            sound=clight*0.11
         end if


         if(cp-cp .ne. 0.0) then
            cp = 1D49
         end if

         if(cv.ne.cv) then
            cv = 1D49
         end if

         cpinf=1.0/0.0

         if(abs(cp).eq.abs(cpinf)) then
            cp = 1D49
         end if

         if(abs(cp).gt.huge(cp)) then
            cp = 1D49
         end if


c         write(*,*) 'mycp',cp
c         if(cp.eq.cp) then
c            write(*,*) 'nextcp',1.0/0.0
c         end if

c     End trying to fix cp,cv,sound
      end if








c         if((abs(dse).lt.1E-5)
c     2      .AND.(abs(dpe).lt.1E-5)
c     3      .AND.(abs(dsp).lt.1E-5)
c     4      .AND.(etot_cgs.lt.0.0)
c     5      .AND.(abs(zbar_row(jhi_eos+1)/abarnum_row(jhi_eos+1) - zbar_row(index)/abarnum_row(index))<1E-5)) then
            
c            if(bu<0) then
c     write(*,*) lsoffset,QEMEV,yelocal,lsoffset*(avo*mev2erg)
c               write(6,100) '\rho,T,Ye,unuc,offset',
c                write(6,100) den_row(index),temp_row(index),zbar_row(index)/abarnum_row(index),utot,bu,lsoffset
c            end if
c         end if







c     Here's where we set Kaz-like quantities for
c     xnuc, npratio, yetot, yefree, yebound
c     That is, the nuclear EOS is overriding whatever the electron EOS thought

ccccccccccccccccccccccccccccccccc
c     "Compute" a few more things




        xnuc = (xnut+xprot)
c        npratiofree = (xnut/mn)/(xprot/mp+1E-30)
c     Shen says never use mp or mn related to xnut or xprot
        npratiofree = (xnut)/(xprot+1E-30)
        yetot = zbar/abarnum
        npratiototal = (1.0-yetot)/yetot
        yefree = 1.0/(1.0+npratiofree)

c     Bottom has 1E-30 to avoid Inf
        npratiobound = (2.0*xalfa+a*(1.0-x)*xh)/
     1       (2.0*xalfa+a*x*xh+1E-30)

        yebound = 1.0/(1.0+npratiobound)

        if(xh>1E-20) then
           abarbound = 1.0/(xh/(a+1E-20) + xalfa/4.0)
           abarbound = abarbound*(xh+xalfa)
        else if(xalfa>1E-20) then
           abarbound = 1.0/(xalfa/4.0)
           abarbound = abarbound*(xh+xalfa)
        else
c     no bound particles
           abarbound = 1.0
        end if
        if(abarbound<1.0) then
           abarbound=1.0
        end if

        kazxheav = xh
        yeheav = x

        rhoblocal = den_cgs_lseos
        nbtotal = (rhoblocal/mb)
        
c     Below as in kazeos_compute.f
c     Below computes some generally-true things needed for Kaz-EOS
        call eos2kazeos_species(nbtotal)



 100  format(40(1pe26.15,' '))








      else
         call zero_cgsnuc_quantities()
      end if






      return
      end













cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Zero-out cgsnuc quantities
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zero_cgsnuc_quantities()
      implicit none

c     Passed quantity
c     Local index
      integer didconverge


c     local variables
      double precision ytot1,zbarxx

c JCM:
c Used to offset total specific energy
      double precision yelocal
      double precision lsoffset


      include 'const.dek'
c..bring in the lattimer-swesty data structures
      include 'eos_m4c.commononly.inc'
      include 'el_eos.inc'
c  JCM:
c  contains single globals like sound, dse, etc.
      include 'vector_eos.single.dek'
      include 'vector_eos.dek'
c Below contains _cgs and _nuc globals
      include 'vector_sneos.dek'




         prad_nuc  =1E-30
         prad_cgs  =1E-30

         pele_nuc  =1E-30
         pele_cgs  =1E-30

         pbulk_nuc =1E-30
         pbulk_cgs =1E-30

         ptot_nuc  =1E-30
         ptot_cgs  =1E-30

         dpdt_nuc  =1E-30
         dpdt_cgs  =1E-30
         dpdd_nuc  =1E-30
         dpdd_cgs  =1E-30
         dpdy_nuc  =1E-30
         dpdy_cgs  =1E-30



c..   photon, e+e-, bulk, and total entropies dimensionless and erg/g/k

         srad_nuc  =1E-30
         srad_cgs  =1E-30

         sele_nuc  =1E-30
         sele_cgs  =1E-30

         sbulk_nuc =1E-30
         sbulk_cgs =1E-30

         stot_nuc  =1E-30
         stot_cgs  =1E-30
         dsdt_nuc  =1E-30
         dsdt_cgs  =1E-30
         dsdd_nuc  =1E-30
         dsdd_cgs  =1E-30
         dsdy_nuc  =1E-30
         dsdy_cgs  =1E-30



c..   photon, e+e-, bulk, and total internal energies in mev/baryon and erg/g

         erad_nuc  =1E-30
         erad_cgs  =1E-30

         eele_nuc  =1E-30
         eele_cgs  =1E-30

         ebulk_nuc =1E-30
         ebulk_cgs =1E-30

         etot_nuc  =1E-30
         etot_cgs  =1E-30
         dedt_nuc  =1E-30
         dedt_cgs  =1E-30
         dedd_nuc  =1E-30
         dedd_cgs  =1E-30
         dedy_nuc  =1E-30
         dedy_cgs  =1E-30



c..   photon, e+e-, bulk, and total helmholtz free energy in mev/baryon and erg/g

         frad_nuc  =1E-30
         frad_cgs  =1E-30

         fele_nuc  =1E-30
         fele_cgs  =1E-30

         fbulk_nuc =1E-30
         fbulk_cgs =1E-30

         ftot_nuc  =1E-30
         ftot_cgs  =1E-30
         dfdt_nuc  =1E-30
         dfdt_cgs  =1E-30
         dfdd_nuc  =1E-30
         dfdd_cgs  =1E-30
         dfdy_nuc  =1E-30
         dfdy_cgs  =1E-30



c..   abar and zbar of the whole mixture
         ytot1  =1E-30
         if (a .ne. 0) ytot1 =1E-30
         zbarxx =1E-30
         abar   =1E-30
         zbar   =1E-30


c..   the temperature and density exponents (c&g 9.81 9.82) 
c..   the specific heat at constant volume (c&g 9.92)
c..   the third adiabatic exponent (c&g 9.93)
c..   the first adiabatic exponent (c&g 9.97) 
c..   the second adiabatic exponent (c&g 9.105)
c..   the specific heat at constant pressure (c&g 9.98) 
c..   and relativistic formula for the sound speed (c&g 14.29)
         zz    =1E-30
         zzi   =1E-30
         chit  =1E-30
         chid  =1E-30
         cv    =1E-30
         xx    =1E-30
         gam3  =1E-30
         gam1  =1E-30
         nabad =1E-30
         gam2  =1E-30
         cp    =1E-30
         zz    =1E-30
         sound =1E-30

c..   maxwell relations; each is zero if the consistency is perfect
         xx  =1E-30
         dse =1E-30
         dpe =1E-30
         dsp =1E-30


         return
         end








      subroutine loadguess
      implicit none
      save

c..loads a table of initial guesses for the lattimer-swesty eos

c..two common blocks for passing around the table 
      integer          nxdata, nydata, nzdata
      parameter        (nxdata=19, nydata=200, nzdata=85)

      double precision xdata(nxdata),ydata(nydata),zdata(nzdata)
      common /gridg/   xdata,ydata,zdata

      double precision pprevdata(nxdata,nydata,nzdata),
     &                 etanguessdata(nxdata,nydata,nzdata),
     &                 etapguessdata(nxdata,nydata,nzdata)
      common /gdata/   pprevdata,etanguessdata,etapguessdata


c..local variables
      integer          i,j,k
      double precision lo,step,lrho,tmev,ye,etanguess,etapguess,yp


c..ye from 0.05 to 0.5 in steps of 0.025
      lo    = 0.05d0
      step  = 0.025d0
      do i=1,nxdata
       xdata(i) = lo + float(i-1)*step
      enddo


c..temp in mev from 0.5 to 100 in steps of 0.5
      lo    = 0.5d0
      step  = 0.5d0
      do i=1,nydata
       ydata(i) = lo + float(i-1)*step
      enddo


c..log density from 7.0 to 15.4 in steps of 0.1
      lo    = 7.0d0
      step  = 0.1d0
      do i=1,nzdata
       zdata(i) = lo + float(i-1)*step
      enddo


c..open the file with the guesses for a 220 compressibility 
      open(unit=2, file='ls220_guesses.dat', status='old')

c      open(unit=3, file='ls220_g2.dat', status='unknown')

c      open(unit=3, file='ls220_g2.dat', status='old')


c..read the file
      do k=1,nzdata
       do j=1,nydata
        do i=1,nxdata
         read(2,*,end=22,err=23) lrho,tmev,ye,etanguess,etapguess,yp
         pprevdata(i,j,k)    = yp
         etanguessdata(i,j,k)= etanguess
         etapguessdata(i,j,k)= etapguess
        enddo
       enddo
      enddo


c..close up shop
 22   close(unit=2)

      return 

c..oops, bad read
 23   stop 'bad read in routine load'
      end







      subroutine lsguess(rho,temp,ye,etanlin,etaplin,yplin)
      implicit none
      save

c..interpolates the lattimer-swesty guesstimates between grid points

c..input is the density rho in g/cm**3, the temperature temp in kelvin, 
c..and the average nuclear charge to weight ratio ye.

c..output are initial guesses for the lattimer-swesty eos;
c..the neutron chemical potential etanlin, proton chemical potential 
c..etaplin, and the exterior proton fraction yplin.


c..declare the pass
      double precision rho,temp,ye,etanlin,etaplin,yplin


c..two common blocks for passing around the table 
      integer          nxdata, nydata, nzdata
      parameter        (nxdata=19, nydata=200, nzdata=85)

      double precision xdata(nxdata),ydata(nydata),zdata(nzdata)
      common /gridg/   xdata,ydata,zdata

      double precision pprevdata(nxdata,nydata,nzdata),
     &                 etanguessdata(nxdata,nydata,nzdata),
     &                 etapguessdata(nxdata,nydata,nzdata)
      common /gdata/   pprevdata,etanguessdata,etapguessdata


c..local variables
      integer          i,j,k
      double precision lrho,tmev,r,s,t


c..lower bounds and strides of the table
      double precision yemin,deltye,dyeinv,tmevmin,dtmev,dtmevinv,
     &                 lrhomin,dlrho,dlrhoinv
      parameter        (yemin   = 0.05d0,  
     &                 deltye   = 0.025d0,  
     &                 dyeinv   = 1.0d0/deltye,
     &                 tmevmin  = 0.5d0, 
     &                 dtmev    = 0.5d0,    
     &                 dtmevinv = 1.0d0/dtmev,
     &                 lrhomin  = 7.0d0, 
     &                 dlrho    = 0.1d0,    
     &                 dlrhoinv = 1.0d0/dlrho)


      include 'const.dek'

c..for unit conversions
c      double precision fm,fm3,avo,kerg,kev,k2mev,mev2erg,clight,light2
c      parameter        (fm      = 1.0d-13,
c     1                  fm3     = fm*fm*fm,
c     2                  avo     = 6.0221367d23,
c     3                  kerg    = 1.380658d-16,
c     4                  kev     = 8.617385d-5,
c     5                  k2mev   = kev * 1.0d-6,
c     6                  mev2erg = 1.602177d-6,
c     7                  clight  = 2.99792458d10,
c     8                  light2  = clight*clight)



c..work with the cgs density in log10 space
c..convert the cgs temperature to mev
      lrho= log10(rho)
      tmev = temp * k2mev


c..hash locate the table indices and bound them
      i = int((ye - yemin)*dyeinv) + 1
      i = max(1,min(i,nxdata-1))
      j = int((tmev-tmevmin)*dtmevinv) + 1
      j = max(1,min(j,nydata-1))
      k = int((lrho - lrhomin)*dlrhoinv) + 1
      k = max(1,min(k,nzdata-1))


c..normalized position within the cube
      r = (ye   - (i*deltye + (yemin   - deltye)))*dyeinv
      s = (tmev - (j*dtmev  + (tmevmin - dtmev)))*dtmevinv
      t = (lrho - (k*dlrho  + (lrhomin - dlrho)))*dlrhoinv


c..linear interpolations
c..proton chemical potential
      etaplin=  (1.-r)*(1.-s)*(1.-t)*etapguessdata(i,j,k)
     &           + r*(1.-s)*(1.-t)*etapguessdata(i+1,j,k)
     &           + r*s*(1.-t)*etapguessdata(i+1,j+1,k)
     &           + (1.-r)*s*(1.-t)*etapguessdata(i,j+1,k)
     &           + (1.-r)*(1.-s)*t*etapguessdata(i,j,k+1)
     &           + r*(1.-s)*t*etapguessdata(i+1,j,k+1)
     &           + r*s*t*etapguessdata(i+1,j+1,k+1)
     &           + (1.-r)*s*t*etapguessdata(i,j+1,k+1)


c..neutron chemical potential
      etanlin=  (1.-r)*(1.-s)*(1.-t)*etanguessdata(i,j,k)
     &           + r*(1.-s)*(1.-t)*etanguessdata(i+1,j,k)
     &           + r*s*(1.-t)*etanguessdata(i+1,j+1,k)
     &           + (1.-r)*s*(1.-t)*etanguessdata(i,j+1,k)
     &           + (1.-r)*(1.-s)*t*etanguessdata(i,j,k+1)
     &           + r*(1.-s)*t*etanguessdata(i+1,j,k+1)
     &           + r*s*t*etanguessdata(i+1,j+1,k+1)
     &           + (1.-r)*s*t*etanguessdata(i,j+1,k+1)


c..external proton fraction
      yplin=  (1.-r)*(1.-s)*(1.-t)*pprevdata(i,j,k)
     &           + r*(1.-s)*(1.-t)*pprevdata(i+1,j,k)
     &           + r*s*(1.-t)*pprevdata(i+1,j+1,k)
     &           + (1.-r)*s*(1.-t)*pprevdata(i,j+1,k)
     &           + (1.-r)*(1.-s)*t*pprevdata(i,j,k+1)
     &           + r*(1.-s)*t*pprevdata(i+1,j,k+1)
     &           + r*s*t*pprevdata(i+1,j+1,k+1)
     &           + (1.-r)*s*t*pprevdata(i,j+1,k+1)

      return
      end

