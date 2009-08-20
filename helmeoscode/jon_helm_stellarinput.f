c
c
c
c compile doing: ifort -cm -w90 -132 -o jon_loop jon_loop.f
c
c
c
      program teos
      implicit none
c
      save
c
      include 'const.dek'
      include 'eosparms.f'
c..bring in the data structure for HELM or TIMMES EOS
      include 'vector_eos.dek'
c..bring in the data structure for KAZ-like equantities
      include 'kazeos.dek'
c     next gets what type of data we are to output
      include 'kazeos.parms.dek'

      real*8 yein,abarin,abarboundin
      integer ionmax
      parameter (ionmax=9) 
      double precision xmass(ionmax)
      integer j,numlines


c     Set diagnostic for total number of converged LSEOS solutions to 0
      totalnumconverged=0
      totalnumgoodconverged=0


      if(whichhcmmethod.ne.1) then
         write(*,*) 'whichhcmmethod must be 1'
         stop
      end if



c     Open data files first time, to be append to later (i.e. not appending to previously existing file)
      open(11,file='eos.dat',status='unknown')
      open(17,file='eosother.dat',status='unknown')
      open(18,file='eoscoulomb.dat',status='unknown')
      open(19,file='eosazbar.dat',status='unknown')
      open(20,file='eosdetails.dat',status='unknown')
      close(11)
      close(17)
      close(18)
      close(19)
      close(120)

c      open(21,file='eosdebug.dat',status='unknown')
c      close(21)
 

c     Setup true and fake nuclear offset
      call setup_lsoffset()


c     Below uses vars set by kazeos.loop*.dek (not globals)
c     Note that most of this data isn't used, but useful to keep SM macros simple
      call outputkazheader()





cccccccccccccccccccccccccccccccccccccccc
c
c     READ IN DATA HEADER
c
cccccccccccccccccccccccccccccccccccccccc

      open (unit=3,file='stellarmodel.head')
      read (3,*) numlines
      close (3)
         

c     Make sure enough memory for storage AND for temporary LSEOS AND temporary HELMEOS space for electron EOS
      if(nrowmax.lt.numlines+nrowextra) then
         write(*,*) 'nrowmax=',nrowmax,'numlines+nrowextra=',numlines+nrowextra
         stop 'nrowmax not large enough'
      end if




c..Setup the pipeline
      jlo_eos = 1
c SUPERDEBUG GODMARK
      jhi_eos = numlines
c      jhi_eos = 2


cccccccccccccccccccccccccccccccccccccccc
c
c     READ IN DATA
c
cccccccccccccccccccccccccccccccccccccccc

      open (unit=3,file='stellarmodel.dat')

ccccccccccccccccccccccccccccccccccccccc
c
c     Convert rhob,tk,X into rhob,tk,u,p,s
c
cccccccccccccccccccccccccccccccccccccccc

c     DEBUG:
c      write(21,*) 'numlines=',numlines,jlo_eos,jhi_eos

c..   set the input vector
      do j=jlo_eos,jhi_eos

c     yein is *total* Ye = n_e / (n_p + n_n).  That is, n_e = n_p, whereas Ye(free) = n_e(free)/(n_p(free)+n_n(free))
c     which is different since n_e(free)=n_p(free)!=n_p
c     read-in rho, temp, X_H, X_He, X_C, X_O, X_N, X_Mg, X_Si, X_Fe
         read (3,*) rhob, tk, tdynorye, abarin, abarboundin, xmass(2)
     1   ,xmass(3), xmass(4), xmass(5), xmass(6), xmass(7), xmass(8)
     1   ,xmass(9)
     1   ,hcm, tdynorynu

c     Below for HELM that fails if outside table
         if(whicheleeos.eq.0) then
            if(tdynorye*rhob.le.2E-12) then
               rhob = 2E-12/tdynorye
            end if
            
            if(tk.le.2E3) then
               tk = 2E3
            end if
         end if


c     Set ??_row independent vars
         hcm_row(j) = hcm+1E-20
         den_row(j)=rhob
         temp_row(j)=tk

c     call azbarset
c     detailed stellar model:
c     Note that inputted A,Z not used if inside nuclear EOS domain and NSE holds
         yein=tdynorye
         if(abarin.lt.0.0) then
            call azbarset2(rhob,tk,yein,ionmax,xmass,abar_row(j),abarnum_row(j),zbar_row(j),abarbound_row(j))
         else
            abar_row(j) = abarin
            abarnum_row(j) = abarin
            zbar_row(j) = yein*abarin
            abarbound_row(j) = abarboundin
         end if
c     simple ye and mutot fits
c     call azbarset1(rhob,tk,yein,abar_row(j),abarnum_row(j),zbar_row(j))

c         write(21,*) 'j=',j,jlo_eos,jhi_eos
c         write(21,*) rhob,tk,yein,abar_row(j),abarnum_row(j),zbar_row(j)

c     Now read in (realize Ynu=1E-4 may be HUGE for optical-depth supressed Ynu)
c         write(*,*) 'abarcompute',rhob,tk,yein,ionmax,xmass,abar_row(j),abarnum_row(j),zbar_row(j),abarbound_row(j)
c         write(*,*) 'abarcompute',abar_row(j),abarnum_row(j)

c     Now set Y_e for Kaz processing parts of jon_lsbox.f
         tdynorye_row(j) = tdynorye

         tdynorynu_row(j) = tdynorynu


      end do
      

c     Close the input file
      close(3)
c      close(21)




ccccccccccccccccccccccccccccc
c      
c     Now that den,temp,abar,abarnum, and zbar are set, can call the EOS
c
ccccccccccccccccccccccccccccc
      call anyeos()


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c..write out the results (if doing multiple calls, appends)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call pretty_eos_out('helm:  ')


c     Report convergence for LSEOS
      write(*,*) 'Total LS converged:',totalnumconverged,'out of',numlines
      write(*,*) 'Total LS good converged:',
     1            totalnumgoodconverged,'out of',numlines




      stop 'normal termination'

03    format(1x,t2,a6,I3,t22,a6,I3)


100   format(40(1pe26.15,' '))

      end   



