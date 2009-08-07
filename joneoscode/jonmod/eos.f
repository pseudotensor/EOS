C======================================================================
      Program NDAF_EOS

C======================================================================
C     This code computes the EOS in NDAF and ADAF
C     The unit is fundamentally Natural Unit.
C     /2005/05/06/ Kazunori Kohri
C======================================================================
      implicit none

c     Bring-in Kaz single global quantities
      include 'kazeos.dek'
      include 'kazeos.parms.dek'
C======================================================================
      integer i,j,k,l           !dummy indeces
C======================================================================
c     Loop quantities
      include 'kazeos.loopvars.dek'
      include 'kazeos.loopparms.dek'
c======================================================================
      integer ntotdatafile      !number of datafile
      parameter(ntotdatafile=14)
C======================================================================
c      integer computespecies


cccccccccccccccccccccccccccccccccc
c
c     Open files for first time (to be appended to later)
c
cccccccccccccccccccccccccccccccccc

      open(11,file='eos.dat',status='unknown')
      open(17,file='eosother.dat',status='unknown')
      close(11)
      close(17)


c     Below uses vars set by kazeos.loop*.dek (not globals)
      call outputkazheader()


c     Below uses vars set by kazeos.loop*.dek (not globals)
      call setupkazloops(intrhob,inttk,inthcm,inttdynorye,inttdynorynu)

      write(*,*) 'nrhob,ntk,nhcm,tdynorye,tdynorynu'
      write(*,*) nrhob,ntk,nhcm,ntdynorye,ntdynorynu


ccccccccccccccccccccccccccccccccccc
c     
c     Start loops
c     
ccccccccccccccccccccccccccccccccccc
      
      do ihcm=1, nhcm
       do itdynorynu=1, ntdynorynu
         do itdynorye=1, ntdynorye
               do itk=1, ntk
                  do irhob=1,nrhob


ccccccccccccccccccccccccccccccccccc
c     
c     Set independent globals
c     
ccccccccccccccccccccccccccccccccccc

                     call setkazindeps(irhob,itk,ihcm,itdynorye,itdynorynu
     1                    ,intrhob,inttk,inthcm,inttdynorye,inttdynorynu
     1                    ,nhcm, rhob,tk,hcm,tdynorye,tdynorynu)


ccccccccccccccccccccccccccccccccccc
c     
c     Call Kaz EOS, which fills Kaz single global quantities
c     
ccccccccccccccccccccccccccccccccccc
c     call kaz_eos(rhob,tk,hcm,tdynorye,tdynorynu)

                     computespecies=1
                     call kaz_eos()

ccccccccccccccccccccccccccccccccccc
c     
c     Write each state's solution to file
c     
ccccccccccccccccccccccccccccccccccc

                     call kaz_output()

                  end do
               end do
            end do
         end do
      end do

c     Final close of files
      close(11)
      close(17)

      stop
      end
