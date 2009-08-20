c     
c     
c     
c     compile doing: ifort -cm -w90 -132 -o jon_loop jon_loop.f
c     
c     Run by doing:
c     rm -rf eosdetails.dat nohup.out eosother.dat eos.dat eoscoulomb.dat eosazbar.dat eos.head
c     ./helmeos.exe
c
c     Or if want to do chunk 1 out of 30 (1 is first chunk), do:
c     
c     echo "1 30" > eoschunk.dat   
c     rm -rf eosdetails.dat nohup.out eosother.dat eos.dat eoscoulomb.dat eosazbar.dat eos.head
c     ./helmeos.exe
c     
c     
      program teos
      implicit none
c     
      save
c     
      include 'const.dek'
c..   bring in the data structure for HELM or TIMMES EOS
      include 'vector_eos.dek'
c     next gets what type of data we are to output
      include 'kazeos.parms.dek'
      include 'kazeos.dek'
c     next 2 are for loop parameters
      include 'kazeos.loopvars.dek'
      include 'kazeos.loopparms.dek'

      integer whichyemethod
      integer truentdynorye
      double precision yein
      integer CHUNK,TOTALCHUNKS,chunkiter,totalchunkiter
      integer chunkiterlow,chunkiterhigh
      real*8 abarlocal,abarboundlocal


c     0 = assign A,Z (rho,T)
c     1 = use real ye from loop values
c     Notice that ye can be treated specially, but not ynu
      whichyemethod=1


c     Set diagnostic for total number of converged LSEOS solutions to 0
      totalnumconverged=0
      totalnumgoodconverged=0



c     Make sure enough memory for storage AND for temporary LSEOS AND temporary HELMEOS space for electron EOS
      if(nrowmax.lt.nrhob+nrowextra) then
         write(*,*) 'nrowmax=',nrowmax,'nrhob+nrowextra=',nrhob+nrowextra
         stop 'nrowmax not large enough'
      end if




c     Open data files first time, to be append to later
      open(11,file='eos.dat',status='unknown')
      open(17,file='eosother.dat',status='unknown')
      open(18,file='eoscoulomb.dat',status='unknown')
      open(19,file='eosazbar.dat',status='unknown')
      open(20,file='eosdetails.dat',status='unknown')
      close(11)
      close(17)
      close(18)
      close(19)
      close(20)


c     Setup true and fake nuclear offset
      call setup_lsoffset()
      

c     Below uses vars set by kazeos.loop*.dek (not globals)
      call outputkazheader()


c     Below uses vars set by kazeos.loop*.dek (not globals)
      call setupkazloops(intrhob,inttk,inthcm,inttdynorye,inttdynorynu)

      if(whichyemethod.eq.0) then
         inttdynorye=0.0
         truentdynorye=1
      else
         truentdynorye=ntdynorye
      end if

c     write(*,*) intrhob,inttk,inthcm,inttdynorye
c     write(*,*) nrhob,ntk,nhcm,ntdynorye


c     Read-in any chunk parameters
      open(680,file='eoschunk.dat',status='unknown')
      read(680,*,END=134) CHUNK,TOTALCHUNKS
      goto 135
 134  CHUNK=1
      TOTALCHUNKS=1
 135  continue

      totalchunkiter=nhcm*ntdynorynu*truentdynorye*ntk
c     Ensure total number of chunks trying to do isn't larger than number that can be done
      if(totalchunkiter.lt.TOTALCHUNKS) then
         write(*,*) 'Limited chunks from/to',TOTALCHUNKS,totalchunkiter
         TOTALCHUNKS=totalchunkiter
      end if

c     Set lower and upper chunkiter range for given CHUNK and TOTALCHUNKS
      chunkiterlow=INT((1.0D0*CHUNK-1.0D0)/(1.0D0*TOTALCHUNKS)*(1.0D0*totalchunkiter)+1.0D0)
      chunkiterhigh=INT((1.0D0*CHUNK)/(1.0D0*TOTALCHUNKS)*(1.0D0*totalchunkiter))

c     Ensure lower and upper bounds grabbed
c     At boundaries within chunks, should consistently get INT()'ed value for any input CHUNK or TOTALCHUNKS
      if(CHUNK.eq.1) then
         chunkiterlow = 1
      end if

      if(CHUNK.eq.TOTALCHUNKS) then
         chunkiterhigh = totalchunkiter
      end if

      write(*,*) 'CHUNK,TOTALCHUNKS',CHUNK,TOTALCHUNKS
      if(TOTALCHUNKS .eq. 1) then
         write(*,*) 'If create eoschunk.dat with 2 integers'
         write(*,*) 'then will do CHUNK out of TOTALCHUNKS'
      end if
      write(*,*) 'Maximum TOTALCHUNKS=',totalchunkiter
      write(*,*) 'chunkiterlow,chunkiterhigh',chunkiterlow,chunkiterhigh



         


cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     loop over temperatures since large ntk and nrhob is too much for continuous memory
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Initialize chunkiter
      chunkiter=1

      do ihcm=1, nhcm
         do itdynorynu=1, ntdynorynu
            do itdynorye=1, truentdynorye
               do itk=1, ntk

c     nrhob loop is an internal subloop later, so can not chunk it, but chunk rest


                  if(chunkiter .ge. chunkiterlow .AND. chunkiter .le. chunkiterhigh) then



c     do inner nrhob subloop
                     do irhob=1,nrhob
                        

c     Set global independent variables: rhob,tk,hcm,tdynorye
                        call setkazindeps(
     1                       irhob,itk,ihcm,itdynorye,itdynorynu
     1                       ,intrhob,inttk,inthcm,inttdynorye,inttdynorynu
     1                       ,rhob,tk,hcm,tdynorye,tdynorynu)
                        

c     define storage mapping function
c     1D input vector (saves memory)
                        index=1+(irhob-1)
                        temp_row(index) = tk
                        den_row (index) = rhob
                        hcm_row (index) = hcm
                        
c     From input vectors of temperature and density, compute fractional species

c     Set ye<0 to tell azbarset to compute using fitting function
c     Tell azbarset to use ye from loop parameters
c     yein = tdynorye




cccccccccccccccccccccccccccccccccccccccc
c     
c     Here yein is determined from A,Z given rho,T
c     
cccccccccccccccccccccccccccccccccccccccc
                        if(whichyemethod.eq.0) then
                           yein=-1
                           call azbarset1(den_row(index),temp_row(index),yein
     1                          ,abar_row(index)
     1                          ,abarnum_row(index)
     1                          ,zbar_row(index)
     1                          ,abarbound_row(index)
     1                          )
                           tdynorye = zbar_row(index)/abarnum_row(index)
                        else

cccccccccccccccccccccccccccccccccccccccc
c     
c     Here yein sets A,Z and nuclear EOS is expected to better handle what A,Z to choose
c     
cccccccccccccccccccccccccccccccccccccccc

                           if(1.eq.0) then
                              yein = tdynorye
c     Need to determine how to set A and Z given Y_e=Z/A
c     For now assume simplest thing, free nucleons with A=1 and effective Z to give ye
c     Assume nuclear EOS will set A,Z for given Ye
                              abar_row(index)=1.0
                              abarnum_row(index)=1.0
                              zbar_row(index)=abarnum_row(index)*yein
                              abarbound_row(index)=1.0
                           end if

c     Use below now so when outside nuclear EOS (or not in NSE) can have stellar-model-motivated A and Z
                           if(1.eq.1) then
                              yein = tdynorye
c     Setup A and Z in case nuclear EOS is out of range or not provided (whichyemethod.eq.1)

c     Get non-nuclear <A> and <Abound>
                              call nonnuclearA(den_row(index),temp_row(index),abarlocal,abarboundlocal)
c     Assign final <A>, <Z>
                              abar_row(index)=abarlocal
                              abarnum_row(index)=abar_row(index)
                              zbar_row(index)=abarnum_row(index)*yein
c     Assign <Abound>
                              abarbound_row(index)=abarboundlocal
                           end if


                        end if

c     Now set Y_e for Kaz processing parts of jon_lsbox.f
                        tdynorye_row (index) = tdynorye




c     Assume doesn't matter what tdynorynu is if whichynumethod==0
                        tdynorynu_row (index) = tdynorynu



                     end do




ccccccccccccccccccccccccccccc
c     
c     Now that den,temp,abar,abarnum, and zbar are set, can call the EOS
c     
ccccccccccccccccccccccccccccc

c     Setup anyeos call
                     jlo_eos = 1
                     jhi_eos = nrhob

                     call anyeos



                     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     
c..   write out the results (if doing multiple calls, appends)
c     
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     call pretty_eos_out('helm:  ')


                  end if
c     First value of chunkiter will be 1 since started at 1, and here we iterate chunkiter
                  chunkiter=chunkiter+1

c     End outer loop
               end do
            end do
         end do
      end do

c     Report convergence for LSEOS
      write(*,*) 'Total LS converged:',totalnumconverged,'out of',nrhob*ntk
      write(*,*) 'Total LS good converged:',
     1     totalnumgoodconverged,'out of',nrhob*ntk


      stop 'normal termination'

 03   format(1x,t2,a6,I3,t22,a6,I3)

 100  format(40(1pe26.15,' '))




      end   




