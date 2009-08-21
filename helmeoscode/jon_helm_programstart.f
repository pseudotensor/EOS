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

      call teos_helm()

      stop 'normal termination'

      end   




