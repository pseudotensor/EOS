

      subroutine event(kflavour, energy, result)
      implicit none
      save

c     include ...

c     Passed:
      integer kflavour
      real*4 energy
c     Returned:
      real*4 result


c      Do something

      result = 1.0d0*(kflavour)*energy;
      

      return
      end


