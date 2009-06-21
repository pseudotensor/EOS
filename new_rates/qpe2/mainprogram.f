C======================================================================
      Program data_taking

C======================================================================
      implicit none

C======================================================================
      integer i,j,k

C======================================================================
      real*8 T11                !temperature in 10^11 K
      integer nT
      parameter(nT=701)
      real*8 intT
      real*8 Tmax
      parameter(Tmax=1.d2)
      real*8 Tmin
      parameter(Tmin=1.d-5)
      
C======================================================================
      real*8 etae               !electron degeneracy parameter \mue/T
      integer neta
      parameter(neta=501)
      real*8 inteta
      real*8 etamax
      parameter(etamax=1.d2)
      real*8 etamin
      parameter(etamin=1.d-3)

C=====================================================================
      real*8 mev3tocc           !1 MeV^3 = (mev3tocc) /cc
      parameter(mev3tocc=1.3014d32)   !1 MeV^3 = (mev3tocc) /cc

C=====================================================================
      real*8 qpe2                 !function

C=====================================================================
      real sint(nT,neta)

C======================================================================
      real*8 tmp(20)

C======================================================================
      open(11,file="qpe21.dat",status="unknown")
      open(12,file="imqpe",form='unformatted',status="unknown")

      intT=log10(Tmax/Tmin)/real(nT-1)
      inteta=log10(etamax/etamin)/real(neta-1)
      do i=1,nT
         T11=Tmin*10.d0**(intT*real(i-1))
         do j=1,neta
            etae=etamin*10.d0**(inteta*real(j-1))
            tmp(1)=qpe2(T11,etae)
            write(11,100) T11, etae, tmp(1)
            sint(i,j)=sngl( log10(tmp(1)) )
            
 100        format(20(1pe11.3,' '))
         end do
      end do

      write(12) nT,neta
      write(12) sint
      close(12)

      close(11)

      stop
      end

