C======================================================================
      Program data_taking

C======================================================================
      implicit none

C======================================================================
      integer i,j,k

C======================================================================
      real*8 T11                !temperature in 10^11 K
      integer nT
      parameter(nT=701)        !check
      real*8 intT
      real*8 Tmax
      parameter(Tmax=1.d2)
      real*8 Tmin
      parameter(Tmin=1.d-5)
      
C======================================================================
      real*8 etae               !electron degeneracy parameter \mue/T
      integer neta
      parameter(neta=501)      !check
c      parameter(neta=6)
      real*8 inteta
      real*8 etamax
      parameter(etamax=1.d2)
      real*8 etamin
      parameter(etamin=1.d-3)

C=====================================================================
      real*8 mev4toecc           !1 MeV^4 = (mev4toecc) erg/cc
      parameter(mev4toecc=2.085d26)

C=====================================================================
      real*8 pe                 !function

C=====================================================================
      real sint(nT,neta)

C======================================================================
      real*8 tmp(20)

C======================================================================
      open(1,file="pe1.dat",status="unknown")
      open(2,file="impe",form='unformatted',status="unknown")

      intT=log10(Tmax/Tmin)/real(nT-1)
      inteta=log10(etamax/etamin)/real(neta-1)
      do i=1,nT
         T11=Tmin*10.d0**(intT*real(i-1))
         do j=1,neta
            etae=etamin*10.d0**(inteta*real(j-1))
            tmp(1)=pe(T11,etae)
            write(1,100) T11, etae, tmp(1)
            sint(i,j)=sngl( log10(tmp(1))+log10(mev4toecc) )
            
 100        format(20(1pe11.3,' '))
         end do
      end do

      write(2) nT,neta
      write(2) sint
      close(2)

      close(1)

      stop
      end

