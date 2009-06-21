C=====================================================================
c
c     Assumes rhob,tk,hcm,tdyn globally set
c
c     This result only depends upon etanu if using whichrnpmethod==0
c     whichrnpmethod==1 sets Y_e directly and doesn't use neutrino processes to determine Y_e so is strictly based upon local charge neutrality
cc
c
cc      subroutine muekT(rho10,T11,amuekT,deg)
c      subroutine muekT(rho10,T11,xnuc,etae,deg)
c      subroutine muekT(rho10,T11,xnuc,etae)
      subroutine muekT(etae,etap,etan,etanu)

C=====================================================================
c     Solves for the quantity eta_e = amuekT = mu_e/kT
C=====================================================================
      implicit none

C======================================================================
c      common/tau/tautel,tautmu,tauttau,tauael,tauamu,tauatau
c      common /variable1/ tdyn,hcm,rhob,tk

      include 'kazeos.parms.dek'
      include 'kazeos1.dek'
      include 'kazeos3.dek'

c      real*8 xnuccalc ! function
      real*8 computefofetae ! function
C======================================================================
      real*8 rho10              !matter density in 10^10 g/cc
      real*8 T10                !temperature in 10^10 K
      real*8 etae,etap,etan,etanu     !degeneracy parameter for nu_e
c      real*8 xnuc               !fraction of free nucleons
C======================================================================
      real*8 T11                !temperature in 10^11 K
C=====================================================================
      real*8 ne              !electrn number density in MeV^3
      real*8 nemin,nemax
C=====================================================================
      real*8 f                  !evaluated_function in bisection method
C=====================================================================
c      real*8 QnpkT              !Q/T
      real*8 etaemin, etaemax   !minimum and maximum of trial values of etae
      real*8 fmin, fmax         !minimum and maximum of the evaluated_function
C=====================================================================
      real*8 eps4               !indicater of judgement for iteration
C=====================================================================
      real*8 deg                !degeneracy-parameter (=etae)
C=====================================================================
c      real*8 lowtempetae ! dimensionless
      integer etaeset
C=====================================================================
      real*8 mutot
      include 'const.dek'


cccccccccccccccccccccccc
c
c Set often-used variables
c
ccccccccccccccccccccccc

      rho10=rhob/1.d10          !in 10^{10} g/cc
      T11 = tk/1D11             !in 10^{10}K

c     Only done once since not dependent upon \eta_e
c     xnuc is not actually used if whichrnpmethod==1
c      xnuc=xnuccalc(rhob,tk)
c     xnuc_small=22.16d0*rho10**(-0.75d0)*T10**1.125d0*exp(-8.209d0/T10)
c     xnuc = min(1.d0, xnuc_small) !check1
      



c     Indicates etae not set yet
      etaeset = 0


      if(whichrnpmethod.eq.0) then
c     Set inital guess for yetot
         yetot=0.5
      else
         yetot=tdynorye
      end if
c     Otherwise yetot known and can compute everything

      



c      QnpkT=0.150d0/T11

c     Start with widely separated initial guess values so as to bracket
c     the solution
c     1E-20 needed for whichrnpmethod==0 rhob=1E2g/cc and T=1E13K
      etaemin=1.d-20
      etaemax=1.d20

c     Get fmin,fmax
      fmin = computefofetae(rho10,T11,etaemin,etap,etan,etanu)
      fmax = computefofetae(rho10,T11,etaemax,etap,etan,etanu)



c      write(*,*) 'fmin=',fmin,'fmax=',fmax
c      write(*,*) 'xnuc=',xnuc


      if (fmin*fmax.gt.0.d0) then
         write (*,*) ' solution not bracketed '
         write(6,108) rhob, tk
 108     format(20(1pe11.3,' '))
         write(*,*) 'etaemin=',etaemin,'etaemax=',etaemax
         write(*,*) 'fmin=',fmin,'fmax=',fmax
         write(*,*) 'nemin=',nemin,'nemax=',nemax

         return
      endif

c      write(*,*) 'etae=',etae



c     Solve for emuekT iteratively.  The equations we solve are:
c             n_p = n_e = (mu_e/hbar c)^3/(3 pi^2)
c             exp(mu_e-Q) = n_n/n_p

      if(etaeset.eq.0) then

         eps4=0.1d0
         do while (eps4.gt.1.d-11)

            etae=sqrt(etaemin*etaemax)

            f = computefofetae(rho10,T11,etae,etap,etan,etanu)

            if (f.gt.0.d0) then
               etaemax=etae
               fmax=f
            else
               etaemin=etae
               fmin=f
            endif

c            write(6,100) fmin,f,fmax,etaemin,etaemax,eps4

c     eps4=(etaemax-etaemin)/abs(etaemin)
            eps4=log10(etaemax/etaemin)/abs(log10(etaemin))
         end do

 20      continue

      end if


c     Solution for etae has converged.  Estimate degeneracy parameter
c     d = 0.178rho10/T11^3.  When d>1, we have a degenerate electron gas.
c      deg=etae

c     Final etae:
      etae=sqrt(etaemin*etaemax)

c      write(*,*) 'GOT HERE1.5',etae,etap,etan,etanu,yetot

 100  format(40(E27.16E5,' '))
c 100  format(1x,6F)
      
      return
      end




      




C=====================================================================
c
c     Assumes rhob,tk,hcm,tdyn globally set
c
      real*8 function computefofetae(rho10,T11,etae,etap,etan,etanu)

C=====================================================================
c     Solves for the quantity eta_e = amuekT = mu_e/kT
C=====================================================================
      implicit none

C======================================================================
c      common/tau/tautel,tautmu,tauttau,tauael,tauamu,tauatau
c      common /variable1/ tdyn,hcm,rhob,tk

      include 'kazeos1.dek'
      include 'kazeos3.dek'
      include 'kazeos.parms.dek'
C======================================================================
c      real*8 rhob,tk
C======================================================================
      real*8 rho10              !matter density in 10^10 g/cc
      real*8 T11                !temperature in 10^11 K
      real*8 etae,etap,etan,etanu     !degeneracy parameter for nu_e
C=====================================================================
      real*8 ne                 !electron number density in MeV^3
C=====================================================================
      real*8 rnp                !function
c      real*8 rnp_current        !substituted rnp
      real*8 yetot_current      ! present value of Y_e [total]
      real*8 nb ! baryon number density
C=====================================================================
c


      include 'const.dek'




      if(whichrnpmethod.eq.0) then

c     effective input of last version of yetot
         call kaz_species_etae(etae,etap,etan,etanu)
         yetot_current=yetot


      else if(whichrnpmethod.eq.1) then

c     Notice that tau's, xnuc, and etae not used in this case
         yetot_current = tdynorye

      end if


      nb = rhob/mb
      computefofetae = 1.0 - yetot_current*nb/(ne(T11,etae)*mev3tocc)
      


      return
      end





