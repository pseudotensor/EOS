C======================================================================

C======================================================================
      real*8 function rnp(tdyn,tauael,taus,T11,etae,etap,etan,etanu)

C======================================================================
C     This computes the neutron to proton ratio
C     /10/24/2004/ K. Kohri
C======================================================================
      implicit none

      include 'kazeos.parms.dek'
c      include 'kazeos1.dek'
c      include 'kazeos3.dek'
      
C=====================================================================
c      common/params/alpha1,em,emdot,r
c      common/tau/tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus
c      common /check6/HoverR0

C======================================================================
c      common /variable1/ tdyn,hcm,rhob,tk

C======================================================================
c      real*8 tdyn               !dynamical timescale in sec
c      real*8 hcm                !disk half thickness in cm
c      real*8 rhob               !baryon denstiy in g/cc
c      real*8 tk                 !temperature in K


C======================================================================
c      real*8 tautel,tautmu,tauttau,tauael,tauamu,tauatau,taus
c      real*8 alpha1,em,emdot,r
c      real*8 HoverR0

C======================================================================
      real*8 T11                !temperature in 10^{11} K
      real*8 etae,etap,etan,etanu     !electron degeneracy parameter
      real*8 tdyn,tauael,taus
C======================================================================
      real*8 QnpkT              !QnpkT = (m_n-m_p)c^2/kT = 1.29MeV/kT

C======================================================================
c      real*8 gammapeint         !rate p + e^- --> n + nue  in /sec
c      real*8 gammaneint         !rate n + e^+ --> p + nuebar in /sec
c      real*8 gammanint          !rate n --> p + e ^- + nuebar in /sec

C=====================================================================
      real*8 OmegaK             !sec^-1
      real*8 tauacc             !accretion timescale in sec
      real*8 Gammaacc           !1/tauacc in sec^-1
      real*8 ratioGG(3)         !Gammareac/Gammaacc

C======================================================================
      real*8 tmp(100)

C======================================================================
      integer ipower(2)
      data ipower /2.d0, 2.d0/

C======================================================================
      real*8 qdotNe, qdotpair, qdotbrem, qdotplasmon

C======================================================================
c      real*8 tauael_ele
c      real*8 tauael_antiele

C======================================================================
      real*8 qminusel_ele
      real*8 qminusel_antiele
      real*8 qdotNe_ele
      real*8 qdotNe_antiele

c======================================================================
c      real*8 taus
c      real*8 tauabs, tauscatt

C======================================================================
c      real*8 qminusel, qminusmu, qminustau

C======================================================================
      real*8 gammape
      real*8 gammane
      real*8 gamman

C======================================================================
      real*8 rho10
C======================================================================
      real*8 ex

C======================================================================
      real*8 expind             !expornent index for smooth connection
      parameter(expind=10.d0)

C======================================================================
      include 'const.dek'
C======================================================================




      if(whichrnpmethod.eq.0) then


         tauacc=tdyn
         QnpkT=0.150d0/T11
         Gammaacc= 1.d0/tauacc

ccccccccccccccc
c
c     Below rates are not corrected for optical depth effects (and only applies to *free* nucleons)
c
ccccccccccccccc

         tmp(1)=dim(gammape(T11,etae),1.d-40)
         tmp(2)=dim(gammane(T11,etae)+gamman(T11,etae),1.d-40)


ccccccccccccccc
c
c  Set ratio of dY_e/dt to accretion timescale
c
ccccccccccccccc

         ratioGG(1)=tmp(1)/Gammaacc

c     ratioGG(2)=tmp(2)/Gammaacc
c     ratioGG(3)= (1.d0/Gammaacc) / (1.d-4*(3.d0*T11)**(-5.d0))

c     if (ratioGG(1).lt.1.d0) then
c     rnp=1.d0
c     go to 20
c     end if

         tmp(30)=2.d0/3.d0      !optical depth

         tmp(10)=min(tauael/tmp(30),(taus)/tmp(30))**real(ipower(1))
c     tmp(10)=min(tauael_ele/tmp(30),tauael_antiele/tmp(30),(taus)/tmp(30))**real(ipower(1))

         tmp(11)=exp(-QnpkT+etae)*tmp(10)/(1.d0+tmp(10))
         tmp(12)=tmp(1)/tmp(2)*1.d0/(1.d0+tmp(10))
         rnp = dim(tmp(11)+tmp(12), 1.d-40)

         if (rnp.gt.1.d3) then
            rnp =1.d3
         end if
c     GODMARK: (what's going on here?)
c     Assume Y_e<1/1.5 strictly since code returns rnp=0 and makes no sense
c         if (rnp.lt.0.5) then
c            rnp = 0.5
c         end if
         if (rnp.lt.1.0) then
            rnp = 1.0
         end if


c      write(*,*) 'rnpissue=',xnuc,yefreelocal,yeboundlocal,rnp


      else if(whichrnpmethod.eq.1) then

         write(*,*) 'Should not be in rnp() with whichrnpmethod==1'
         stop

      end if
      



 20   continue

      return
      end























c     GODMARK: since depends upon \tau(H), depends upon H.
c     So in HARM, must store gamma??globals and optical depth and reconstruct
c     this function and the above dyedtcalc() function
c     This takes into account free vs. bound protons/neutrons
C======================================================================
      subroutine dyedtcalc(tdyn,yelocal,gamman2p,gammap2n,dyedt,thermalye)
C======================================================================
C     This computes the term dY_e/dt *naively* without optical depth effects generally accounted for
C     JCM
C======================================================================
      implicit none

      include 'kazeos.parms.dek'
c      include 'kazeos1.dek' ! we use global gamma's

c     Passed:
      real*8 tdyn
      real*8 yelocal            ! Y_e total
      real*8 gammap2n,gamman2p

c     Returned:
      real*8 dyedt,thermalye

c     Local:
      real*8 thermalnpratio
C======================================================================



      if(whichrnpmethod.eq.0) then

         thermalye=yelocal

c     If perfectly thermalized, then dY_e/dt = 0 locally
c     Kaz actually tries to compute some incomplete thermalization using tdyn
c     But I don't understand what reference (starting) Y_e is for that
c     Probably Y_e=0.5 is starting point and assume accrete/fix Y_e
c     on the timescale of tdyn
c         dyedtcalc = 0.0
c     So in reality this is not 0 for Kaz, it's really
c     Assume rnp already computd so thermalye computed
         dyedt = (thermalye-0.5)/tdyn


      else if(whichrnpmethod.eq.1) then

c     Assume optical depths have been stored globally if here
c     Only electron neutrino optical depth used to determine electron exchanges
c     Assume expential supression of those rates dependent on large optical depth

c     Note that for Y_e=0 that gamman2p only increases Y_e, as required
c     Assume dyedtcalc is used in such a way that new Y_e>=0 or some minye
         dyedt = gamman2p - (gamman2p + gammap2n)*yelocal

c         write(*,*) dyedtcalc,gamman2pglobal,gammap2nglobal,yelocal

c     Below does not assume neutrinos are perfectly thermalized since
c     only case if \tau_{\nu}>>1
         thermalnpratio = gammap2n/gamman2p
c     This is some kind of global thermal value, not an absolute thermal value
c     that would be obtained by computing \eta_e assuming the thermal Y_e
         thermalye = 1.0/(1.0+thermalnpratio)

      end if





      return
      end








c     This function does not depend on optical depths and so hasn't accounted for them yet
c     this function and the above dyedtcalc() function
c     This takes into account free vs. bound protons/neutrons
C======================================================================
      subroutine computegammarates(T11,etae,etap,etan
     1     ,aheav,zheav,npheav,nnheav
     1     ,npfree,nnfree,nptotal,nntotal
     1     ,gamman2p,gammap2n
     1     ,gammape_val,gammaAe_val,gammapnu_val,gammapenu_val,gamman_val,gammane_val,gammannu_val)
C======================================================================
C======================================================================
      implicit none

      include 'kazeos.parms.dek'
c      include 'kazeos1.dek'

c     Passed:
      real*8 T11                !temperature in 10^{11} K
      real*8 etae,etap,etan     !electron degeneracy parameter
      real*8 aheav,zheav,npheav,nnheav
     1     ,npfree,nnfree,nptotal,nntotal
     1     ,gamman2p,gammap2n
     1     ,gammape_val,gammaAe_val,gammapnu_val,gammapenu_val,gamman_val,gammane_val,gammannu_val
C======================================================================
      real*8 gammannu,gammapnu,gammapenu,gammane,gamman,gammape,gammaAe ! functions
C======================================================================
      real*8 taulocal
      real*8 ex ! function
C======================================================================
      real*8 expind             !exponent index for smooth connection
      parameter(expind=1.d0)
c      real*8 nbtotal
      real*8 nA
C======================================================================

c     GODMARK: What to do?
c         taulocal = ddim(ex(-expind/(tautel+1.d-40)),0.0)
c     GODMARK: How should I combine the optical depths?
c         taulocal = ddim((1.0-ex(-(tautel+tauael))),0.0)

c     JCM: Below are per baryon rates, so processes that only operate on free nucleons or only bound nuclei need to be accounted for
c         nbtotal = rhob/mb

      
         

c     Terms related to gammap2n
         gammape_val   = (npfree/nptotal)*ddim(gammape(T11,etae),1.d-40)
c     JCM: below takes into account fact that gammaAe() returns rate per nuclei with atomic mass A
c     gammaAe() n_A = +dn_n/dt contribution and -dn_p/dt contribution
c     Write as per proton in nuclei since each unit converts 1 proton/sec to 1 neutron/sec
         if(aheav>=0.5) then
            nA = ((npheav+nnheav)/aheav)
            gammaAe_val   = (nA/nptotal)*dim(gammaAe(T11,etae,etap,etan,aheav,zheav),1.d-40)
         else
            gammaAe_val = 0.0
         end if

         gammapnu_val  = (npfree/nptotal)*ddim(gammapnu(T11,etae),1.d-40)
         gammapenu_val = (npfree/nptotal)*ddim(gammapenu(T11,etae),1.d-40)

c         gammap2n = gammape_val+gammaAe_val + taulocal*(gammapnu_val+gammapenu_val)
c     Use only optically thin rates
         gammap2n = gammape_val+gammaAe_val

c     Terms related to gamman2p
         gamman_val   =  (nnfree/nntotal)*ddim(gamman(T11,etae),1.d-40)
         gammane_val  =  (nnfree/nntotal)*ddim(gammane(T11,etae),1.d-40)
         gammannu_val =  (nnfree/nntotal)*ddim(gammannu(T11,etae),1.d-40)

c         gamman2p = gamman_val + gammane_val + taulocal*(gammannu_val)
c     Use only optically thin rates
         gamman2p = gamman_val + gammane_val

c         write(*,*) 'gammas=',gammape_val,gammaAe_val,etae,etap,etan,kazabar,kazzbar


         return
         end

