c$$$ From Kaz email "my subroutines" from 11/1/2007
c$$$ However, please be careful to use them because of the following
c$$$ reasons,
c$$$ 
c$$$  i) They are slow to be performed.
c$$$ 
c$$$ ii) In calculations of the reaction rates and the emission rates
c$$$     where neutrinos appear in initial states such as n + \nu_e --> p + e^- 
c$$$     (not final states such as n + e^+ --> p + \bar{\nu}_e ), the
c$$$     subroutines assume the perfect Fermi-Dirac distribution of neutrinos
c$$$     and anti-neutrinos.
c$$$ 
c$$$ So, you cannot use the reaction rates of (\Gamma_{n + \nu_e --> p +
c$$$ e^-} and \Gamma_{n + \nu_e + e^+--> p } and \Gamma_{p + \bar{\nu}_e 
c$$$ --> n + e^+ } when the neutrinos are optically thin.  That was because
c$$$ I did not solve the energy transfer Boltzmann equations of neutrinos
c$$$ which has been main crucial topics among experts of the core-collapsed 
c$$$ SNe and difficult to calculate them. I have never tried to do
c$$$ that. Instead of solving the Boltzmann equations correctly, I simply
c$$$ adopted the approximate method to treat both the optically thick and
c$$$ thin cases by introducing TDYN and Hcm and so on, and then including 
c$$$ and excluding the initial state neutrinos although we know that
c$$$ approximations should be correct in perfectly-optically thick and thin
c$$$ limit cases, respectively. As I told you that this idle approximation
c$$$ was discussed in Appendix B of Kohri-Narayan-Piran (05). 
c
c Jon interpretation:
c
c Rates below are not to be used in optically thin regime:
c 
c gammannu  = \Gamma_{n + \nu_e --> p + e^- }
c
c gammapnu  = \Gamma_{p + \bar{\nu}_e --> n + e^+ }
c
c gammapenu = \Gamma_{p + e^- + \bar{\nu}_e  --> n } (strong even though 3-body)
c
c ????????? = \Gamma_{n + \nu_e + e^+ --> p } (rate is 3-body and weak)
c
c
c gammane   = \Gamma_{n + e^+ --> p + \bar{\nu}_e }
c 
c gamman    = \Gamma_{n --> p + e^- + \bar{\nu}_e }
c
c gammape   = \Gamma_{p + e^- --> n + \nu_e }
c
c
c gammap2n = gammape + (if optically thick)*(gammapnu+gammapenu)
c
c gamman2p = gamman + gammane + (if optically thick)*(gammannu)
c
c



C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammannu(T11,etae,etanu,extrablock)

      implicit none

C----------PARAMTER.
C======================================================================
      real*8 T11,etae,etanu,extrablock
C======================================================================
      real*8 dum(100)            !dummy parameters
C======================================================================
      integer ratepow
C======================================================================
      real*8 gammannu_gen ! function
C======================================================================


      include 'const.dek'

      ratepow=2
      gammannu = gammannu_gen(ratepow,T11,etae,etanu,extrablock)

      return   !check2
      end




C======================================================================
c      Program reaction_rate_p_electron  !check2
      real*8 function gammannu_gen(ratepow,T11,etae,etanu,extrablock)
c GODMARK: This rate doesn't account for neutrino energy
c
c
c

C======================================================================
C     This function computes the reaction rate through p + electron --> n
C     + nue , as functions of temperature and the degeneracy parameter
C     of electrons 
C     /04/10/4/ K. Kohri
C======================================================================
      implicit none

C=====================================================================
      common /dist1/etae1,metilde,Qtilde,etanu1,ratepow1,extrablock1

C=====================================================================
C----------PARAMTER.
      integer   iter
      PARAMETER (iter=50)          !Number of gaussian quads.
      real*8 xintd              !function
      real*8 distintegrate      !function
C======================================================================
      real*8 T11                !temperature in 10^11 K
      real*8 mue11              !electron chemical potential in 10^11 K
      real*8 muemev             !electron chemical potential in MeV

C======================================================================
      real*8 tmev               !temperature in MeV
      real*8 etae,etanu              !electron degeneracy parameter \mue/T
      real*8 extrablock
      integer ratepow
C======================================================================
      real*8 ne                 !check1

C======================================================================
c      real*8 gammannu                !number deinsty of e- in MeV^3 !check2
      real*8 npositron          !number density of e+ in MeV^3

C======================================================================
c      real*8 memev              !electron mass in MeV
c      parameter(memev=0.511d0)
c      parameter(memev=0.d0)  !check1
      real*8 metilde           !me/T

C======================================================================
c      real*8 Qmev                  !qvalue of weak interaction between N,P
c      parameter(Qmev=1.29d0)
      real*8 Qtilde
      real*8 ratepow1
C======================================================================
c      real*8 mev2K              ! 1MeV = mev2K  K
c      parameter(mev2K=1.1605d10) ! 1MeV = mev2K  K
c      real*8 PI
c      parameter(PI=3.14159265d0)
c      real*8 mevtosec           !1 MeV = (mevtosec) /sec
c      parameter(mevtosec=1.519d21) !1 MeV = (mevtosec) /sec

C======================================================================
c      real*8 ga                 !axial vector coupling
c      parameter(ga=1.39d0)
c      real*8 GF                 !Fermi coupling constant
c      parameter(GF=1.1664d-11)  !in MeV^-2

C======================================================================
      real*8 dum(100)            !dummy parameters
      real*8 etae1,etanu1              !electron degeneracy parameter \mue/T
      real*8 extrablock1
C======================================================================
      real*8 np                 !proton number density in MeV^3

C======================================================================
      real*8 neint                 !function of electron number density in MeV^3

C======================================================================
      external integgammannu

C======================================================================
c      read(*,*) T11,etae        !check2



      include 'const.dek'


      etae1=etae
      etanu1=etanu
      extrablock1=extrablock
      ratepow1=dble(ratepow)

      tmev = T11*1.d11/mev2K
      Qtilde=Qmev/tmev
      metilde=memev/tmev

c      dum(1)=xintd(0.d0,1.d0,integgammannu,iter)
c      dum(2)=xintd(1.d0,10.d0,integgammannu,iter)
c      dum(3)=xintd(10.d0,100.d0,integgammannu,iter)
c      dum(4)=xintd(100.d0,1000.d0,integgammannu,iter)
c      dum(5) = dum(1) + dum(2) + dum(3) + dum(4)
c JCM:
      dum(5) = distintegrate(integgammannu)


      dum(10)=GF**2/(2.d0*PI**3)*(1.d0+3.d0*ga**2)*tmev**5 !in MeV
      if(ratepow.eq.2) then
         dum(11) =  dum(10) * mevtosec ! in 1/sec
      else if(ratepow.eq.3) then
         dum(11) =  dum(10) * tmev * mev2toergs ! in erg/sec
      end if

      gammannu_gen = dum(11) * dum(5)

c      write(6,100) tmev, etae, gammannu  !check2
 100  format(20(1pe11.3,' '))

c      stop
      return   !check2
      end

C======================================================================
      real*8 function integgammannu(x)

C======================================================================
c     These are functions of distributions of electrons and positrons
C     /04/10/03/   K. Kohri
C======================================================================
      implicit none

C======================================================================
      common /dist1/etae1,metilde,Qtilde,etanu1,ratepow1,extrablock1

C======================================================================
      real*8 x                  !variable of integral 

C======================================================================
      real*8 etae1,etanu1              !electron degeneracy parameter
      real*8 extrablock1
      real*8 metilde           !me/T
      real*8 Qtilde             !Q/T
      real*8 ratepow1
C======================================================================
c     real*8 integgammannu       !integrand of electrons distribution function

C======================================================================
      real*8 tmp(20)

C======================================================================
      real*8 ex                 !function

C======================================================================
C----------------------------------------------------------------------
C     JCM: Note: x = pc/(kb*T)
      tmp(12) = ddim((x+Qtilde)**2 -metilde**2,0.0d0)
      tmp(1)= dsqrt( tmp(12)) * (x+Qtilde)* x**(ratepow1)
c     Pauli Blocking factor
      tmp(2)=dexp(x-etanu1) + 1.d0

c     JCM: Below should be (according to KNP05) (correct:)
      tmp(3)=1.d0/(  1.d0 + dexp(-(x+Qtilde - etae1)))
c     Was:
c      tmp(3)=1.d0/(  1.d0 + dexp(-(x+metilde - etae1)))

c     As KNP05, but with extra blocking factor
      tmp(4)=1.d0*extrablock1/(  1.d0 + dexp(x+Qtilde - etae1))
      tmp(3)=1.d0-tmp(4)

      integgammannu=tmp(1) / tmp(2)*tmp(3)

      RETURN

      END


