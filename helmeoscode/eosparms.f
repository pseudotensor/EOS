
c
ccccccccccccccccccccccccccccccc
c
c Determine which basic electron and nuclear EOSs to use
c
c  whichnucleareos:
c
c  0 = HELM EOS
c  1 = LSEOS (and then internally use whicheleeos)
c  2 = TIMMES EOS
c  3 = SHEN EOS
c  4 = KAZ EOS


c  Choose electron (and radiation) EOS ONLY if whichnucleareos=1, otherwise HELM EOS is used for electron EOS
c
c  whicheleeos:
c
c  0 = HELM EOS internally called (limited to tablulated range of density/temperature)
c  1 = LSEOS's original ele EOS (internally called)
c  2 = TIMMES EOS internally called
c  3 = KAZ EOS internally called
c
c Shen EOS assumes whicheleeos is either 0 or 2
c

c..Sets which electron EOS to use
      integer whicheleeos,whichnucleareos,OUTPUTDETAILS
c whicheleeos,whichnucleareos,
      parameter (whicheleeos = 2)
c GODMARK: Below normally "3"
      parameter (whichnucleareos = 3)
c     whether to output details (alot of info, so only for debug usually)
      parameter (OUTPUTDETAILS = 0)

c     Which type of output to generate (see jon_helm_outputstyle.f)
c      0=pure HELM, 1=Kaz-like parts, 2=PWF99 all 3 = HELM+Kaz entropy(causes HARM to fail alot)
c      always use HELM+KAZ2002 for neutrino rates
      integer OUTPUTTYPE
      parameter (OUTPUTTYPE = 0)


c     Whether to use Kaz-based neutrino calculation even with HELM/TIMMES/etc. (i.e. non-Kaz codes without neutrino physics)
      integer kazlikeneutrinos
      parameter (kazlikeneutrinos = 1)




c.. Convergence criterion for LSEOS
      real*8 LSTOL
      parameter (LSTOL=5E-2)



c	Parameters for Shen EOS and general code for determing Y_e
      double precision azmin,xmin,aheavtrust,zheavtrust,xhtrust
     1 	,yetolerance,xchecktolerance
c      common /floorsye/ azmin,xmin,aheavtrust,zheavtrust,xhtrust
      parameter (azmin = 1E-10
     1   	,xmin=1E-10
     1          ,aheavtrust=1.0
     1          ,zheavtrust=azmin
     1          ,xhtrust=xmin
     1	        ,yetolerance=1E-10
     1	        ,xchecktolerance=1E-10
     1 	        )




  
c     local parameter
c     whichtable=0 : Shen original table
c     whichtable=1 : Matlab interpolated Shen table
      integer whichtable
      parameter (whichtable = 1)




      double precision lYpfloor,ltempfloor,Ypfloor,tempfloor,
     1  nbmin,mstarmin,fakefloor
c      common /floors1/ lYpfloor,ltempfloor,Ypfloor,tempfloor,
c     1	nbmin,mstarmin,fakefloor



cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Same values used in matlab's shen_interp.m
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      parameter lYpfloor=lypminin-7
      parameter (lYpfloor=-9.0d0)
c     T[MeV]=10^(-9) MeV is used instead of T= 0 MeV (i.e. far from minimum
c     T in sheneos.t00 table
c      parameter ltempfloor=ltkminin-7
      parameter (ltempfloor=-8.0d0)
      
c      parameter (Ypfloor=10.0d0**(lYpfloor))
c      parameter (tempfloor=10.0d0**(ltempfloor))
c     Below used for MEX (g77) use
      parameter (Ypfloor=1.0D-9)
      parameter (tempfloor=1.0D-8)
 
      parameter (nbmin=1.0D-30)
      parameter (mstarmin=1.0D-10)

      parameter (fakefloor=1.0D-30)







