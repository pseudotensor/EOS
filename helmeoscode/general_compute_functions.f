

c     Assumes xnut, xprot, xalfa, xh are set globally
      subroutine enforce_x_consistency(ye_inp, xnut, xprot, xalfa, xh, a, x, xconsistent)
      implicit none
      save

      include 'eosparms.f'
c     Constants
      include 'const.dek'

c     Local parameters
      logical ifxh,ifxnut,ifxprot,ifxalfa,ifnonelargest
c     Local variables
      double precision xcheckafter
      
c     Passed quantity
      double precision ye_inp, xnut, xprot, xalfa, xh, a, x
      integer xconsistent


      double precision xhorig,xnutorig,xprotorig,xalfaorig

      integer didfixx
      integer iii

ccccccccccccccccccccccccc
c     
c     xcheck
c     
c     consistency of fractions relation #2 in Shen EOS guide.ps
c     
ccccccccccccccccccccccccc

      didfixx=0

      xhorig=xh
      xnutorig=xnut
      xprotorig=xprot
      xalfaorig=xalfa


ccccccccccccccccccccccccc
c     
c     Determine which species is largest
c     
ccccccccccccccccccccccccc
c     xcheck = 1.0 - (xnut + xprot + xalfa + xh)
c     likely that error is mostly in xh since reaches large values and can suddenly drop, so put error there
c     However, in general want to correct largest fractions that will absorb smallest relative error
      ifxh=(((xh>xnut) .AND. (xh>xprot)) .AND. (xh>xalfa) )
      ifxnut=(((xnut>xh) .AND. (xnut>xprot)) .AND. (xnut>xalfa) )
      ifxprot=(((xprot>xh) .AND. (xprot>xnut)) .AND. (xprot>xalfa) )
      ifxalfa=(((xalfa>xh) .AND. (xalfa>xnut)) .AND. (xalfa>xprot) )
      ifnonelargest=((ifxh).AND.(ifxnut)).AND.((ifxprot).AND.(ifxalfa))



ccccccccccccccccccccccccc
c     
c     Correct largest x
c     
ccccccccccccccccccccccccc

c     Make 3 passes
      do iii=1,3

         if(ifxh.OR.ifnonelargest) then
            xh    = 1.0 - (xnut  + xprot + xalfa)
            if(xh.lt.xhtrust .OR. xh.gt.1.0) then
c     If xh got messed up, then try xnut
               xh=xhorig
               ifxh=.FALSE.
               ifnonelargest=.FALSE.
               ifxnut=.TRUE.
            else
               didfixx=1
            end if
         end if

         if(ifxnut) then
            xnut = 1.0 - (xh    + xprot + xalfa)
            if(xnut.lt.xmin .OR. xnut.gt.1.0) then
c     If xnut got messed up, then try xprot
               xnut=xnutorig
               ifxnut=.FALSE.
               ifxprot=.TRUE.
            else
               didfixx=1
            end if
         end if

         if(ifxprot) then
            xprot = 1.0 - (xnut + xh    + xalfa)
            if(xprot.lt.xmin .OR. xprot.gt.1.0) then
c     If xprot got messed up, then try xalfa
               xprot=xprotorig
               ifxprot=.FALSE.
               ifxalfa=.TRUE.
            else
               didfixx=1
            end if
         end if

         if(ifxalfa) then
            xalfa = 1.0 - (xnut + xprot + xh)
            if(xalfa.lt.xmin .OR. xalfa.gt.1.0) then
c     If xprot got messed up, then try xh
               xalfa=xalfaorig
               ifxalfa=.FALSE.
               ifxh=.TRUE.
            else
               didfixx=1
            end if
         end if

         if(didfixx.eq.1) then
            go to 652
         end if

      end do

 652  continue

ccccccccccccccccccccccccc
c     
c     Enforce minimums after correction
c     
ccccccccccccccccccccccccc
      if(xnut<xmin) then
         xnut=xmin
      end if
      if(xprot<xmin) then
         xprot=xmin
      end if
      if(xalfa<xmin) then
         xalfa=xmin
      end if
      if(xh<xhtrust) then
         xh=xhtrust
      end if


c     *0.99999's are to avoid machine issues
c     *10.0 on xchecktolerance is for Matlab
      xcheckafter = 1.0 - (xnut + xprot + xalfa + xh)
      if(
     1     (abs(xcheckafter)>xchecktolerance*10.0)
     1     .OR.(xnut<xmin*0.99999 .OR. xnut>1.0)
     1     .OR.(xprot<xmin*0.99999 .OR. xprot>1.0)
     1     .OR.(xalfa<xmin*0.99999 .OR. xalfa>1.0)
     1     .OR.(xh<xhtrust*0.99999 .OR. xh>1.0)
     1     ) then
         write(*,*) 'Problem with xcheckafter(1)=',xcheckafter,xnut,xprot,xalfa,xh
         write(*,*) 'Problem with xcheckafter(2)=',xmin,xhtrust
         xconsistent=0
      else
         xconsistent=1
      end if


c Some matlab issue with precision:

cBegin mex test Problem with xcheckafter(1)= -1.000000290907188E-010  0.800000000000000
c  0.100000000000000       1.000000013351432E-010  0.100000000000000
c Problem with xcheckafter(2)=  1.000000013351432E-010  1.000000013351432E-010




      return
      end








c     Assumes ye_inp, xnut, xprot, xalfa, xh, x are set globally
      subroutine enforce_ye_consistency(ye_inp, xnut, xprot, xalfa, xh, a, x, yeconsistent)
      implicit none
      save

      include 'eosparms.f'
c     Constants
      include 'const.dek'

c     passed variables
      integer yeconsistent
      double precision ye_inp, xnut, xprot, xalfa, xh, a, x
      

c     Local parameters
c     Local variables
      double precision abarnum,abar,zbar,yelocal
c     double precision yelocal

      logical ifoutputfixtype

      double precision xcheck

      double precision xhorig,xnutorig,xprotorig,xalfaorig
      double precision xorig



      xhorig=xh
      xnutorig=xnut
      xprotorig=xprot
      xalfaorig=xalfa
      xorig=x

      ifoutputfixtype=.FALSE.
c      ifoutputfixtype=.TRUE.
      yeconsistent=0



      if((xnut.gt.xprot).AND.(xnut.gt.xalfa).AND.(xnut.gt.xh)) then
         xnut = -ye_inp + 1.0 - 0.5*xalfa + xh*(x-1.0)
         xprot = -0.5*xalfa-x*xh+ye_inp
         if(xprot.lt.xmin .OR. xprot.gt.1.0 .OR. xnut.lt.xmin .OR. xnut.gt.1.0 ) then
            xnut = xnutorig
            xprot = xprotorig
            xnut = (0.5*xalfa+xprot-x*(xprot+xalfa-1.0)-ye_inp)/x
            xh = (ye_inp-xprot-0.5*xalfa)/x
c     If making correction on x, then fix a as well so Y_e calculation consistent
            if(a.lt.aheavtrust) then
               a = 1.0
            end if
            if(xnut.lt.xmin .OR. xnut.gt.1.0 .OR. xh.lt.xhtrust .OR. xh.gt.1.0 ) then
               write(*,*) 'Fix did not work 1b',xnut,xh
               xnut = xnutorig
               xh = xhorig
            else
               if(ifoutputfixtype) write(*,*) 'type1bfix',xnut,xh
               yeconsistent=1
            end if
         else
            if(ifoutputfixtype) write(*,*) 'type1afix',xnut,xprot
            yeconsistent=1
         end if
      else if((xprot.gt.xnut).AND.(xprot.gt.xalfa).AND.(xprot.gt.xh)) then
         xprot = ye_inp - 0.5*xalfa - xh*x
         xnut = 1.0 - (xprot+xalfa+xh)
         if(xprot.lt.xmin .OR. xprot.gt.1.0 .OR. xnut.lt.xmin .OR. xnut.gt.1.0 ) then
c     Then assume couldn't fix
            write(*,*) 'Fix did not work 2',xprot,xnut
            xprot = xprotorig
            xnut = xnutorig
         else
            if(ifoutputfixtype) write(*,*) 'type2fix',xprot,xnut
            yeconsistent=1
         end if
      else if((xh.gt.xnut).AND.(xh.gt.xalfa).AND.(xh.gt.xprot)) then
         x = (xalfa+2.0*xprot-2.0*ye_inp)/(2.0*(-1.0+xalfa+xnut+xprot))
         xh = 1.0 - (xprot+xalfa+xnut)
c     If making correction on x, then fix a as well so Y_e calculation consistent
         if(a.lt.aheavtrust) then
            a = 1.0
         end if
         if(xh.lt.xhtrust .OR. xh.gt.1.0 .OR. x.lt.zheavtrust/aheavtrust) then
            x = xorig
            xh = xhorig
            xh = -(0.5*xalfa+xprot-ye_inp)/x
            xnut = (0.5*xalfa+xprot-x*(-1.0+xalfa+xprot)-ye_inp)/x
            if(xh.lt.xhtrust .OR. xh.gt.1.0 .OR. xnut.lt.xmin .OR. xnut.gt.1.0 ) then
               write(*,*) 'Fix did not work 3',xh,xnut
               xh = xhorig
               xnut = xnutorig
            else
               if(ifoutputfixtype) write(*,*) 'type3bfix',x,xh
               yeconsistent=1
            end if
         else
            if(ifoutputfixtype) write(*,*) 'type3afix',x,xh
            yeconsistent=1
         end if
      else if((xalfa.gt.xnut).AND.(xalfa.gt.xprot).AND.(xalfa.gt.xh)) then
         xalfa = -2.0*(x*xh+xprot-ye_inp)
         xnut = 1.0 + (2.0*x-1)*xh+xprot-2.0*ye_inp
         if(xalfa.lt.xmin .OR. xalfa.gt.1.0 .OR. xnut.lt.xmin  .OR. xnut.gt.1.0 ) then
            xalfa = xalfaorig
            xnut = xnutorig
            xalfa = 2.0*(x-1.0)*xh-2.0*(xnut+ye_inp-1.0)
            xprot = -1+xh-2.0*x*xh+xnut+2.0*ye_inp
            if(xalfa.lt.xmin .OR. xalfa.gt.1.0 .OR. xprot.lt.xmin .OR. xprot.gt.1.0 ) then
c     Then assume couldn't fix
               write(*,*) 'Fix did not work 4b',xalfa,xprot
               xalfa = xalfaorig
               xprot = xprotorig
            else
               if(ifoutputfixtype) write(*,*) 'type4bfix',xalfa,xprot
               yeconsistent=1
            end if
         else
            if(ifoutputfixtype) write(*,*) 'type4afix',xalfa,xnut
            yeconsistent=1
         end if
      else
         write(*,*) 'yebad',xnut,xprot,xalfa,xh,x
      end if





      if(yeconsistent.eq.0) then
         write(*,*) 'Did not make Y_e consistent',ye_inp,xnut,xprot,xalfa,xh,x,a
      end if

      
c     Make sure x's are consistently kept satisfying below
      xcheck = 1.0 - (xnut + xprot + xalfa + xh)

      if(abs(xcheck)>xchecktolerance) then
         write(*,*) 'Thought fixed Y_e, but messed up xcheck in process',ye_inp,xnut,xprot,xalfa,xh,x
      end if



c Matlab problem:
c Fix did not work 1  0.705138549764483       6.759689526642490E-006
c Did not make Y_e consistent  6.744462874835755E-002  0.705138549764483
c  6.759689526642490E-006  5.802894130654267E-005  0.294796661604684
c  0.228687623657101
cSTILL problem with Y_e(1): (yeconsistent=0, errorye=5.55493e-05): 167 10 146
cSTILL problem with Y_e(2): 0.0674446 0.705139 6.75969e-06 5.80289e-05 0.294797 152.578 0.228688
cSTILL problem with Y_e(3): 1.41424 1.41424 0.0953937 0.0674521







      return
      end








cccccccccccccccccccccccc
c     
c     Compute Y_e (and abar,abarnum,zbar) for "inputs" of:
c     xnut,xprot,xalfa,xh,a,x
c     
ccccccccccccccccccccccc
      subroutine compute_nuclear_azbar(xnut, xprot, xalfa, xh, a, x,
     1     abarnum,abar,zbar
     1     ,yelocal)
      implicit none


c..   local variables
      double precision ytot1,zbarxx
      double precision  aheavmin,xheavmin

c     Passed
      double precision xnut, xprot, xalfa, xh, a, x

c     Variables to be returned
      double precision abarnum,abar,zbar
      double precision yelocal



c     include 'eos_m4c.commononly.inc'
c     JCM: below used to store electron EOS properties needed by LSEOS to generate final solution
c     include 'el_eos.inc'
c     JCM:
      include 'const.dek'
      include 'eosparms.f'
c     include 'vector_eos.dek'
c     include 'vector_sneos.dek'



c..   abar and zbar of the whole mixture
c     Since LSEOS assumes nothing beyond alphas, then this overwrites choices made at first
c..   JCM: Therefore whicheleeos=10,12 don't make sense since don't know what abar,zbar will be until here
      ytot1   = xnut  + xprot + 0.25d0*xalfa 
c     if (a .ne. 0.0) ytot1 = ytot1 + xh/a
c     Assume that below this one really means aheav=0 even if xheav very small
c     Just artifact of doing log-interpolation on a and x in Matlab
c     aheavtrust
c     xhtrust
      if ((a .ge. aheavtrust).AND.(x*a.ge.zheavtrust)) ytot1 = ytot1 + xh/a
      zbarxx  = xprot + 0.5d0*xalfa
      if ((a .ge. aheavtrust).AND.(x*a.ge.zheavtrust)) zbarxx = zbarxx + x*xh

      abarnum    = 1.0d0/ytot1
      abar=abarnum
      zbar    = zbarxx * abar
      yelocal = zbar/abarnum



      return
      end







