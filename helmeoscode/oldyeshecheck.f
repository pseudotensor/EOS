

c     Assumes ye_inp, xnut, xprot, xalfa, xh, a, x are set globally
      subroutine enforce_ye_consistency()
      implicit none
      save

c     Below includes xnut,xprot, etc. used for ye computation
      include 'eos_m4c.inc'
c     contains single globals (also has abar,abarnum,zbar)
c      include 'vector_eos.single.dek'
      include 'vector_sneos.dek'
c     contains global parameters and Y_e calculation limits
      include 'vector_eos.dek'
c     Constants
      include 'const.dek'

c     Local parameters
      logical ifaheavlarger,ifxh,ifxneut,ifxprot,ifxalfa,ifnonelargest
      logical doxcheck,doyecheck,dothermo1check,dothermo2check
c     Local variables
      double precision abarnum,abar,zbar,yelocal
c      double precision yelocal

c     Logicals
      logical trustheavy,xhdominates,changeaheav,changezheav,changexprot,changexalfa,changexh
      logical changeazheav

      double precision heavyterm


      double precision newxprot,newxalfa,newxh





         if(0) then




c     again, least trustable quantities are related to interpolations of
c     heavy nucleus quantities.  Already fixed xh, now choose to put error
c     equally in aheav and zheav since similar type of quantities
c     for now try only putting into aheav
c     see mathematica Sheneos.nb

c     below returns back good (interpolated) yp
c     correct either a or z, whatever is largest and put minimum on answer

c     First make sure aheav and zheav are reasonable
            if(shenaheav_sing<aheavtrust) then
               shenaheav_sing=1.0
               shenzheav_sing=1.0
            end if

c     Below line (trustheavy) in Matlab is near Y_e calculation, but here Y_e is computed using a function that doesn't return this, so need to compute it
c     trustheavy=(shenxh_sing>xhtrust).AND.(shenaheav_sing>aheavtrust)
            trustheavy=(shenzheav_sing>zheavtrust).AND.(shenaheav_sing>aheavtrust)
            xhdominates=(ifxh.OR.ifnonelargest)
     1           .OR.(ifxneut.AND.(shenxh_sing>shenxprot_sing).AND.(shenxh_sing>shenxalfa_sing))
            ifaheavlarger=(shenaheav_sing>shenzheav_sing)
c     Change aheav if trust and large, or if don't trust
            changeaheav = xhdominates.AND.(( ifaheavlarger.eq..TRUE.).AND.trustheavy).OR.(shenaheav_sing<aheavtrust)
            changezheav = (changeaheav.eq..FALSE.).AND.xhdominates.AND.(trustheavy.OR.(shenzheav_sing<zheavtrust))
            changeazheav=(changeaheav.OR.changezheav)

c     if( (trustheavy.eq..FALSE.).AND.(xhdominates.eq..TRUE.)) then
c     write(*,*) 
c     1           'Bad situation where Xh dominates but can not trust',trustheavy,xhdominates
c     end if

            if(changeaheav.eq..TRUE.) then
               shenaheav_sing = (shenxh_sing*shenzheav_sing)/(yetrue-0.5*shenxalfa_sing-shenxprot_sing)
            end if
            if(changezheav.eq..TRUE.) then
               shenzheav_sing = shenaheav_sing*(yetrue-0.5*shenxalfa_sing-shenxprot_sing)/shenxh_sing
            end if

c     xhange xprot if not chosen by xcheck, xh didn't dominate (so could change aheav or zheav), but larger than one of the other species
c     Figure out what x is second (i.e. between charged species)
            if(changeazheav.eq..FALSE.) then
c     Assume Xh didn't dominate if didn't change, since should have trusted if dominated
c     So must consider changing xprot, xalfa, and xh
c     If Xh dominates, then already changed Xh for xcheck, so can only change xprot or xalfa




               if(trustheavy) then
                  heavyterm = shenxh_sing*shenzheav_sing/shenaheav_sing
               else
                  heavyterm = 0.0
               end if
               newxprot = yetrue - 0.5*shenxalfa_sing - heavyterm

               
               if(trustheavy) then
                  heavyterm = shenxh_sing*shenzheav_sing/shenaheav_sing
               else
                  heavyterm = 0.0
               end if
               newxalfa = 2.0*(yetrue - shenxprot_sing - heavyterm)


               if(trustheavy) then
                  newxh = shenaheav_sing*(yetrue-0.5*shenxalfa_sing-shenxprot_sing)/shenzheav_sing
               else
c     heavyterm = 0.0
                  write(*,*)'Bad situation where
     1                 changing Xh but trustheavy is false'
               end if



               changexprot=(ifxprot.eq..FALSE.
     1              .AND.(shenxprot_sing>shenxh_sing .OR. ifxh .OR. (newxh<xmin) )
     1              .AND.(shenxprot_sing>shenxalfa_sing .OR. ifxalfa .OR. (newxalfa<xmin) )
     1              )

               changexalfa=(ifxalfa.eq..FALSE. .AND. changexprot.eq..FALSE.
     1              .AND.(shenxalfa_sing>shenxh_sing .OR. ifxh .OR. (newxh<xmin) )
     1              .AND.(shenxalfa_sing>shenxprot_sing .OR. ifxprot .OR. (newxprot<xmin) )
     1              )

               changexh=(ifxh.eq..FALSE. .AND. changexprot.eq..FALSE. .AND. changexalfa.eq..FALSE.
     1              .AND.(shenxh_sing>shenxprot_sing .OR. ifxprot .OR. (newxprot<xmin) )
     1              .AND.(shenxh_sing>shenxalfa_sing .OR. ifxalfa .OR. (newxalfa<xmin) )
     1              )

               if(changexprot.eq..TRUE.) then
                  shenxprot_sing = newxprot
               end if
               if(changexalfa.eq..TRUE.) then
                  shenxalfa_sing = newxalfa
               end if
               if(changexh.eq..TRUE.) then
                  shenxh_sing = newxh
               end if

               
            else
               changexprot=.FALSE.

               changexalfa=.FALSE.

               changexh=.FALSE.
               
            end if
            




c     Need to recorrect xcheck
            if((changexprot.eq..TRUE.).OR.(changexalfa.eq..TRUE.).OR.(changexh.eq..TRUE.)) then
               if(ifxh.OR.ifnonelargest) then
                  shenxh_sing    = 1.0 - (shenxneut_sing  + shenxprot_sing + shenxalfa_sing)
               else if(ifxneut) then
                  shenxneut_sing = 1.0 - (shenxh_sing    + shenxprot_sing + shenxalfa_sing)
               else if(ifxprot) then
                  shenxprot_sing = 1.0 - (shenxneut_sing + shenxh_sing    + shenxalfa_sing)
               else if(ifxalfa) then
                  shenxalfa_sing = 1.0 - (shenxneut_sing + shenxprot_sing + shenxh_sing)
               end if
            end if



c     End if old Y_e check method
         end if




         if(1) then

            








         end if


