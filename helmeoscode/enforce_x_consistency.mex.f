








C
C YPRIMEG.FOR - Gateway function for YPRIME.FOR
C
C This is an example of the FORTRAN code required for interfacing
C a .MEX file to MATLAB.
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C Copyright 1984-2004 The MathWorks, Inc.
C $Revision: 1.9.2.1 $
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      INTEGER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      INTEGER MXCREATEDOUBLEMATRIX, MXGETPR


C-----------------------------------------------------------------------
C

      INTEGER MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C

cccccccccccc THINGS TO READ IN (pointer)
c      INTEGER YPP, TP, YP
      INTEGER*8 ye_inp_p, xnut_p, xprot_p
     1     , xalfa_p, xh_p, a_p, x_p


cccccccccccc THINGS TO SENT BACK (pointer)
      INTEGER*8 s_xnut_p, s_xprot_p
     1     , s_xalfa_p, s_xh_p, s_a_p, s_x_p
      INTEGER*8 s_xconsistent_p

C-----------------------------------------------------------------------
C

c   NORMAL VALUES OF THINGS
      INTEGER M, N
c      REAL*8 RYPP(4), RTP, RYP(4)
      real*8 r_ye_inp, r_xnut, r_xprot
     1     , r_xalfa, r_xh, r_a, r_x
      integer*4 i_xconsistent

C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 7) THEN
        CALL MEXERRMSGTXT('requires 7 input arguments')
      ELSEIF (NLHS .NE. 7) THEN
        CALL MEXERRMSGTXT('requires 7 output argument')
      ENDIF
C
C CHECK THE DIMENSIONS OF Y.  IT CAN BE 4 X 1 OR 1 X 4.
C
      do iii=1,NRHS

         M = MXGETM(PRHS(iii))
         N = MXGETN(PRHS(iii))
C
         IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
            CALL MEXERRMSGTXT('All required to be a 1 x 1 vector')
         ENDIF

      end do
C
C CREATE A MATRIX FOR RETURN ARGUMENT

c      write(*,*) 'mark1'
C
c      PLHS(1) = MXCREATEDOUBLEMATRIX(1,1,0)
c      PLHS(1) = MXCREATEDOUBLESCALAR(0.0)
      PLHS(1) = 
     1 mxcreatenumericarray(1,1,mxClassIDFromClassName('int32'),0)
      do iii=2,NLHS
         PLHS(iii) = MXCREATEDOUBLESCALAR(0.0)
      end do



C
c      write(*,*) 'mark2'
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
c      YPP = MXGETPR(PLHS(1))
      s_xconsistent_p = MXGETPR(PLHS(1))
      s_xnut_p = MXGETPR(PLHS(2))
      s_xprot_p = MXGETPR(PLHS(3))
      s_xalfa_p = MXGETPR(PLHS(4))
      s_xh_p = MXGETPR(PLHS(5))
      s_a_p = MXGETPR(PLHS(6))
      s_x_p = MXGETPR(PLHS(7))
C
c      write(*,*) 'mark3'
c mxGetScalar
c      TP = MXGETPR(PRHS(1))
c      YP = MXGETPR(PRHS(2))
      ye_inp_p = MXGETPR(PRHS(1))
      xnut_p = MXGETPR(PRHS(2))
      xprot_p = MXGETPR(PRHS(3))
      xalfa_p = MXGETPR(PRHS(4))
      xh_p = MXGETPR(PRHS(5))
      a_p = MXGETPR(PRHS(6))
      x_p = MXGETPR(PRHS(7))
C
c      write(*,*) 'mark4'

C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      CALL MXCOPYPTRTOREAL8(ye_inp_p, r_ye_inp, 1)
c      write(*,*) 'mark4a'
      CALL MXCOPYPTRTOREAL8(xnut_p, r_xnut, 1)
c      write(*,*) 'mark4b'
      CALL MXCOPYPTRTOREAL8(xprot_p, r_xprot, 1)
c      write(*,*) 'mark4c'
      CALL MXCOPYPTRTOREAL8(xalfa_p, r_xalfa, 1)
c      write(*,*) 'mark4d'
      CALL MXCOPYPTRTOREAL8(xh_p, r_xh, 1)
c      write(*,*) 'mark4e'
      CALL MXCOPYPTRTOREAL8(a_p, r_a, 1)
c      write(*,*) 'mark4f'
      CALL MXCOPYPTRTOREAL8(x_p, r_x, 1)
C
c      write(*,*) 'mark5'
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  
C
      call enforce_x_consistency(r_ye_inp, r_xnut, r_xprot, r_xalfa,
     1     r_xh, r_a, r_x, i_xconsistent)
c      CALL YPRIME(RYPP,RTP,RYP)
C
c      write(*,*) 'mark6'
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
c      CALL MXCOPYREAL8TOPTR(RYPP, YPP, 4)
      CALL MXCOPYINTEGER4TOPTR(i_xconsistent, s_xconsistent_p, 1)
      CALL MXCOPYREAL8TOPTR(r_xnut,s_xnut_p , 1)
      CALL MXCOPYREAL8TOPTR(r_xprot,s_xprot_p , 1)
      CALL MXCOPYREAL8TOPTR(r_xalfa,s_xalfa_p , 1)
      CALL MXCOPYREAL8TOPTR(r_xh,s_xh_p , 1)
      CALL MXCOPYREAL8TOPTR(r_a,s_a_p , 1)
      CALL MXCOPYREAL8TOPTR(r_x,s_x_p , 1)

c      write(*,*) 'mark7'
C
      RETURN
      END

