MODULE mod_matadd

      USE felib_globals,only : wp
      USE felib_utils,only : setopt

      PRIVATE
      PUBLIC matadd

CONTAINS

SUBROUTINE matadd(a,b,m,n,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      MATADD adds the matrix A to B, storing the result in A
!
! HISTORY
!
!      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.1  29 Oct 1979 (IMS)
!
! ARGUMENTS in
!      A       array of dimension (IA, JA)
!      B       array of dimension (IB, JB)
!      M       number of rows of A and B to be added
!      N       number of columns of A and B to be added
!      ITEST   error checking option
!
! ARGUMENTS out
!      A       on exit, contains sum of A and B
!
! ROUTINES called
!      ERRMES
!
!***********************************************************************
IMPLICIT NONE

! Dummy arguments

INTEGER, optional :: itest , m , n
REAL(wp) , pointer, DIMENSION(:,:) :: a
REAL(wp) , pointer, DIMENSION(:,:) :: b

INTENT (IN)  m , n
INTENT (INOUT) itest

INTEGER :: mm, nn, jtest


! Optinal parameters

      CALL SETOPT(mm,size(a,1),m)
      CALL SETOPT(nn,size(a,2),n)
      CALL SETOPT(jtest,0,itest)

!     Main body

       a(1:mm,1:nn) = a(1:mm,1:nn) + b(1:mm,1:nn)

END SUBROUTINE matadd
END MODULE mod_matadd
