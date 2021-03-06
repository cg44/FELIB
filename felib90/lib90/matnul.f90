
    MODULE mod_matnul

!-----------------------------------------------------------------------
! PURPOSE
!      MATNUL creates and sets matrix A to the null matrix

! HISTORY

!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 Jul 2000 (CG)

! ARGUMENTS in
!      M       number of rows of A set to zero
!      N       number of columns of A set to zero
!      ITEST   error checking option 

! ARGUMENTS out
!      A       array set to zeros

! ROUTINES called

!***********************************************************************

      USE felib_globals,only : wp
      USE mod_space,only : create
      USE mod_setopt,only : setopt

      PRIVATE
      PUBLIC matnul

      INTERFACE matnul
        MODULE PROCEDURE complex_matnul
        MODULE PROCEDURE real_matnul
        MODULE PROCEDURE integer_matnul
      END INTERFACE

      CHARACTER (5) :: srname = 'MATNUL'    

    CONTAINS

       SUBROUTINE complex_matnul(a,m,n,itest)

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: m, n, itest
        COMPLEX (wp), POINTER :: a(:,:)

        INTENT (IN) :: m,n
        INTENT (INOUT) :: itest

! Local variables

        INTEGER :: i, ierror, j, mm, nn
        CHARACTER (5) :: srname = 'MATNUL'

! Check optional arguments and association of A

        CALL setopt(mm,size(a,1),m)
        CALL setopt(nn,size(a,2),n)

! Create A if necessary then initialise

        CALL create(a,mm,nn)

        a(1:mm,1:nn)=(0,0)

      END SUBROUTINE complex_matnul
     
       SUBROUTINE real_matnul(a,m,n,itest)

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: m, n, itest
        REAL (wp), POINTER :: a(:,:)

        INTENT (IN) :: m,n
        INTENT (INOUT) :: itest

! Local variables

        INTEGER :: i, ierror, j, mm, nn
        CHARACTER (5) :: srname = 'MATNUL'

! Check optional arguments and association of A

        CALL setopt(mm,size(a,1),m)
        CALL setopt(nn,size(a,2),n)

! Create A if necessary then initialise

        CALL create(a,mm,nn)

        a(1:mm,1:nn)=0

      END SUBROUTINE real_matnul

       SUBROUTINE integer_matnul(a,m,n,itest)

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: m, n, itest
        INTEGER, POINTER :: a(:,:)

        INTENT (IN) :: m,n
        INTENT (INOUT) :: itest

! Local variables

        INTEGER :: i, ierror, j, mm, nn

! Check optional arguments and association of A

        CALL setopt(mm,size(a,1),m)
        CALL setopt(nn,size(a,2),n)

! Create A if necessary then initialise

        CALL create(a,mm,nn)

        a(1:mm,1:nn)=0

      END SUBROUTINE integer_matnul

    END MODULE mod_matnul
