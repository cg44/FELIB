
    MODULE mod_vecnul

!-----------------------------------------------------------------------
! PURPOSE
!      VECNUL creates and sets vector A to the null vector

! HISTORY

!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 Jul 2000 (CG)

! ARGUMENTS in
!      M       number of rows of A set to zero
!      ITEST   error checking option 

! ARGUMENTS out
!      A       array set to zeros

! ROUTINES called

!***********************************************************************

      USE felib_globals,only : wp
      USE mod_space,only : create

PRIVATE
PUBLIC vecnul

      INTERFACE vecnul
        MODULE PROCEDURE real_vecnul
        MODULE PROCEDURE integer_vecnul
      END INTERFACE

      CHARACTER (6) :: srname = 'VECNUL'


    CONTAINS

      SUBROUTINE real_vecnul(a,m,itest)

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: m, itest
        REAL (wp), POINTER :: a(:)

        INTENT (IN) :: m

! Local variables

        INTEGER :: i, mm

! Check optional arguments and association of A

        IF (present(m)) THEN
          mm = m
        ELSE
          mm = size(a,1)
        END IF

! Create A if necessary then initialise

        CALL create(a,mm)

        a(1:mm) = 0

      END SUBROUTINE real_vecnul

      SUBROUTINE integer_vecnul(a,m,itest)

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: m, itest
        INTEGER, POINTER :: a(:)

        INTENT (IN) :: m

! Local variables

        INTEGER :: i, mm

! Check optional arguments and association of A

        IF (present(m)) THEN
          mm = m
        ELSE
          mm = size(a,1)
        END IF

! Create A if necessary then initialise

        CALL create(a,mm)

        a(1:mm) = 0

      END SUBROUTINE integer_vecnul

    END MODULE mod_vecnul
