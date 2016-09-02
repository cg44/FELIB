
    MODULE mod_matvec

      USE felib_globals,only : wp
      USE mod_space,only : create

      PRIVATE
      PUBLIC matvec

    CONTAINS

      SUBROUTINE matvec(a,v,w,m,n,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      MATVEC post-multiplies the matrix A by the vector V, storing
!      the result in the vector W

! HISTORY

!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 Jul 2000 (CG)

! ARGUMENTS in
!      A       array 
!      V       vector 
!      M       number of rows of A to be used in the
!              multiplication
!      N       number of columns of A and the number of
!              elemenets of V to be used in the multiplication
!      ITEST   error checking option!
! ARGUMENTS out
!      W       vector which contains the result of
!              the operation W=A*V

! ROUTINES called

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments
        INTEGER, OPTIONAL :: m, n, itest
        REAL (wp), POINTER :: a(:,:), v(:), w(:)

! Local variables

        INTEGER :: mm, nn, i, j
        REAL (wp) :: x
        CHARACTER (6) :: srname = 'MATVEC'

! Check optional arguments
        IF (present(m)) THEN
          mm = m
        ELSE
          mm = size(a,1)
        END IF
        IF (present(n)) THEN
          nn = n
        ELSE
          nn = size(a,2)
        END IF

! Get space for W if necessary

        CALL create(w,mm)

!     Main body
        DO i = 1, mm
          x = 0.0D0
          DO j = 1, nn
            x = x + a(i,j)*v(j)
          END DO
          w(i) = x
        END DO

      END SUBROUTINE matvec
    END MODULE mod_matvec
