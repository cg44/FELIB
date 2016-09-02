
    MODULE mod_matran

      USE felib_globals,only : wp
      USE mod_space,only: create

      PRIVATE
      PUBLIC matran

    CONTAINS

      SUBROUTINE matran(a,b,m,n,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      MATRAN forms the transpose of the matrix A

! HISTORY

!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 Jul 2000 (CG)

! ARGUMENTS in
!      A       array to be transposed
!      M       number of rows of A to be transposed
!      N       number of columns of A to be transposed
!      ITEST   error checking option 

! ARGUMENTS out
!      B       array B(J, I) = A(I, J) for
!              I=1(1)M and J=1(1)N

! ROUTINES called

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: m, n, itest
        REAL (wp), POINTER :: a(:,:), b(:,:)

! Local variables

        INTEGER :: i, ierror, j, mm, nn
        CHARACTER (6) :: srname = 'MATRAN'

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

! Create b if necessary

        CALL create(b,nn,mm)

!     Main body

        DO i = 1, mm
          DO j = 1, nn
            b(j,i) = a(i,j)
          END DO
        END DO

      END SUBROUTINE matran
    END MODULE mod_matran
