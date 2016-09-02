
    MODULE mod_matmul

      USE felib_globals,only : wp
      USE mod_space,only : create

      PRIVATE
      PUBLIC matmul

    CONTAINS

      SUBROUTINE matmul(a,b,c,l,m,n,itest)

        USE felib_globals

!-----------------------------------------------------------------------
! PURPOSE
!      MATMUL pre-multiplies matrix B by A, storing the result in C

! HISTORY

!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 July 2000 (CG)

! ARGUMENTS in
!      A       real array
!      B       real array
!      L       number of rows of A to be used in multiplication
!      M       number of columns of A and number of rows of B
!              to be used in multiplication
!      N       number of columns of B to be used in
!              multiplication
!      ITEST   error checking option

! ARGUMENTS out
!      C       contains result of matrix multiplication (C=A*B)

! ROUTINES called

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments
        INTEGER, OPTIONAL :: l, m, n, itest
        REAL (wp), POINTER :: a(:,:), b(:,:), c(:,:)
       
        INTENT (IN) :: l,m,n
        INTENT (INOUT) itest

! Local variables

        INTEGER :: i, ierror, j, k, mm, nn, ll
        REAL (wp) :: x

        CHARACTER (6) :: srname = 'MATMUL'

! Check optional arguments

        IF (present(l)) THEN
          ll = l
        ELSE
          ll = size(a,1)
        END IF
        IF (present(m)) THEN
          mm = m
        ELSE
          mm = size(b,1)
        END IF
        IF (present(n)) THEN
          nn = n
        ELSE
          nn = size(b,2)
        END IF

! Allocate space for C if necessary

        CALL create(c,ll,nn)

!     Main body

        DO i = 1, ll
          DO j = 1, nn
            x = 0.0
            DO k = 1, mm
              x = x + a(i,k)*b(k,j)
            END DO
            c(i,j) = x
          END DO
        END DO

      END SUBROUTINE matmul
    END MODULE mod_matmul
