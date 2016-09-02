
    MODULE mod_matinv

      USE felib_globals,only : wp
      USE mod_errmes,only : errmes
      USE mod_space,only : create

      PRIVATE
      PUBLIC matinv

    CONTAINS

      SUBROUTINE matinv(a,b,det,n,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      MATINV forms the inverse of matrix A in B

! HISTORY

!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 July 2000 (CG)

! ARGUMENTS in
!      A       array to be inverted

! ARGUMENTS out
!      B       array containing inverse of A
!      DET     contains value of determinant of matrix A

! ROUTINES called

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: n, itest
        REAL (wp), POINTER :: a(:,:), b(:,:)
        REAL (wp) :: det

        INTENT (IN) :: n
        INTENT (OUT) :: det
        INTENT (INOUT) :: itest

! Local variables

        INTEGER :: ln, m, k, l, ktest, jtest
        REAL (wp) :: x
        CHARACTER (6) :: srname = 'MATINV'

        INTRINSIC abs, present

! Check optional arguments

        IF (present(itest)) THEN
          ktest = itest
        ELSE
          ktest = 0
        END IF
        IF (present(n)) THEN
          ln = n
        ELSE
          ln = size(a,1)
        END IF

        m = size(a,2)

! Get space for inverse

        CALL create(b,ln,ln)

! Main body

        m = ln - 1
        IF (m==2) THEN

!     Code for 3x3 matrix

          det = a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
          det = det - a(1,2)*(a(2,1)*a(3,3)-a(3,1)*a(2,3))
          det = det + a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))

          b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
          b(2,1) = -a(2,1)*a(3,3) + a(3,1)*a(2,3)
          b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
          b(1,2) = -a(1,2)*a(3,3) + a(3,2)*a(1,3)
          b(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
          b(3,2) = -a(1,1)*a(3,2) + a(3,1)*a(1,2)
          b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
          b(2,3) = -a(1,1)*a(2,3) + a(2,1)*a(1,3)
          b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)

          DO k = 1, 3
            DO l = 1, 3
              b(k,l) = b(k,l)/det
            END DO
          END DO

        ELSE IF (m==1) THEN

!     Code for 2x2 matrix

          det = a(1,1)*a(2,2) - a(1,2)*a(2,1)

!     Check that determinant is not near zero

!IF ( jtest/=-1 ) THEN
!   ierror = 0
!   IF ( abs(det)<unflo() ) ierror = 4
!   itest = errmes(jtest,ierror,srname)
!   IF ( itest/=0 ) RETURN
!END IF

          b(1,1) = a(2,2)
          b(1,2) = -a(1,2)
          b(2,1) = -a(2,1)
          b(2,2) = a(1,1)
          DO k = 1, 2
            DO l = 1, 2
              b(k,l) = b(k,l)/det
            END DO
          END DO
        END IF
        RETURN

      END SUBROUTINE matinv
    END MODULE mod_matinv
