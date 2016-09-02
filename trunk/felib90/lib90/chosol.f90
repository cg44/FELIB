
    MODULE mod_chosol

      USE felib_globals, ONLY : wp
      USE mod_errmes

      PRIVATE
      PUBLIC chosol

    CONTAINS

      SUBROUTINE chosol(a,r,n,hband,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      CHOSOL solves a set of real symmetric positive definite banded
!      equations with a single right hand side by choleski decomposition.
!      Only the lower band and diagonal are stored in a rectangular
!      array A.

! HISTORY

!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 July 2002 (CG)

! ARGUMENTS in
!      A       on entry contains lower half of pd symmetric
!              band matrix stored as a rectangular array
!      R       contains elements of right hand side
!      N       order of matrix A
!      HBAND   semi-bandwidth of A (includes diagonal)
!      ITEST

! ARGUMENTS out
!      A       on exit, contains lower triangular reduced
!      R       matrix solution vector

! ROUTINES called
!      ERRMES

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: n, hband, itest
        REAL (wp),POINTER :: a(:,:), r(:)

        INTENT (IN) :: n, hband
        INTENT (INOUT) :: itest

! Local variables

        INTEGER :: i, ierror, ij, ik, j, jtest, k, l, la, lb, lk, w, ln, &
          ktest, ir, ia, ja, lhband
        REAL (wp) :: x

        CHARACTER (6) :: srname = 'CHOSOL'

        INTRINSIC sqrt

! Checking optional arguments

        IF (present(n)) THEN
          ln = n
        ELSE
          ln = size(a,1)
        END IF

        IF (present(hband)) THEN
          lhband = hband
        ELSE
          lhband = size(a,2)
        END IF

        IF (present(itest)) THEN
          ktest=itest
        ELSE
          ktest=0
        END IF
        

! Collecting array sizes

        ir = size(r,1)
        ia = size(a,1)
        ja = size(a,2)

!     Parameter checking

        jtest = ktest
        IF (jtest/=-1) THEN
          ierror = 0
          IF (ir<ln) ierror = 3
          IF (ia<ln .OR. ja<lhband) ierror = 2
          IF (ln<=0 .OR. lhband<=0) ierror = 1
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) RETURN
        END IF

!     Main body

        w = lhband - 1
        DO i = 1, ln
          x = 0.0D0
          DO j = 1, w
            x = x + a(i,j)*a(i,j)
          END DO

          a(i,w+1) = sqrt(a(i,w+1)-x)

          DO k = 1, w
            x = 0.0D0
            IF (i+k<=ln) THEN
              IF (k/=w) THEN
                l = w - k
                DO
                  ik = i + k
                  lk = l + k
                  x = x + a(ik,l)*a(i,lk)
                  l = l - 1
                  IF (l == 0) EXIT
                END DO
              END IF
              la = i + k
              lb = w - k + 1
              a(la,lb) = (a(la,lb)-x)/a(i,w+1)
            END IF
          END DO
        END DO
        r(1) = r(1)/a(1,w+1)
        DO i = 2, ln
          x = 0.0D0
          k = 1
          IF (i<=w+1) k = w - i + 2
          DO j = k, w
            ij = i + j - w - 1
            x = x + a(i,j)*r(ij)
          END DO
          r(i) = (r(i)-x)/a(i,w+1)
        END DO
        r(ln) = r(ln)/a(ln,w+1)
        i = ln - 1
        DO
          x = 0.0D0
          l = i + w
          IF (i>ln-w) l = ln
          lhband = i + 1
          DO j = lhband, l
            ij = w + i - j + 1
            x = x + a(j,ij)*r(j)
          END DO
          r(i) = (r(i)-x)/a(i,w+1)
          i = i - 1
        IF (i == 0) EXIT
        END DO
      END SUBROUTINE chosol

    END MODULE mod_chosol
