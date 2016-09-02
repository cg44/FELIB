
    MODULE mod_asful

! .. Use Statements ..
      USE felib_globals, ONLY : wp
      USE mod_errmes
! .. Default Accessibility ..
      PRIVATE
! .. Public Statements ..
      PUBLIC :: asful

    CONTAINS

      SUBROUTINE asful(sysk,elk,steer,dofel,itest)
!-----------------------------------------------------------------------
! PURPOSE
!      ASFUL assembles full real system matrix

! HISTORY

!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  1 April 2002 (CG)

! ARGUMENTS in
!      SYSK    contains system matrix prior to addition of
!              current elemnt matrix contribution
!      ELK     element matrix
!      STEER   contains freedom numbers associated with element
!              matrix contributions to system matrix
!      DOFEL   maximum number of degrees of freedom associated
!              with element type (optional). If ommitted DOFEL=SIZE(steer)
!      ITEST   error checking option (optional). If ommitted ITEST=0

! ARGUMENTS out
!      SYSK    system matrix

! ROUTINES called
!      ERRMES

!***********************************************************************!

! Dummy arguments

! Local variables

! .. Scalar Arguments ..
        INTEGER :: dofel, itest
! .. Array Arguments ..
        REAL (wp) :: elk(:,:), sysk(:,:)
        INTEGER :: steer(:)
! .. Local Scalars ..
        INTEGER :: dofel1, i, ierror, jtest, ktest
        CHARACTER (6) :: srname = 'ASFUL '
! .. Intrinsic Functions ..
        INTRINSIC present, size, ubound
! .. Optional Statements ..
        OPTIONAL :: dofel, itest
! .. Intent Statements ..
        INTENT (IN) :: dofel, elk, steer
        INTENT (INOUT) :: itest, sysk
! Checking optional arguments

        IF (present(itest)) THEN
          ktest = itest
        ELSE
          ktest = 0
        END IF
        IF (present(dofel)) THEN
          dofel1 = dofel
        ELSE
          dofel1 = size(steer)
        END IF

!     Parameter checking

        jtest = ktest
        IF (jtest/=-1) THEN
          ierror = 0
          IF (ubound(steer,1)<dofel1) ierror = 3
          IF (ubound(elk,1)<dofel1 .OR. ubound(elk,2)<dofel1) ierror = 2
          IF (dofel1<0) ierror = 1
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) RETURN
        END IF

!     Main loops

!        DO i = 1, dofel
!          steeri = steer(i)
!          IF (steeri>0) THEN
!            DO j = 1, dofel
!              steerj = steer(j)
!              IF (steerj>0) THEN

!     Range checking on STEERI and STEERJ

!                IF (jtest/=-1) THEN
!                  ierror = 0
!                  IF (ubound(sysk,1)<steeri .OR. ubound(sysk,2)<steerj) &
!                    ierror = 4
!                  itest = errmes(jtest,ierror,srname)
!                  IF (itest/=0) RETURN
!                END IF

!                sysk(steeri,steerj) = sysk(steeri,steerj) + elk(i,j)
!              END IF
!            END DO

!          END IF
!        END DO
        IF (jtest/=-1) THEN
          ierror = 0
          DO i = 1, dofel1
            IF (ubound(sysk,1)<steer(i) .OR. ubound(sysk,2)<steer(i)) &
              ierror = 4
          END DO
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) RETURN
        END IF

        sysk(steer,steer) = sysk(steer,steer) + elk(1:dofel1,1:dofel1)

      END SUBROUTINE asful

    END MODULE mod_asful
