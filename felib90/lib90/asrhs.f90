    MODULE mod_asrhs

      USE felib_globals, ONLY : wp
      USE mod_errmes

      PRIVATE
      PUBLIC asrhs

    CONTAINS

      SUBROUTINE asrhs(rhs,value,steer,dofel,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      ASRHS the routine adds into the right-hand side of a system
!      the values contianed in an element vector, thus
!      forming the right-hand side.

! HISTORY

!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   1 Apr 2002 (CG)

! ARGUMENTS in
!      RHS     the right-hand side of the system
!      VALUE   the element vector of the current element to
!              be added into the right-hand side
!      STEER   the steering vector containing the freedom
!              numbers of the freedoms associated with the
!              current element in the local order
!      DOFEL   the maximum number of degrees of freedom on
!              an element of the current type
!      ITEST   error checking option (optional). If omitted ITEST=0

! ROUTINES called


!***********************************************************************


        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: dofel, itest
        INTEGER :: steer(:)
        REAL (wp) :: rhs(:), value(:)
        INTENT (IN) :: value, steer, dofel
        INTENT (INOUT) :: rhs, itest

! Local variables

        INTEGER :: ierror, jtest, k, irhs, ivalue, isteer, dofel1, &
          ktest
        CHARACTER (6) :: srname = 'ASRHS'

! Check optional arguments

        IF (present(dofel)) THEN
          dofel1 = dofel
        ELSE
          dofel1 = size(steer)
        END IF
        IF (present(itest)) THEN
          ktest = itest
        ELSE
          ktest = 0
        END IF

! Size of rhs, value and steer

        irhs = size(rhs)
        ivalue = size(value)
        isteer = size(steer)

! Check on arguments

        jtest = ktest
        IF (jtest/=-1) THEN
          ierror = 0
          IF (isteer<dofel1) ierror = 3
          IF (ivalue<dofel1) ierror = 2
          IF (dofel1<=0) ierror = 1
          ktest = errmes(jtest,ierror,srname)
        END IF


! Check values is steer

         IF (jtest/=-1) THEN
           ierror=0
           DO k = 1, dofel1
              IF(steer(k)>size(rhs) .OR. steer(k)<=0) ierror=4
           END DO
           ktest=errmes(jtest,ierror,srname)
           IF(ktest/=0) RETURN
         END IF
        
         rhs(steer(1:dofel1)) = rhs(steer(1:dofel1)) + value(1:dofel1)

      END SUBROUTINE asrhs
    END MODULE mod_asrhs
