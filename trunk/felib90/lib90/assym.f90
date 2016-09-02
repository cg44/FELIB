
    MODULE mod_assym

      USE felib_globals, ONLY : wp
      USE mod_errmes

      PRIVATE
      PUBLIC assym

    CONTAINS

      SUBROUTINE assym(sysk,elk,steer,hband,dofel,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      ASSYM adds the contribution from an element matrix to
!      a real symmetric system matrix

! HISTORY

!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   2 July 2002 (CG)

! ARGUMENTS in
!      SYSK    contains system matrix prior to addition of
!              current element matrix contribution
!      ELK     element matrix
!      STEER   contains freedom numbers associated with element
!              matrix contributions to system matrix
!      HBAND   semi-bandwidth of system matrix, including
!              diagonal (optional). If ommitted hband=size(sysk,2).
!      DOFEL   maximum degrees of freedom associated with this
!              element in the mesh (optintal). If ommitted dodel=size(steer)

! ARGUMENTS out
!      SYSK    system matrix

! ROUTINES called

!***********************************************************************


        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: dofel, hband, itest
        INTEGER , POINTER :: steer(:)
        REAL (wp), POINTER :: elk(:,:), sysk(:,:)

        INTENT (IN)  hband, dofel
        INTENT (INOUT) itest

! Local variables

        INTEGER :: i, ierror, j, jtest, steeri, steerj, offset
        INTEGER :: dofel1, hband1, ktest, isteer, ielk, jelk, jsysk, isysk
        CHARACTER (6) :: srname = 'ASSYM'

!     Check optional arguments

        IF (present(hband)) THEN
          hband1 = hband
        ELSE
          hband1 = size(sysk,2)
        END IF
        IF (present(dofel)) THEN
          dofel1 = dofel
        ELSE
          dofel1 = size(steer,1)
        END IF

        IF (present(itest)) THEN
          ktest = itest
        ELSE
          ktest = 0
        END IF

! Set sizes of arrays for checking

        isteer = size(steer)
        ielk = size(elk,1)
        jelk = size(elk,2)
        isysk = size(sysk,1)
        jsysk = size(sysk,2)


!     Parameter checking

        jtest = ktest
        IF (jtest/=-1) THEN
          ierror = 0
          IF (isteer<dofel1) ierror = 4
          IF (ielk<dofel1 .OR. jelk<dofel1) ierror = 3
          IF (jsysk<hband1) ierror = 2
          IF (hband1<=0 .OR. dofel1<=0) ierror = 1
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) RETURN
        END IF

!     Main loops

        DO i = 1, dofel1
          steeri = steer(i)
          IF (steeri>0) THEN
            DO j = 1, dofel1
              steerj = steer(j)
              IF (steerj>0) THEN
                offset = steerj - steeri + hband1
                IF (offset<=hband1) THEN

!     Range checking on STEERI

                  IF (jtest/=-1) THEN
                    ierror = 0
                    IF (isysk<steeri) ierror = 4
                    ktest = errmes(jtest,ierror,srname)
                    IF (ktest/=0) RETURN
                  END IF

                 sysk(steeri,offset) = sysk(steeri,offset) + elk(i,j)
                END IF
              END IF
            END DO
          END IF
        END DO

      END SUBROUTINE assym
    END MODULE mod_assym
