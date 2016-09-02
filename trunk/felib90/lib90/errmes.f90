    MODULE mod_errmes

      USE mod_adunit
      USE mod_erunit

      PRIVATE
      PUBLIC errmes

    CONTAINS

      INTEGER FUNCTION errmes(itest,ierror,srname)
!-----------------------------------------------------------------------
! PURPOSE
!      ERRMES returns the value of IERROR or terminates the program,
!      printing a failure message

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  21 Jan 2003 (CG)

! ARGUMENTS in
!      ITEST   contains either 0 (hard fail) or 1 (soft fail).
!              any other entry gives hard fail.
!      IERROR  contains the number of the detected error
!      SRNAME  contains up to 8 characters - usually a library
!              routine name

! ARGUMENTS out
!      ERRMES  routine name, contains the value of IERROR

! ROUTINES called
!      can call auxiliary routine in some versions of library
!      ERUNIT and ADUNIT

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER :: ierror, itest
        CHARACTER (6) :: srname

        INTENT (IN) ierror, itest, srname

! Local variables

        INTEGER :: get, jtest, unit


        DATA get/0/, jtest/1/

!     Hard failure

        IF (itest==-99) THEN

!     To return Release message

          CALL adunit(unit,get,jtest)
          WRITE (unit,fmt=99002)

        ELSE IF (itest==1 .OR. ierror==0) THEN

!     Soft failure

          errmes = ierror
          IF (itest/=0 .AND. ierror/=0) THEN
            CALL adunit(unit,get,jtest)
            WRITE (unit,fmt=99001) srname, ierror
          END IF

        ELSE

          CALL erunit(unit,get,jtest)
          WRITE (unit,fmt=99001) srname, ierror
          STOP

        END IF
99001   FORMAT ('*** Error detected by Level 0 Library routine ',A,' - ITEST = ', &
          I5,//)
99002   FORMAT (' FELIB90 Release 1.0  -  1 JAN 2003')

      END FUNCTION errmes

    END MODULE mod_errmes
