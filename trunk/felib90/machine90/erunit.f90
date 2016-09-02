
    MODULE mod_erunit

      USE felib_globals, only: erout
      USE felib_utils, only: setopt

      PRIVATE
      PUBLIC erunit

    CONTAINS

      SUBROUTINE erunit(nunit,ikey,itest)
!-----------------------------------------------------------------------
! PURPOSE
!      ERUNIT returns the value of the current error message UNIT
!      number or sets the current error message UNIT number
!      to a new value

!      *********************************************
!      **********MACHINE DEPENDENT ROUTINE**********
!      *********************************************

! HISTORY

!     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                          Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  18 Jan 2003 (CG)

! ARGUMENTS in
!      NUNIT   new UNIT number (IKEY = 0)
!      IKEY    action key (= 0 set NUNIT to UNIT - enquire,
!                          = 1 set UNIT to NUNIT - reset)
!      ITEST   standard error check flag

! ARGUMENTS out
!      NUNIT   old UNIT number
!      ITEST   error conditions returned or printed
!              - 1: invalid IKEY.
!              - 2: NUNIT not in valid range.

!***********************************************************************

!     For unix based machines

        IMPLICIT NONE

! Dummy arguments

        INTEGER :: ikey, nunit, zero
        INTEGER, OPTIONAL :: itest
        INTENT (IN) ikey
        INTENT (INOUT) itest, nunit

! Local variables

        INTEGER :: jtest, ktest, range1, range2
        CHARACTER (6) :: srname = 'ERUNIT'

        DATA range1/0/, range2/100/, zero/0/

        call setopt(jtest,zero,itest)

        IF (ikey==1) THEN

!   Reset channel number  

          IF ((nunit>=range1) .AND. (nunit<=range2)) THEN
            erout = nunit
          ELSE
            ktest = 2
          END IF

        ELSE IF (ikey==0) THEN

!   Enquiry

          nunit = erout

        ELSE

!   Invalid ikey

          ktest = 1

        END IF

!     Exit sequence if error is detected

        IF (jtest==1) THEN
          jtest = ktest
        ELSE IF (jtest==0 .AND. ktest/=0) THEN
          WRITE (erout,fmt=99001) srname, jtest
          STOP
        END IF

99001   FORMAT (' *** Error detected by Level 0 Library routine ',A8, &
          ' - ITEST = ',I5,//)

      END SUBROUTINE erunit

    END MODULE mod_erunit
