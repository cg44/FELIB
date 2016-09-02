
    MODULE mod_adunit

      USE felib_globals, ONLY: adout, erout
      USE felib_utils, ONLY: setopt

      PRIVATE
      PUBLIC adunit

    CONTAINS

      SUBROUTINE adunit(nunit,ikey,itest)
!-----------------------------------------------------------------------
! PURPOSE
!      ADUNIT returns the value of the current advisory message UNIT
!      number or sets the current advisory message UNIT number
!      to a new value

!      *********************************************
!      **********MACHINE DEPENDENT ROUTINE**********
!      *********************************************

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0 January 2003 (CG)

! ARGUMENTS in
!      NUNIT   new UNIT number (IKEY = 0)
!      IKEY    action key (= 0 set NUNIT to UNIT - inquire,
!                          = 1 set UNIT to NUNIT - reset)
!      ITEST   standard error check flag

! ARGUMENTS out
!      NUNIT   old UNIT number
!      ITEST   error conditions returned or printed
!              - 1: invalid IKEY.
!              - 2: NUNIT not in valid range.

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER :: ikey, nunit
        INTEGER, OPTIONAL :: itest

        INTENT (IN) ikey
        INTENT (INOUT) itest, nunit

! Local variables

        INTEGER :: jtest, ktest, range1, range2, zero
        CHARACTER (6) :: srname = 'ADUNIT'

        DATA range1/5/, range2/10/, zero/0/

        call setopt(jtest,zero,itest)

        IF (ikey==1) THEN

!   Reset channel number  

          IF ((nunit>=range1) .AND. (nunit<=range2)) THEN
             adout = nunit
          ELSE
            ktest = 2
          END IF

        ELSE IF (ikey==0) THEN

!   Inquiry

          nunit = adout

        ELSE

!   Invalid ikey

          ktest = 1

        END IF

!     Exit sequence if error is detected

        IF (jtest==1) THEN
          jtest = ktest
        ELSE IF (jtest==0 .AND. ktest/=0) THEN
          WRITE (erout,fmt=99001) srname, ktest
          STOP
        END IF

99001   FORMAT (' *** Error detected by Level 0 Library routine ',A8, &
          ' - ITEST = ',I5,//)

      END SUBROUTINE adunit

    END MODULE mod_adunit
