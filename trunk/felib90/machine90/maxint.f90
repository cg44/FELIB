    MODULE mod_maxint

      PRIVATE
      PUBLIC maxint

    CONTAINS

      INTEGER FUNCTION maxint()
!-----------------------------------------------------------------------
! PURPOSE
!      MAXINT returns the value of the largest integer capable of
!      being stored by the computer

!      *********************************************
!      **********MACHINE DEPENDENT ROUTINE**********
!      *********************************************

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.1  29 Jan 2003 (CG)

!************************************************************************

        IMPLICIT NONE

        INTRINSIC huge

        maxint = huge(maxint)

      END FUNCTION maxint

    END MODULE mod_maxint

