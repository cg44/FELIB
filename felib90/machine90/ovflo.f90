    MODULE mod_ovflo

      USE felib_globals, ONLY : wp

      PRIVATE
      PUBLIC ovflo

    CONTAINS

      REAL (wp) FUNCTION ovflo()
!-----------------------------------------------------------------------
! PURPOSE
!      OVFLO returns the value of the largest positive real floating-
!      point number representable on the computer

!      *********************************************
!      *********MACHINE DEPENDENT ROUTINE***********
!      *********************************************

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  18 Jan 2003 (CG)

!***********************************************************************

! for UNIX based machine

        IMPLICIT NONE

        INTRINSIC huge

        ovflo = huge(ovflo)

      END FUNCTION ovflo

    END MODULE mod_ovflo
