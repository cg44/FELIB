    MODULE mod_veps

      USE felib_globals, ONLY : wp

      PRIVATE
      PUBLIC veps

    CONTAINS

      REAL (wp) FUNCTION veps()
!-----------------------------------------------------------------------
! PURPOSE
!      VEPS returns the value of the smallest number VEPS such that
!      1.D0+VEPS > 1.D0 on this computer

!      *********************************************
!      **********MACHINE DEPENDENT ROUTINE**********
!      *********************************************

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  18 Jan 2003 (CG)

!***********************************************************************


        IMPLICIT NONE

        INTRINSIC epsilon

        veps = epsilon(veps)

      END FUNCTION veps

    END MODULE mod_veps
