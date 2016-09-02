    MODULE mod_vtol

      USE felib_globals, ONLY : wp
      USE mod_unflo,     ONLY : unflo
      USE mod_veps,      ONLY : veps

      PRIVATE
      PUBLIC vtol

    CONTAINS

      REAL (wp) FUNCTION vtol()
!-----------------------------------------------------------------------
! PURPOSE
!      VTOL returns the ratio of the smallest number which does not
!      cause underfow to the number VEPS which is defined by
!      1.D0+VEPS > 1.D0

!      ***************************************************
!      **********USES MACHINE DEPENDENT ROUTINES**********
!      **********   CODE IS MACHINE DEPENDENT   **********
!      ***************************************************

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  18 Jan 2003 (CG)

! ROUTINES called
!      UNFLO   VEPS

!************************************************************************


        IMPLICIT NONE

        vtol = unflo()/veps()

      END FUNCTION vtol

    END MODULE mod_vtol
