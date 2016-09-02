    MODULE mod_unflo

      USE felib_globals, ONLY : wp

      PRIVATE
      PUBLIC unflo

    CONTAINS

      REAL (wp) FUNCTION unflo()
!-----------------------------------------------------------------------
! PURPOSE
!      UNFLO returns the value of the smallest positive real floating-
!      point number exactly representable on the computer

!      *********************************************
!      **********MACHINE DEPENDENT ROUTINE**********
!      *********************************************

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  18 Jan 2003 (CG)

!***********************************************************************

        IMPLICIT NONE

        INTRINSIC tiny

        unflo = tiny(unflo)

      END FUNCTION unflo

    END MODULE mod_unflo
