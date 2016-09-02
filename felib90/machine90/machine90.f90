    MODULE machine90

!-----------------------------------------------------------------------
! PURPOSE
!      The MACHINE90 modules provides a framework for the library
!      of machine dependent routines. These routines mirror those
!      in the Fortran 77 library. However it should be noted that
!      some are not strictly necessary as their function is provide
!      through Fortran 95 intrinsic functions.

!      FELIB90 provides them as capatability routines.

!      Each module subprogramme is USEd by MACHINE90 to build a
!      complete machine constant library.

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0 January 2003 (CG)

!***********************************************************************


      USE mod_adunit, ONLY : adunit
      USE mod_erunit, ONLY : erunit
      USE mod_maxint, ONLY : maxint
      USE mod_ovflo,  ONLY : ovflo
      USE mod_unflo,  ONLY : unflo
      USE mod_veps,   ONLY : veps
      USE mod_vtol,   ONLY : vtol

    END MODULE machine90

