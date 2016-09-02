!
MODULE mod_scaprd

USE FELIB_GLOBALS,only : wp

PRIVATE
PUBLIC scaprd

contains

SUBROUTINE scaprd(v,w,prodct,n,itest)

!-----------------------------------------------------------------------
!
! PURPOSE
!      SCAPR forms the scalar product of two vectors V and W, storing
!      the result in PRODCT
!
! HISTORY
!
!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.0   2 July 2000 (CG)
!
! ARGUMENTS in
!      V       vector of length IV
!      W       vector of length IW
!
! ARGUMENTS out
!      PRODCT  contains: sigma I = 1 to N of (V(I)*W(I))
!
! ROUTINES called
!
!************************************************************************

IMPLICIT NONE

! Dummy arguments
INTEGER, OPTIONAL :: n, itest
REAL (wp) :: prodct
REAL (wp) :: v(:), w(:)

INTENT (IN) :: n
INTENT (INOUT) :: itest

! Local variables

CHARACTER(6) :: srname

!
DATA srname/'SCAPRD'/
!
!     Parameter checking
!
!     Main body
!
   prodct = DOT_PRODUCT(v,w)
!
END SUBROUTINE scaprd
END MODULE mod_scaprd
