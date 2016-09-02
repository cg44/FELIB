!
MODULE mod_quam4

USE FELIB_GLOBALS,only : wp
use mod_space,only : create

PRIVATE
PUBLIC quam4

contains

SUBROUTINE quam4(fun,der,xi,eta,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      QUAM4 returns the values of shape functions and their
!      derivatives at a specified point for a 4-noded c0
!      continuous quadrilateral element.  The approximated
!      function will be continuous across element boundaries
!
! HISTORY
!
!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.0   2 Jul 2000 (CG)
!
! ARGUMENTS in
!      XI      first local coordinate value
!      ETA     second local coordinate value
!
! ARGUMENTS out
!      FUN     vector of length IFUN. FUN(I) contains the
!              value of the i'th shape function at the point
!              (XI, ETA)
!      DER     DER(I, J) contains the value of the derivative of 
!              the j'th shape function with respect to the i'th
!              coordinate at the point (XI, ETA)
!
! ROUTINES called
!
!***********************************************************************

IMPLICIT NONE

! Dummy arguments

INTEGER, optional :: itest 
REAL (wp), POINTER :: der(:,:), fun(:)
REAL (wp) :: eta , xi

INTENT (IN) :: eta, xi
INTENT (INOUT) itest

! Local variables

INTEGER, PARAMETER :: dimen=2, nodel=4

INTEGER :: ierror
REAL (wp) :: dummy , etam , etap , val , xim , xip
CHARACTER(6) :: srname

INTRINSIC abs
!
DATA srname/'QUAM4'/
!
!     Parameter checking
!
CALL create(fun,nodel)
CALL create(der,dimen,nodel)

!     Main body
!
etam = 0.250D0*(1.0D0-eta)
etap = 0.250D0*(1.0D0+eta)
xim = 0.250D0*(1.0D0-xi)
xip = 0.250D0*(1.0D0+xi)
!
fun(1) = 4.0D0*xim*etam
fun(2) = 4.0D0*xim*etap
fun(3) = 4.0D0*xip*etap
fun(4) = 4.0D0*xip*etam
!
der(1,1) = -etam
der(2,1) = -xim
der(1,2) = -etap
der(2,2) = xim
der(1,3) = etap
der(2,3) = xip
der(1,4) = etam
der(2,4) = -xip
!
END SUBROUTINE quam4
END MODULE mod_quam4
