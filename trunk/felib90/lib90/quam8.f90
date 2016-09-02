!
MODULE mod_quam8

USE FELIB_GLOBALS,only : wp
use mod_space,only : create

PRIVATE
PUBLIC quam8

contains

SUBROUTINE quam8(fun,der,xi,eta,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      QUAM8 returns the values of shape functions and their
!      derivatives at a specified point for an 8-noded c0
!      continuous quadrilateral element. The approximated
!      function will be continuous across element boundaries
!
! HISTORY
!
!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.0  29 Oct 2000 (CG)
!
! ARGUMENTS in
!      XI      first local coordinate
!      ETA     second local coordinate
!      ITEST   error checking option
!
! ARGUMENTS out
!      FUN     vector of length IFUN. FUN(I) contains the
!              value of the i'th shape function at (XI, ETA)
!              for i=1(1)8
!      DER     array of dimension (IDER, JDER). DER(I, J)
!              contains the value of the derivative of the j'th
!              shape function with respect to the i'th
!              coordinate, for i=1(1)2 and j=1(1)8
!
! ROUTINES called
!
!***********************************************************************

IMPLICIT NONE

! Dummy arguments

INTEGER , optional :: itest
REAL (wp) :: eta , xi
REAL (wp), POINTER :: der(:,:), fun(:)

INTENT (IN) :: xi, eta
INTENT (INOUT) :: itest

! Local variables

INTEGER :: ierror
REAL (wp) :: dummy , etam , etap , val , xim , xip
CHARACTER(6) :: srname


INTRINSIC abs
!
DATA srname/'QUAM8'/
!
!     Parameter checking
!

!
!     Main body
!
etam = 0.250D0*(1.0D0-eta)
etap = 0.250D0*(1.0D0+eta)
xim = 0.250D0*(1.0D0-xi)
xip = 0.250D0*(1.0D0+xi)
!
der(1,1) = etam*(2.0D0*xi+eta)
der(1,2) = -8.0D0*etam*etap
der(1,3) = etap*(2.0D0*xi-eta)
der(1,4) = -4.0D0*etap*xi
der(1,5) = etap*(2.0D0*xi+eta)
der(1,6) = 8.0D0*etap*etam
der(1,7) = etam*(2.0D0*xi-eta)
der(1,8) = -4.0D0*etam*xi
der(2,1) = xim*(xi+2.0D0*eta)
der(2,2) = -4.0D0*xim*eta
der(2,3) = xim*(2.0D0*eta-xi)
der(2,4) = 8.0D0*xim*xip
der(2,5) = xip*(xi+2.0D0*eta)
der(2,6) = -4.0D0*xip*eta
der(2,7) = xip*(2.0D0*eta-xi)
der(2,8) = -8.0D0*xim*xip
!
fun(1) = 4.0D0*etam*xim*(-xi-eta-1.0D0)
fun(2) = 32.0D0*xim*etam*etap
fun(3) = 4.0D0*etap*xim*(-xi+eta-1.0D0)
fun(4) = 32.0D0*xim*xip*etap
fun(5) = 4.0D0*xip*etap*(xi+eta-1.0D0)
fun(6) = 32.0D0*xip*etap*etam
fun(7) = 4.0D0*xip*etam*(xi-eta-1.0D0)
fun(8) = 32.0D0*xim*xip*etam
!
END SUBROUTINE quam8
END MODULE mod_quam8
