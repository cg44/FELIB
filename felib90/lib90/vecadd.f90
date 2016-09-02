!
MODULE mod_vecadd

      USE felib_globals,only : wp
      USE mod_space,only : create

PRIVATE
PUBLIC vecadd

CONTAINS

SUBROUTINE vecadd(v,w,n,itest)
!-----------------------------------------------------------------------
! PURPOSE
!      VECADD adds the vector V to the vector W, storing the result
!      in V
!
! HISTORY
!
!      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.1  29 Oct 1979 (IMS)
!
! ARGUMENTS in
!      V       vector of length IV. On entry, contains the
!              values of one of the vectors to be added
!              in elements 1 to N
!      W       vector of length IW. Contains the values of the
!              second vector to be added in elements 1 to N
!      N       number of elements of V and W to be added
!      ITEST   error checking option
!
! ARGUMENTS out
!      V       vector of length IV. On exit, V(I) contains
!              the sum of V(I) and W(I) for I=1(1)N
!
!***********************************************************************
!
IMPLICIT NONE

! Dummy arguments

INTEGER , optional :: itest , n
real(wp) , DIMENSION(:) :: v
real(wp) , DIMENSION(:) :: w

INTENT (IN)  n 
INTENT (INOUT) itest 

! Local variables

INTEGER :: i , ierror
CHARACTER(6) :: srname

!
DATA srname/'VECADD'/
!
!     Parameter checking
!
!     Main body
!
   v = v + w
!
END SUBROUTINE vecadd
END MODULE mod_vecadd
