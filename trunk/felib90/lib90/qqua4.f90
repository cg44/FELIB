!
MODULE mod_qqua4

USE FELIB_GLOBALS,only : wp
USE mod_SPACE,only : create
USE mod_setopt, only : setopt

PRIVATE
PUBLIC qqua4

contains

SUBROUTINE qqua4(wght,abss,nqp,itest)

IMPLICIT NONE

! Dummy arguments

INTEGER :: nqp
INTEGER, OPTIONAL :: itest
REAL (wp), POINTER  :: wght(:), abss(:,:)

INTENT (INOUT) :: itest
INTENT (OUT) :: nqp

INTEGER :: dimen, ktest
REAL (wp) :: area
CHARACTER(6) :: srname

DATA srname/'QQUA4'/
INTRINSIC sqrt
!
nqp = 4
dimen =2 
!
call setopt(ktest,0,itest)

CALL create(wght,nqp)
CALL create(abss,dimen,nqp)

!
!     Main body
!
area = 4.0D0
wght(1) = area*0.25D+0
wght(2) = wght(1)
wght(3) = wght(1)
wght(4) = wght(1)
!
abss(1,1) = -1.0D0/sqrt(3.0D0)
abss(1,2) = abss(1,1)
abss(1,3) = +1.0D0/sqrt(3.0D0)
abss(1,4) = abss(1,3)
abss(2,1) = -1.0D0/sqrt(3.0D0)
abss(2,2) = +1.0D0/sqrt(3.0D0)
abss(2,3) = abss(2,2)
abss(2,4) = abss(2,1)
!
END SUBROUTINE QQUA4
END MODULE mod_qqua4
