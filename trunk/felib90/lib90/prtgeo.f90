
MODULE mod_prtgeo

  USE FELIB_GLOBALS,only : wp, stdout

  PRIVATE
  PUBLIC prtgeo

contains

   SUBROUTINE prtgeo(coord,totnod,dimen,nout,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      PRTGEO prints out element geometry in a standard format
!
! HISTORY
!
!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.0  May 2002 (CG)
!
! ARGUMENTS in
!      TOTNOD  total number of nodes in the mesh
!      DIMEN   dimensionality of the geometric data
!      COORD   array of dimension (ICOORD, JCOORD) containing
!              global coordinates of the nodes
!      ICOORD  first dimension of COORD (.GE. TOTNOD)
!      JCOORD  second dimension of COORD (.GE. DIMEN)
!      NOUT    fortran unit number
!      ITEST   error checking option
!
! ROUTINES called
!      ERRMES
!
!***********************************************************************

IMPLICIT NONE

! Dummy arguments

INTEGER , optional :: dimen , nout , totnod, itest
REAL (KIND=wp) , POINTER :: coord(:,:)

INTENT (IN) dimen , nout , totnod
INTENT (INOUT) itest

! Local variables

INTEGER :: i , ierror , j, dimen1, totnod1, ktest, chan

CHARACTER(6) :: srname='PRTGEO'

!     Parameter checking
if(present(dimen))then
dimen1=dimen
else
dimen1=size(coord,2)
endif
if(present(totnod))then
totnod1=totnod
else
totnod1=size(coord,1)
endif
if(present(nout)) then
chan=nout
else
chan=stdout
endif
!
!     Main body
!
WRITE (chan,FMT=99001)
WRITE (chan,FMT=99002) totnod1
WRITE (chan,FMT=99003) dimen1
IF ( dimen1==1 ) WRITE (chan,FMT=99004)
IF ( dimen1==2 ) WRITE (chan,FMT=99005)
IF ( dimen1==3 ) WRITE (chan,FMT=99006)
IF ( dimen1==4 ) WRITE (chan,FMT=99007)
!
DO i = 1 , totnod1
   WRITE (chan,FMT=99008) i , (coord(i,j),j=1,dimen1)
END DO
!

99001 FORMAT (//'Nodal Geometry'/)
99002 FORMAT ('Number of nodes = ',I3)
99003 FORMAT ('Dimension      = ',I3)
99004 FORMAT (/' ',2X,'Node',9X,'x1',/' ')
99005 FORMAT (/' ',2X,'Node',9X,'x1',12X,'x2',/' ')
99006 FORMAT (/' ',2X,'Node',9X,'x1',12X,'x2',12X,'x3',/' ')
99007 FORMAT (/' ',2X,'Node',9X,'x1',12X,'x2',12X,'x3',12X,'x4',/' ')
99008 FORMAT (' ',2X,I3,5X,4(D12.4,2X))

END SUBROUTINE prtgeo
END MODULE mod_prtgeo
