!
MODULE mod_prttop

USE FELIB_GLOBALS,only : stdout

PRIVATE
PUBLIC prttop

contains

SUBROUTINE prttop(eltop,totels,nout,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      PRTTOP prints element topologies in a standard format
!
! HISTORY
!
!      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.1  29 Oct 1979 (CG)
!      Commented    14 Oct 1980 (KR)
!      Release 4.0   2 Oct 1996 (CG)
!
! ARGUMENTS in
!      TOTELS  total number of elements in the mesh
!      ELTOP   integer array of dimension (IELTOP, JELTOP)
!              containing element topologies, element type, and
!              number of nodes on the element
!      NOUT    fortran unit number
!      ITEST   error checking option
!
!***********************************************************************
IMPLICIT NONE

! Dummy arguments

INTEGER , OPTIONAL :: nout, totels, itest
INTEGER , POINTER :: eltop(:,:)

INTENT (IN) nout , totels
INTENT (INOUT) itest

! Local variables

INTEGER :: i , ierror , j , jtest , k , l, totels1, chan

CHARACTER(6) :: srname='PRTTOP'

!     Parameter checking

if(present(totels))then
totels1=totels
else
totels1=size(eltop,1)
endif
if(present(nout)) then
chan=nout
else
chan=stdout
endif

! Main Body

WRITE (chan,FMT=99001)
WRITE (chan,FMT=99002) totels1
WRITE (chan,FMT=99003)

DO i = 1 , totels1
   l = eltop(i,2)
   k = l + 2
   WRITE (chan,FMT=99004) i , (eltop(i,j),j=1,k)
END DO

99001 FORMAT (//'Element Topology'/)
99002 FORMAT (/'Numer of elements = ',I3)
99003 FORMAT (/' ',2X,'Elem',4X,'Eltyp',4X,'Nodel',4X,'Nodes',/' ')
99004 FORMAT (' ',2X,I4,4X,I3,7X,I3,4X,9(I4,2X),/' ',27X,9(I4,2X))
!
END SUBROUTINE prttop
END MODULE mod_prttop
