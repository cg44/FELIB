
    MODULE mod_gettop

      USE felib_globals, only : wp, stdin
      USE mod_errmes,only : errmes
      USE mod_space,only : create

      PRIVATE
      PUBLIC gettop

    CONTAINS

      SUBROUTINE gettop(eltop,totels,nin,itest)


!-----------------------------------------------------------------------
! PURPOSE
!      GETTOP reads element topologies in a standard format

! HISTORY

!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  May 2002 (CG)

! ARGUMENTS in
!      TOTELS  total number of elements in the mesh
!      ELTOP   integer array of dimension (IELTOP, JELTOP)
!              containing element topologies, element type, and
!              number of nodes on the element
!      NIN    fortran unit number (optional)
!      ITEST   error checking option (optional)

! ROUTINES called
!      ERRMES


!      SUBROUTINE PRTTOP(TOTELS,ELTOP,IELTOP,JELTOP,NOUT,ITEST)
!***********************************************************************
        IMPLICIT NONE

! Dummy arguments

        INTEGER :: totels
        INTEGER, OPTIONAL :: itest, nin
        INTEGER , POINTER :: eltop(:,:)

        INTENT (IN) nin
        INTENT (OUT) totels
        INTENT (INOUT) itest

! Local variables

        INTEGER :: i, j, k, l, nodel, elnum, eltyp, chan, ktest
        CHARACTER (6) :: srname = 'GETTOP'

        IF (present(nin)) THEN
          chan = nin
        ELSE
          chan = stdin
        END IF
        IF (present(itest)) THEN
           ktest=itest
        ELSE
           ktest=0
        END IF

!     Get size of topology data

        READ (chan,fmt=99001) eltyp, totels, nodel

! Allocate space eltop(1:totels,1:node+2)

        CALL create(eltop,totels,nodel+2)

        DO i = 1, totels
          READ (chan,fmt=99001) elnum, (eltop(elnum,j+2),j=1,nodel)

! Should really check elnum against totels

          eltop(elnum,1) = eltyp
          eltop(elnum,2) = nodel

        END DO

99001   FORMAT (10I5)

      END SUBROUTINE gettop
    END MODULE mod_gettop
