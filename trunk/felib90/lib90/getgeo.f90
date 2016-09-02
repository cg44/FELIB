
    MODULE mod_getgeo

      USE felib_globals, ONLY : wp , stdin
      USE mod_space
      USE mod_setopt

      PRIVATE
      PUBLIC getgeo

    CONTAINS

      SUBROUTINE getgeo(coord,totnod,dimen,nin,f1,f2,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      GETGEO read the element geometry in a standard format 
!      from the channel NIN.

!      The File must be structured:
!         TOTNOD, DIMEN      -  FORMAT(2I3)
!         NODNUM X Y [Z]     -  FORMAT(I3,3F10.5)
!      COORD will be allocated if necessary

! HISTORY

!      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0 May 2003 (CG)

! ARGUMENTS in
!      TOTNOD  total number of nodes in the mesh
!      DIMEN   dimensionality of the geometric data
!      COORD   array of dimension (ICOORD, JCOORD) containing
!              global coordinates of the nodes
!      NIN    fortran unit number (optional)
!      ITEST   error checking option (optional)

! ROUTINES called
!      ERRMES

!***********************************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER :: dimen, totnod
        INTEGER, OPTIONAL :: nin, itest
        REAL (wp), POINTER :: coord(:,:)
        CHARACTER (*), OPTIONAL :: f1,f2

        INTENT (IN) nin
        INTENT (OUT) dimen, totnod
        INTENT (INOUT) itest

! Local variables

        INTEGER :: i, j, k, chan, ktest
        CHARACTER (80) :: ff1, ff2
        CHARACTER (*), PARAMETER :: FORMAT1='(2I5)'
        CHARACTER (*), PARAMETER :: FORMAT2='(I5,3F10.5)'
        CHARACTER (*), PARAMETER :: srname = 'GETGEO'

!     Main body

        CALL SETOPT(chan,stdin,nin)
        CALL SETOPT(ktest,0,itest)
        CALL SETOPT(ff1,FORMAT1,f1)
        CALL SETOPT(ff2,FORMAT2,f2)
        
! Read size of coord data

        READ (chan,ff1) totnod, dimen

!     Space allocation

        CALL create(coord,totnod,dimen)

        DO i = 1, totnod
          READ (chan,fmt=ff2) k, (coord(k,j),j=1,dimen)
        END DO

      END SUBROUTINE getgeo
    END MODULE mod_getgeo
