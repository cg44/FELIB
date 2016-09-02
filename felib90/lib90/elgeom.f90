
    MODULE mod_elgeom

      USE felib_globals, only : wp
      USE mod_errmes, only : errmes
      USE mod_space,only : create

      PRIVATE
      PUBLIC elgeom

    CONTAINS

      SUBROUTINE elgeom(nele,eltop,coord,geom,dimen,itest)

!-----------------------------------------------------------------------!
! PURPOSE
!      ELGEOM constructs the element geometry array for the specified
!      element in the local node order

! HISTORY

!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  May 2002 (CG)

! ARGUMENTS in
!      NELE    element number for which the nodal coordinate
!              array is required
!      ELTOP   contains element type, number of nodes, and
!              element topologies for all the elements
!      COORD   COORD(I, J) contains the j'th global coordinate
!              of the i'th node
!      DIMEN   dimensionality of the problem (Optional). Taken as
!              size of COORD.
!      ITEST   error checking flag (Optional). Taken as ITEST=0

! ARGUMENTS out
!      GEOM    contains global coordinates of the nodes
!              associated with element NELE in the local node
!              numbering order

! ROUTINES called
!      ERRMES
!***********************************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER, OPTIONAL :: dimen, itest
        INTEGER :: nele
        INTEGER, POINTER :: eltop(:,:)
        REAL (wp), POINTER :: coord(:,:), geom(:,:)

! Local variables

        INTEGER :: ierror, inod, jnod, jtest, ktest, nodel,  &
          ieltop, jeltop, icoord, jcoord, ldimen, igeom, jgeom

        CHARACTER (6) :: srname = 'ELGEOM'

! Checking optional arguments

        IF (present(itest)) THEN
          ktest = itest
        ELSE
          ktest = 0
        END IF
        IF (present(dimen)) THEN
          ldimen = dimen
        ELSE
          ldimen = size(coord,2)
        END IF

! Collecting parammeter data

        icoord = size(coord,1)
        jcoord = size(coord,2)
        ieltop = size(eltop,1)
        jeltop = size(eltop,2)

!     Parameter checking

        jtest = ktest
        IF (jtest/=-1) THEN
          ierror = 0
          IF (jcoord<ldimen) ierror = 3
          IF (ieltop<nele) ierror = 2
          IF (nele<=0 .OR. ldimen<=0) ierror = 1
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) THEN
            IF (present(itest)) THEN
              itest = ktest
              RETURN
            END IF
          END IF
        END IF

        nodel = eltop(nele,2)

! Create space

        CALL create(geom,nodel,ldimen)

        igeom = size(geom,1)
        jgeom = size(geom,2)

!     Range checking on NODEL : NODEL+2 <= JELTOP and NODEL <= IGEOM

        IF (jtest/=-1) THEN
          ierror = 0
          IF (jgeom<ldimen) ierror = 4
          IF (jeltop<nodel+2) ierror = 5
          IF (igeom<nodel) ierror = 7
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) THEN
            IF (present(itest)) THEN
              itest = ktest
              RETURN
            END IF
          END IF
        END IF

        DO jnod = 1, nodel
          inod = eltop(nele,jnod+2)

!     Range checking on INOD : INOD <= ICOORD

          IF (jtest/=-1) THEN
            ierror = 0
            IF (icoord<inod) ierror = 6
            ktest = errmes(jtest,ierror,srname)
            IF (ktest/=0) THEN
              IF (present(itest)) THEN
                itest = ktest
                RETURN
              END IF
            END IF
          END IF

!     Construct geometry array GEOM

            geom(jnod,1:ldimen) = coord(inod,1:ldimen)
        
        END DO

      END SUBROUTINE elgeom
    END MODULE mod_elgeom
