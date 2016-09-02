
    MODULE mod_bndwth

      USE felib_globals
      USE mod_errmes
      USE mod_maxint

      PRIVATE
      PUBLIC bndwth

    CONTAINS

      SUBROUTINE bndwth(eltop,nf,hband,dofnod,totels,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      BNDWTH calculates the maximum freedom number difference over
!      all the elements in the mesh.

! HISTORY

!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0  May 2002 (CG)

! ARGUMENTS in
!     ELTOP        ELTOP(I,1) = element type of element I
!                  ELTOP(I,2) = number of nodes on element I
!                  ELTOP(I,J+2), J=1(1)NUMBER of nodes on element,
!                  contains the nodes associated with element I
!      NF        NF(I,J) contains the freedom numbers associated
!                  with node I
!      OPT_DOFNOD  number of degrees of freedom per node on the
!                  element
!      OPT_TOTELS  number of elements in the mesh
!      OPT_ITEST   error checking option

! ARGUMENTS out
!      HBAND   semi bandwidth

!*********************************************************

        IMPLICIT NONE

! Dummy arguments

        INTEGER :: hband
        INTEGER, OPTIONAL :: dofnod, totels, itest
        INTEGER  :: eltop(:,:), nf(:,:)

        INTENT (IN) eltop, nf, dofnod, totels
        INTENT (OUT) hband
        INTENT (INOUT) itest

! Local variables

        INTEGER :: i, ideg, iele, ierror, inod, j, jtest, maxd, mind, nodel, &
          totnod, ltotels, ldofnod, ktest, jeltop
        CHARACTER (6) :: srname

        DATA srname/'BNDWTH'/

!     Initialisation

        hband = 0

! Checking optional arugments

        IF (present(itest)) THEN
          ktest = itest
        ELSE
          ktest = 0
        END IF
        IF (present(dofnod)) THEN
          ldofnod = dofnod
        ELSE
          ldofnod = size(nf,2)
        END IF
        IF (present(totels)) THEN
          ltotels = totels
        ELSE
          ltotels = size(eltop,1)
        END IF

!     Main scanning loops

        totnod = size(nf,1)
        jeltop = size(eltop,2)

        DO iele = 1, ltotels
          maxd = 0
          mind = maxint()
          nodel = eltop(iele,2)

!     Range checking on NODEL
          jtest=ktest
          IF (jtest/=-1) THEN
            ierror = 0
            IF (jeltop<nodel+2) ierror = 4
            ktest = errmes(jtest,ierror,srname)
            IF (ktest/=0) RETURN
          END IF

!     Loop over nodes

          DO i = 1, nodel
            inod = eltop(iele,i+2)

!     Range checking on INOD
            jtest=ktest
            IF (jtest/=-1) THEN
              ierror = 0
              IF (totnod<inod) ierror = 5
              ktest = errmes(jtest,ierror,srname)
              IF (ktest/=0) RETURN
            END IF

!     Loop over freedoms

            DO j = 1, ldofnod
              ideg = nf(inod,j)
              IF (ideg/=0) THEN
                maxd = max(ideg,maxd)
                mind = min(ideg,mind)
              END IF
            END DO
          END DO

!     Maximum freedom number difference

          hband = max(hband,maxd-mind)
        END DO

!     Semi band width

        hband = hband + 1

      END SUBROUTINE bndwth
    END MODULE mod_bndwth
