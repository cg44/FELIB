
    MODULE mod_direct

      USE mod_errmes
      USE mod_space

      PRIVATE
      PUBLIC direct

    CONTAINS

      SUBROUTINE direct(nele,eltop,nf,steer,dofnod,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      DIRECT constructs the steering vector to direct the assembly of a
!      system matrix

! HISTORY

!      Copyright (C) 2002 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX

!      Release 1.0   1 June 2002 (CG)

! ARGUMENTS in
!      NELE    element number
!      ELTOP   2D array containing element type, number of
!              nodes on the element, and the element topologies
!      NF      contains freedom numbers associated with each
!              node
!      DOFNOD  number of degrees of freedom at each node

! ARGUMENTS out
!      STEER   vector containing freedom numbers associated
!              with element NELE, arranged in local node order

! ROUTINES called
!      ERRMES

!***********************************************************************
        IMPLICIT NONE

! Dummy arguments

        INTEGER :: nele
        INTEGER, OPTIONAL :: dofnod, itest
        INTEGER, POINTER  :: eltop(:,:), nf(:,:), steer(:)

        INTENT (IN) nele, dofnod
        INTENT (INOUT) itest

! Local variables

        INTEGER :: ideg, ierror, inod, jnod, jtest, k, nodel, ldofnod, ktest
        INTEGER :: inf, jnf, ieltop, jeltop

        CHARACTER (6) :: srname = 'DIRECT'

! Checking optional arguments

        IF (present(dofnod)) THEN
          ldofnod = dofnod
        ELSE
          ldofnod = size(nf,2)
        END IF
        IF (present(itest)) THEN
          ktest = itest
        ELSE
          ktest = 0
        END IF

! Collect array sizes

        jnf = size(nf,2)
        inf = size(nf,1)
        ieltop = size(eltop,1)
        jeltop = size(eltop,2)

!     Parameter checking

        jtest = ktest
        IF (jtest/=-1) THEN
          ierror = 0
          IF (jnf<ldofnod) ierror = 3
          IF (ieltop<nele) ierror = 2
          IF (nele<=0 .OR. ldofnod<=0) ierror = 1
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) then
             if(present(itest)) itest=ktest
             RETURN
          endif
        END IF

! Get element size

        nodel = eltop(nele,2)


!     Range check on NODEL : number of node in element
        
        jtest=ktest
        IF (jtest/=-1) THEN
          ierror = 0
          IF (jeltop<nodel+2) ierror = 4
          ktest = errmes(jtest,ierror,srname)
          IF (ktest/=0) RETURN
        END IF

! Allocate space for steer if necessary

        CALL create(steer,nodel*ldofnod)

!     Main loops

        k = 1
        DO inod = 1, nodel
          jnod = eltop(nele,inod+2)

!     Construct steering vector

          DO ideg = 1, ldofnod
            steer(k) = nf(jnod,ideg)
            k = k + 1
          END DO
        END DO

      END SUBROUTINE direct
    END MODULE mod_direct
