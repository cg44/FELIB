!*****************************************************************

    PROGRAM seg3p1

!*****************************************************************

!     Copyright (C) 2003 : CLRC, Rutherford Appleton Laboratory
!     Chilton, Didcot, Oxfordshire OX11 0QX

! N.B. The working precision of the current library is held
!      in the variable wp. This must be used in all REAL
!      declarations of variables used by FELIB90.

!      The program also uses the standard FELIB90 values for 
!      nin and nout.

!*****************************************************************

      USE felib90 ! Use FELIB90 all routines

      USE def3p1 ! Use standard SEG3P1 definitions

      IMPLICIT NONE

! Parameters

      REAL (wp), PARAMETER :: scale = 1.0E+10

! Local variables

      INTEGER :: bndnod, dimen, dofel, dofnod, hband, i, iquad, itest, j, &
        nele, nodel, nqp, totdof, totels, totnod
      REAL (wp) :: det, eta, quot, strgth, x, xi, y

! Allocatable arrays - mesh size dependent

      INTEGER, POINTER :: bnode(:), nf(:,:), eltop(:,:)
      REAL (wp), POINTER :: bval(:), rhs(:), coord(:,:), sysk(:,:)

! Intrinsic functions

      INTRINSIC abs

!     Initialisation of POINTERS to main arrays

      NULLIFY (bnode,nf,eltop)
      NULLIFY (bval,rhs,coord,sysk)

      itest = 0

!     **********************
!     *                    *
!     * Input Data Section *
!     *                    *
!     **********************

!     Input of nodal geometry

      CALL getgeo(coord,totnod,dimen)
      CALL prtgeo(coord)

!     Input of element topology

      CALL gettop(eltop,totels)
      CALL prttop(eltop)

!     Input of permeabilities, construction of permeability matrix P
!     and source strength

      CALL matnul(p,dimen,dimen)
      WRITE (nout,'(/A)') 'Permeabilities'
      READ (nin,'(2F10.0)') (p(i,i),i=1,dimen)
      WRITE (nout,'(2F10.5)') (p(i,i),i=1,dimen)

      WRITE (nout,'(/A)') 'Source Strength'
      READ (nin,'(F10.0)') strgth
      WRITE (nout,'(F10.5)') strgth

!     Input of number of degrees of freedom per node, input of
!     boundary conditions and construction of nodal freedom array NF

      WRITE (nout,'(/A)') 'Degrees of freedom per node (DOFNOD)'
      READ (nin,'(I5)') dofnod
      WRITE (nout,'(I5)') dofnod

!    Input boundary conidtions

      WRITE (nout,'(/A)') 'Boundary Conditions'
      READ (nin,'(I5)') bndnod
      WRITE (nout,'(I5)') bndnod

      CALL vecnul(bnode,bndnod)
      CALL vecnul(bval,bndnod)
      DO i = 1, bndnod
        READ (nin,'(I5,F10.0)') bnode(i), bval(i)
        WRITE (nout,'(I5,F10.5)') bnode(i), bval(i)
      END DO

! Setup nodel freedom array


      CALL matnul(nf,totnod,dofnod)
      totdof = 0
      DO i = 1, totnod
        DO j = 1, dofnod
          totdof = totdof + 1
          nf(i,j) = totdof
        END DO
      END DO

!     Calculation of semi-bandwidth

      CALL bndwth(eltop,nf,hband)

!  ************************************
!  *                                  *
!  * System Stiffness Matrix Assembly *
!  *                                  *
!  ************************************


! System matrices setup and initalise : rhs, sysk

      CALL matnul(sysk,totdof,hband)
      CALL vecnul(rhs,totdof) 

! Setup quadrature

      CALL qqua4(wght,abss,nqp)

! Begin main element loop

      DO nele = 1, totels !Loop over all elements

        nodel = eltop(nele,2)
        dofel = dofnod*nodel
        CALL elgeom(nele,eltop,coord,geom)

!     Integration loop for element stiffness using NQP quadrature
!     points

        CALL matnul(elk,dofel,dofel) 
        CALL vecnul(elq,dofel)
        CALL vecnul(scvec,dofel)

        DO iquad = 1, nqp !Numerical integration

!     Form linear shape function and space derivatives in the local
!     corrdinates. Transform local derivatives to global coordinate
!     system

          xi = abss(1,iquad)
          eta = abss(2,iquad)
          CALL quam4(fun,lder,xi,eta)

          CALL matran(geom,geomt)
          CALL matvec(geomt,fun,xy) 
          x = xy(1)
          y = xy(2)

          CALL matmul(lder,geom,jac) 
          CALL matinv(jac,jacin,det)
          CALL matmul(jacin,lder,gder)

!     Formation of element stiffness ELK

          CALL matmul(p,gder,pd)
          CALL matran(gder,gdert)
          CALL matmul(gdert,pd,dtpd)

          quot = abs(det)*wght(iquad)
          dtpd = dtpd*quot
          scvec = fun*src(x,y,strgth)*quot

          CALL matadd(elk,dtpd)
          CALL vecadd(elq,scvec)

        END DO !Loop over quadrature points - iquad

!     Assembly of system stiffness matrix

        CALL direct(nele,eltop,nf,steer)
        CALL assym(sysk,elk,steer)
!       CALL asusm(sysk,elk,steer)
        CALL asrhs(rhs,elq,steer)

      END DO !Loop over elements - nele

!     *********************
!     *                   *
!     * Equation Solution *
!     *                   *
!     *********************

!     Modification of stiffness matrix and right-hand side to
!     implement boundary conditions

      DO i = 1, bndnod
        j = bnode(i)
        sysk(j,hband) = sysk(j,hband)*scale
        rhs(j) = sysk(j,hband)*bval(i)
      END DO

!     Solution of system matrix for the nodal values of the
!     potential

       CALL chosol(sysk,rhs)
!     CALL gausol(sysk,rhs)
      WRITE (nout,'(/A)') 'Nodal Potentials'
      CALL prtval(rhs,nf)

      STOP

    CONTAINS

!*****************************************************************
! Source function

      FUNCTION src(x,y,strgth)

!       USE felib90

        IMPLICIT NONE

! Dummy arguments

        REAL (wp) :: src
        REAL (wp) :: strgth, x, y
        INTENT (IN) strgth, x, y

        src = 0.0D0
        IF ((x>1.0D0) .AND. (x<2.0D0) .AND. (y>1.0D0) .AND. (y<2.0D0)) &
          src = strgth

      END FUNCTION src

    END PROGRAM seg3p1
