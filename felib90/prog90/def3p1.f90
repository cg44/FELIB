    MODULE def3p1

! All arrays within FELIB programs are defined as POINTERs. This
! enables dynamic allocation within FUNCTIONS and SUBROUTINES. Many
! routines automatically ALLOCATE space of a current set of intermediate
! variables for the solution process.

! These variables are created through a library of space management routines
! which includes garbage collection. These routines do not inhibit the user
! defining his own arrays locally or by using the space creataion routines.

! This modules contains a standard set of variable definitions
! often found in basic FELIB programs.

      USE felib_globals, ONLY : wp


! Allocatable arrays - dependent on element types and problem 
!                      dimensionality

      INTEGER, POINTER ::     &
        steer(:)   => null()    ! Element steering vector
      REAL (wp), POINTER ::   & 
        elq(:)     => null(), & ! Element vector                              
        fun(:)     => null(), & ! Shape function vector
        xy(:)      => null(), & ! Global coordinate vector
        geom(:,:)  => null(), & ! Local geometry array
        elk(:,:)   => null(), & ! Element stiffness arrat
        lder(:,:)  => null(), & ! Local derivatives of shape funtions
        jac(:,:)   => null(), & ! Transformation jacobian
        jacin(:,:) => null(), & ! Inverse of JAC
        geomt(:,:) => null(), & ! Element geometry transposed
        wght(:)    => null(), & ! Quadrature weights
        abss(:,:)  => null(), & ! Quadrature abssise
        p(:,:)     => null(), & ! Permeitvity array P
        pd(:,:)    => null(), & ! P transposed
        scvec(:)   => null(), & ! Element source vector
        gder(:,:)  => null(), & ! Global derivatives of shape functions
        dtpd(:,:)  => null(), & ! Element matrix
        gdert(:,:) => null()    ! Tranpose of GDERT


    END MODULE def3p1
