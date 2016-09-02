    MODULE felib90

! MODULE felib90 is the main defining modules of FELIB90. All
! user callable routines are included here.

! System

      USE felib_globals, ONLY : wp, nin=>stdin, nout=>stdout
      USE mod_space

! Routines

      USE mod_bndwth,   ONLY : bndwth
      USE mod_qqua4,    ONLY : qqua4
      USE mod_elgeom,   ONLY : elgeom
      USE mod_quam4,    ONLY : quam4
      USE mod_quam8,    ONLY : quam8
      USE mod_scaprd,   ONLY : scaprd
      USE mod_matmul,   ONLY : matmul
      USE mod_matnul,   ONLY : matnul
      USE mod_matran,   ONLY : matran
      USE mod_matvec,   ONLY : matvec
      USE mod_prtval,   ONLY : prtval
      USE mod_asrhs,    ONLY : asrhs
      USE mod_assym,    ONLY : assym
!     USE mod_gausol,   ONLY : gausol
      USE mod_chosol,   ONLY : chosol
      USE mod_direct,   ONLY : direct
      USE mod_matinv,   ONLY : matinv
      USE mod_getgeo,   ONLY : getgeo
      USE mod_gettop,   ONLY : gettop
      USE mod_prtgeo,   ONLY : prtgeo
      USE mod_prttop,   ONLY : prttop
      USE mod_errmes,   ONLY : errmes
      USE mod_asful,    ONLY : asful
      USE mod_vecnul,   ONLY : vecnul
      USE mod_vecadd,   ONLY : vecadd
      USE mod_matadd,   ONLY : matadd
      USE mod_setopt,   ONLY : setopt

      USE mod_matmul_intrinsic, ONLY : matrix_multiply => matmul

    END MODULE felib90

