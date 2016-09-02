    MODULE felib_utils

! MODULE felib_utils is the main defining modules of FELIB90. All
! user callable routines are included here.

! System

      USE felib_globals, ONLY : wp, nin=>stdin, nout=>stdout

! Routines

      USE mod_space,   ONLY : create, destroy
      USE mod_setopt,   ONLY : setopt

    END MODULE felib_utils 

