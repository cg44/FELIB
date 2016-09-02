
! Private values for FELIB

    MODULE felib_globals

      IMPLICIT NONE

      INTRINSIC kind

      INTEGER, PARAMETER :: dp = kind(1.0D0)
      INTEGER, PARAMETER :: sp = kind(1.0E0)

      INTEGER, PARAMETER :: wp = dp

      INTEGER, PARAMETER :: stdin = 5, stdout = 6

      INTEGER :: erout = 6, adout = 6


    END MODULE felib_globals
