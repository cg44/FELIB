C***********************************************************************

      DOUBLE PRECISION FUNCTION VTOL(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE RATIO OF THE SMALLEST NUMBER WHICH DOES NOT
C      CAUSE UNDERFOW TO THE NUMBER VEPS WHICH IS DEFINED BY
C      1.DO+VEPS > 1.D0
C
C      ***************************************************
C      **********USES MACHINE DEPENDENT ROUTINES**********
C      **********  CODE IS MACHINE INDEPENDENT  **********
C      ***************************************************
C
C HISTORY
C
C    COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                         CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    23 OCT 1980 (KR)
C
C ARGUMENTS IN
C      X       DUMMY ARGUMENT
C
C ROUTINES CALLED
C      UNFLO   VEPS
C
C
C     DOUBLE PRECISION FUNCTION VTOL(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE
      DOUBLE PRECISION UNFLO, VEPS
C
      VTOL = UNFLO(VALUE)/VEPS(VALUE)
      RETURN
      END
