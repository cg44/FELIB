C
      DOUBLE PRECISION FUNCTION VTOL(X)
C-----------------------------------------------------------------------
C PURPOSE
C      VTOL returns the ratio of the smallest number which does not
C      cause underfow to the number VEPS which is defined by
C      1.D0+VEPS > 1.D0
C
C      ***************************************************
C      **********USES MACHINE DEPENDENT ROUTINES**********
C      **********   CODE IS MACHINE DEPENDENT   **********
C      ***************************************************
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Release 3.0  18 Jan 1984 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      X       dummy argument
C
C ROUTINES called
C      UNFLO   VEPS
C
C      DOUBLE PRECISION FUNCTION VTOL(X)
C************************************************************************
C
      DOUBLE PRECISION X,VALUE
C
      DOUBLE PRECISION UNFLO,VEPS
      EXTERNAL UNFLO,VEPS
C
      VTOL = UNFLO(VALUE)/VEPS(VALUE)
C
      RETURN
      END
