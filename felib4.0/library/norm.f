C
      DOUBLE PRECISION FUNCTION NORM(RHS,IRHS,TOTDOF,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      NORM computes an L2 norm of a vector for use
C      in terminating non-linear iterations
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 3.0  29 Jun 1986 (CJH)
C      Release 4.0   2 Oct 1996 (CG)
C
C ARGUMENTS in
C      RHS     vector containing values
C      IRHS    dimension of vector RHS (IRHS .GE. TOTDOF)
C      TOTDOF  number of entries in RHS
C      ITEST   error checking option
C
C ARGUMENTS out
C      NORM    the value of the NORM (function value)
C
C ROUTINES called
C      ERRMES
C
C     DOUBLE PRECISION FUNCTION NORM(RHS,IRHS,TOTDOF,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,IRHS,ITEST,TOTDOF
      CHARACTER*6 SRNAME
      DOUBLE PRECISION RHS
      DIMENSION RHS(IRHS)
      DATA SRNAME/'NORM'/
C
C     Parameter checking
C
      NORM = 0.0D0
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (TOTDOF.GT.IRHS) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Compute NORM
C
      DO 1000 I = 1,TOTDOF
         NORM = NORM + RHS(I)**2
 1000 CONTINUE
      NORM = DSQRT(NORM)
C
      END
