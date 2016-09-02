C
      SUBROUTINE BQBRK(ABSS,IABSS,JABSS,WORK,IWORK,JWORK,NQP,FACNUM,
     *                 COEF,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      BQBRK forms an equivalent two-dimensional quadrature rule of
C      a given three-dimensional rule for integration over the
C      specified face of a brick element.
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
C      Release 2.0   1 Feb 1981 (CG)
C      Commented    24 Jun 1981 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      ABSS    array holding the abscissae of the three-
C              dimensional quadrature rule to be used
C      IABSS   first dimension of array ABSS
C      JABSS   second dimension of array ABSS
C      IWORK   first dimension of array WORK
C      JWORK   second dimension of array WORK
C      NQP     number of quadrature points in the rule
C      FACNUM  the face number for which the two-
C              dimensional rule is required
C      ITEST   error checking option
C
C ARGUMENTS out
C      WORK    array containing the abscissae of the
C              equivalent one-dimensional rule
C      COEF    a multiplier of the rule
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE BQBRK(ABSS,IABSS,JABSS,WORK,IWORK,JWORK,NQP,FACNUM,
C    *                 COEF,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,FACNUM,I,IABSS,IERROR,ITEST,IWORK,J,JABSS,JWORK,K,
     *        L,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,COEF,VAL,WORK
      DIMENSION ABSS(IABSS,JABSS),WORK(IWORK,JWORK)
C
      EXTERNAL ERRMES
C
      DATA SRNAME/'BQBRK'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (NQP.EQ.0) IERROR = 1
         IF (FACNUM.LE.0 .OR. FACNUM.GT.6) IERROR = 2
         IF (IWORK.LT.2 .OR. JWORK.LT.NQP) IERROR = 3
         IF (IABSS.LT.3 .OR. JABSS.LT.NQP) IERROR = 4
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      GO TO (1000,1010,1020,1030,1040,1050) FACNUM
      GO TO 1020
C
C     Face number 1 : zeta = -1
C
 1000 CONTINUE
      I = 3
      J = 1
      K = 2
      VAL = -1.0D0
      COEF = 1.0D0
      GO TO 1060
C
C     Face number 2 : zeta = +1
C
 1010 CONTINUE
      I = 3
      J = 1
      K = 2
      VAL = 1.0D0
      COEF = 1.0D0
      GO TO 1060
C
C     Face number 3 : xi = -1
C
 1020 CONTINUE
      I = 1
      J = 2
      K = 3
      VAL = -1.0D0
      COEF = 1.0D0
      GO TO 1060
C
C     Face number 4 : eta = +1
C
 1030 CONTINUE
      I = 2
      J = 1
      K = 3
      VAL = 1.0D0
      COEF = 1.0D0
      GO TO 1060
C
C     Face number 5 : xi = +1
C
 1040 CONTINUE
      I = 1
      J = 2
      K = 3
      VAL = 1.0D0
      COEF = 1.0D0
      GO TO 1060
C
C     Face number 6 : eta = -1
C
 1050 CONTINUE
      I = 2
      J = 1
      K = 3
      VAL = -1.0D0
      COEF = 1.0D0
C
 1060 CONTINUE
      DO 1070 L = 1,NQP
         ABSS(I,L) = VAL
         ABSS(J,L) = WORK(1,L)
         ABSS(K,L) = WORK(2,L)
 1070 CONTINUE
C
      END
