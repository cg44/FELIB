C
      SUBROUTINE BQQUA(ABSS,IABSS,JABSS,WORK,IWORK,NQP,SIDNUM,COEF,
     *                 ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      BQQUA forms an equivalent one-dimensional quadrature rule of
C      a given two-dimensional rule for integration along the
C      specified side of a rectangular element.
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
C      ABSS    array holding the abscissae of the two-
C              dimensional quadrature rule to be used
C      IABSS   first dimension of array ABSS
C      JABSS   second dimension of array ABSS
C      IWORK   dimension of array WORK
C      NQP     number of quadrature points in the rule
C      SIDNUM  the side number for which the one-
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
C     SUBROUTINE BQQUA(ABSS,IABSS,JABSS,WORK,IWORK,NQP,SIDNUM,COEF,
C    *                 ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IABSS,IERROR,ITEST,IWORK,J,JABSS,K,NQP,SIDNUM
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,COEF,VAL,WORK
      DIMENSION ABSS(IABSS,JABSS),WORK(IWORK)
C
      EXTERNAL ERRMES
C
      DATA SRNAME/'BQQUA'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IABSS.LT.2 .OR. JABSS.LT.NQP) IERROR = 4
         IF (IWORK.LT.NQP) IERROR = 3
         IF (SIDNUM.LE.0 .OR. SIDNUM.GT.4) IERROR = 2
         IF (NQP.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      COEF = 1.0D0
      GO TO (1000,1010,1020,1030) SIDNUM
C
C     Side number 1 : xi = -1
C
 1000 CONTINUE
      I = 1
      J = 2
      VAL = -1.0D0
      GO TO 1040
C
C     Side number 2 : xi = +1
C
 1010 CONTINUE
      I = 2
      J = 1
      VAL = 1.0D0
      GO TO 1040
C
C     Side number 3 : eta = -1
C
 1020 CONTINUE
      I = 1
      J = 2
      VAL = 1.0D0
      GO TO 1040
C
C     Side number 4 : eta = +1
C
 1030 CONTINUE
      I = 2
      J = 1
      VAL = -1.0D0
C
 1040 CONTINUE
      DO 1050 K = 1,NQP
         ABSS(I,K) = VAL
         ABSS(J,K) = WORK(K)
 1050 CONTINUE
C
      END
