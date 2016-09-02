C
      SUBROUTINE BQTRI(ABSS,IABSS,JABSS,WORK,IWORK,NQP,SIDNUM,COEF,
     *                 ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      BQTRI forms an equivalent one-dimensional quadrature rule of
C      a given two-dimensional rule for integration along the
C      specified side of a triangular element.
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
C     SUBROUTINE BQTRI(ABSS,IABSS,JABSS,WORK,IWORK,NQP,SIDNUM,COEF,
C    *                 ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IABSS,IERROR,ITEST,IWORK,JABSS,NQP,SIDNUM
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,COEF,ETA,L1,L2,L3,LL1,LL2,LL3,WORK,XI
      DIMENSION ABSS(IABSS,JABSS),WORK(IWORK)
C
      EXTERNAL ERRMES
C
      DATA SRNAME/'BQTRI'/
C
C     Statement functions
C
      XI(LL1,LL2,LL3) = 1.0D0 - 3.0D0/2.0D0* (LL2+LL3)
      ETA(LL1,LL2,LL3) = DSQRT(3.0D0)/2.0D0* (LL3-LL2)
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IABSS.LE.2 .OR. JABSS.LE.NQP) IERROR = 4
         IF (IWORK.LE.NQP) IERROR = 3
         IF (SIDNUM.LE.0 .OR. SIDNUM.GT.3) IERROR = 2
         IF (NQP.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      COEF = 0.5D0
      GO TO (1000,1020,1040) SIDNUM
C
C     Side number 1 : L3 = 0
C
 1000 CONTINUE
      L3 = 0.0D0
      DO 1010 I = 1,NQP
         L1 = 0.5D0* (1.D0+WORK(I))
         L2 = 1.0D0 - L1
         ABSS(1,I) = XI(L1,L2,L3)
         ABSS(2,I) = ETA(L1,L2,L3)
 1010 CONTINUE
      RETURN
C
C     Side number 2 : L1 = 0
C
 1020 CONTINUE
      L1 = 0.0D0
      DO 1030 I = 1,NQP
         L2 = 0.5D0* (1.D0+WORK(I))
         L3 = 1.0D0 - L2
         ABSS(1,I) = XI(L1,L2,L3)
         ABSS(2,I) = ETA(L1,L2,L3)
 1030 CONTINUE
      RETURN
C
C     Side number 3 : L2 = 0
C
 1040 CONTINUE
      L2 = 0.0D0
      DO 1050 I = 1,NQP
         L3 = 0.5D0* (1.D0+WORK(I))
         L1 = 1.0D0 - L3
         ABSS(1,I) = XI(L1,L2,L3)
         ABSS(2,I) = ETA(L1,L2,L3)
 1050 CONTINUE
      END
