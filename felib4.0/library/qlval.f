C
      SUBROUTINE QLVAL(DIAG,IDIAG,SUB,ISUB,N,EPS,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QLVAL calculates all the eigenvalues of a real symmetric
C      tridiagonal matrix or of a full real symmetric matrix
C      that has been reduced to tridiagonal form using HOUSE
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
C      Release 1.1  29 Oct 1979 (IMS)
C      Commented    15 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      DIAG    vector of length IDIAG.  On entry, contains the
C              diagonal elements of the tridiagonal matrix
C      IDIAG   dimension of vector DIAG (.GE. N)
C      SUB     vector of dimension ISUB.  On entry, contains
C              SUB-DIAGONAL elements of tridiagonal matrix in
C              elements SUB(2) to SUB(N). SUB(1) May be
C              arbitrary. Contents destroyed during execution
C              of HOUSE
C      ISUB    dimension of SUB (.GE. N)
C      N       order of tridiagonal matrix
C      EPS     smallest positive number such that 1.+EPS>1.
C      ITEST   error checking option
C
C ARGUMENTS out
C      DIAG    vector of dimension IDIAG.  on exit, contains
C              eigenvalues in ascending order
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QLVAL(DIAG,IDIAG,SUB,ISUB,N,EPS,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,I1,IDIAG,IERROR,II,ISUB,ITEST,J,JTEST,L,M,M1,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION B,C,DIAG,EPS,F,G,H,P,R,S,SUB
      DIMENSION DIAG(IDIAG),SUB(ISUB)
      DATA SRNAME/'QLVAL'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IDIAG.LT.N .OR. ISUB.LT.N) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      IF (N.NE.1) THEN
         DO 1000 I = 2,N
            SUB(I-1) = SUB(I)
 1000    CONTINUE
      END IF
C
      SUB(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
C
      DO 1080 L = 1,N
         J = 0
         H = EPS* (DABS(DIAG(L))+DABS(SUB(L)))
C
C     Look for small SUB-DIAGONAL element
C
         IF (B.LT.H) B = H
         DO 1010 M = L,N
            IF (DABS(SUB(M)).LE.B) GO TO 1020
 1010    CONTINUE
C
 1020    CONTINUE
         IF (M.NE.L) THEN
 1030       CONTINUE
            IF (J.EQ.30) THEN
               GO TO 1090
            ELSE
C
C     Form shift
C
               J = J + 1
               G = DIAG(L)
               H = DIAG(L+1) - G
               IF (DABS(H).GE.DABS(SUB(L))) THEN
C
                  P = 2.0D0*SUB(L)/H
                  R = DSQRT(P*P+1.0D0)
                  DIAG(L) = SUB(L)*P/ (1.0D0+R)
               ELSE
                  P = H*0.5D0/SUB(L)
                  R = DSQRT(P*P+1.0D0)
                  H = P + R
                  IF (P.LT.0.0D0) H = P - R
                  DIAG(L) = SUB(L)/H
               END IF
               H = G - DIAG(L)
               I1 = L + 1
               IF (I1.LE.N) THEN
C
                  DO 1040 I = I1,N
                     DIAG(I) = DIAG(I) - H
 1040             CONTINUE
               END IF
C
C     Ql transformation
C
               F = F + H
               P = DIAG(M)
               C = 1.0D0
               S = 0.0D0
               M1 = M - 1
               DO 1050 II = L,M1
                  I = M1 - II + L
                  G = C*SUB(I)
                  H = C*P
                  IF (DABS(P).LT.DABS(SUB(I))) THEN
                     C = P/SUB(I)
                     R = DSQRT(C*C+1.0D0)
                     SUB(I+1) = S*SUB(I)*R
                     S = 1.0D0/R
                     C = C/R
                  ELSE
                     C = SUB(I)/P
                     R = DSQRT(C*C+1.0D0)
                     SUB(I+1) = S*P*R
                     S = C/R
                     C = 1.0D0/R
                  END IF
                  P = C*DIAG(I) - S*G
                  DIAG(I+1) = H + S* (C*G+S*DIAG(I))
 1050          CONTINUE
C
               SUB(L) = S*P
               DIAG(L) = C*P
               IF (DABS(SUB(L)).GT.B) GO TO 1030
            END IF
         END IF
C
C     Order eigenvalue
C
         P = DIAG(L) + F
         IF (L.NE.1) THEN
            DO 1060 II = 2,L
               I = L - II + 2
               IF (P.GE.DIAG(I-1)) THEN
                  GO TO 1070
               ELSE
                  DIAG(I) = DIAG(I-1)
               END IF
 1060       CONTINUE
         END IF
         I = 1
 1070    CONTINUE
         DIAG(I) = P
 1080 CONTINUE
C
      ITEST = 0
      RETURN
C
 1090 CONTINUE
      IF (JTEST.NE.-1) THEN
         IERROR = 3
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
C
      END IF
      END
