C
      SUBROUTINE QLVEC(DIAG,IDIAG,SUB,ISUB,T,IT,JT,N,EPS,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QLVEC calculates all the eigenvalues and eigenvectors of a
C      real symmetric tridiagonal matrix, or of a full real
C      symmetric matrix reduced to tridiagonal form by HOUSE
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
C      Release 1.1  29 Oct 1979 (ims)
C      Commented    16 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      DIAG    vector of length IDIAG. On entry, contains
C              diagonal elements of tridiagonal matrix
C      IDIAG   dimension of vector DIAG (.GE.N)
C      SUB     vector of length ISUB. On entry, contains
C              SUB-DIAGONAL elements of tridiagonal matrix in
C              SUB(2) to SUB(N). On exit, contents of SUB are
C              destroyed
C      ISUB    dimension of vector SUB (.GE. N)
C      T       array of dimension (IT, JT).  There are two modes
C              of operation:
C                (I)  eigen vectors of tridiagonal matrix
C                     on entry T should contain the identity
C                     matrix of order N
C                (II) eigenvectors of full symmetric matrix
C                     on entry T should contain the orthogonal
C                     matrix obtained from HOUSE
C      IT      first dimension of array T (.GE. N)
C      JT      second dimension of T (.GE. N)
C      N       order of tridiagonal matrix
C      EPS     smallest positive number such that 1.+EPS>1.
C      ITEST   error checking option
C
C ARGUMENTS out
C      DIAG    vector of dimension IDIAG.  on exit, contains
C              the eigenvalues in ascending order
C      T       array of dimension (IT, JT).  on exit, T contains
C              the normalised eigenvectors , with T(I, J),
C              I=1(1)N containing the eigenvector corresponding
C              to the eigenvalue in DIAG(J)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QLVEC(DIAG,IDIAG,SUB,ISUB,T,IT,JT,N,EPS,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,I1,IDIAG,IERROR,II,ISUB,IT,ITEST,J,JT,JTEST,K,L,
     *        M,M1,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION B,C,DIAG,EPS,F,G,H,P,R,S,SUB,T
      DIMENSION DIAG(IDIAG),SUB(ISUB),T(IT,JT)
      DATA SRNAME/'QLVEC'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IT.LT.N .OR. JT.LT.N) IERROR = 3
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
      DO 1070 L = 1,N
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
               GO TO 1110
            ELSE
C
C     Form shift
C
               J = J + 1
               G = DIAG(L)
               H = DIAG(L+1) - G
               IF (DABS(H).GE.DABS(SUB(L))) THEN
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
               DO 1060 II = L,M1
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
C
C     Form vector
C
                  DIAG(I+1) = H + S* (C*G+S*DIAG(I))
                  DO 1050 K = 1,N
                     H = T(K,I+1)
                     T(K,I+1) = S*T(K,I) + C*H
                     T(K,I) = C*T(K,I) - S*H
 1050             CONTINUE
 1060          CONTINUE
C
               SUB(L) = S*P
               DIAG(L) = C*P
               IF (DABS(SUB(L)).GT.B) GO TO 1030
            END IF
         END IF
         DIAG(L) = DIAG(L) + F
 1070 CONTINUE
C
C     Order eigenvalues and eigenvectors
C
      DO 1100 I = 1,N
         K = I
         P = DIAG(I)
         I1 = I + 1
         IF (I1.LE.N) THEN
            DO 1080 J = I1,N
               IF (DIAG(J).LT.P) THEN
                  K = J
                  P = DIAG(J)
               END IF
 1080       CONTINUE
         END IF
C
         IF (K.NE.I) THEN
            DIAG(K) = DIAG(I)
            DIAG(I) = P
            DO 1090 J = 1,N
               P = T(J,I)
               T(J,I) = T(J,K)
               T(J,K) = P
 1090       CONTINUE
         END IF
 1100 CONTINUE
C
      ITEST = 0
      RETURN
C
 1110 CONTINUE
      IF (JTEST.NE.-1) THEN
         IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
      END IF
C
      END
