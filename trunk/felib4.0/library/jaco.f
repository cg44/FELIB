C
      SUBROUTINE JACO(A,IA,JA,DIAG,IDIAG,SUB,ISUB,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      JACO reduces a real symmetric band matrix to tridiagonal form
C      using jacobi rotations, for use with QLVAL or QLVEC.
C      the lower triangle of the matrix is stored in a
C      rectangular array.
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
C      Commented    26 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       contains the elements of the lower triangle of
C              the positive definite band matrix
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. HBAND)
C      IDIAG   dimension of vector DIAG (.GE. N)
C      ISUB    dimension of vector SUB (.GE. N)
C      N       order of matrix A
C      HBAND   semi-bandwidth of matrix A (includes diagonal)
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       destroyed on succesful exit
C      DIAG    contains the diagonal elements of the
C              tridiagonal matrix
C      SUB     contains the N-1 off-diagonal elements of the
C              tridiagonal matrix stored in SUB(2) to SUB(N),
C              with SUB(1)=0
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE JACO(A,IA,JA,DIAG,IDIAG,SUB,ISUB,N,HBAND,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IDIAG,IERROR,IR,IRR,ISUB,ITEST,IUGL,J,
     *        J2,JA,JL,JM,K,KR,L,M,MAXL,MAXR,N,N2
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,B,C,C2,CS,DIAG,G,S,S2,SUB,U,U1,VEPS,X
      DIMENSION A(IA,JA),DIAG(IDIAG),SUB(ISUB)
C
      EXTERNAL ERRMES,VEPS
C
      DATA SRNAME/'JACO'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (ISUB.LT.N) IERROR = 4
         IF (IDIAG.LT.N) IERROR = 3
         IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      G = 0.0D0
      DO 1010 I = 1,N
         DO 1000 J = 1,HBAND
            L = HBAND + 1 - J
            K = I + 1 - L
            IF (K.GT.0) A(K,L) = A(I,J)
 1000    CONTINUE
 1010 CONTINUE
C
      KR = HBAND - 1
      DO 1030 I = 1,KR
         M = HBAND - I
         DO 1020 J = 1,M
            K = N + 1 - I
            L = HBAND + 1 - J
            A(K,L) = 0.0D0
 1020    CONTINUE
 1030 CONTINUE
C
      M = HBAND - 1
      N2 = N - 2
      IF (N2.GE.1) THEN
         DO 1080 K = 1,N2
            MAXR = M
            IF (N-K.LT.M) MAXR = N - K
            DO 1070 IRR = 2,MAXR
               IR = 2 + MAXR - IRR
               KR = K + IR
               DO 1060 J = KR,N,M
                  IF (J.EQ.KR) THEN
                     IF (ABS(A(K,IR+1)).LT.VEPS(X)) THEN
                        GO TO 1070
                     ELSE
                        B = -A(K,IR)/A(K,IR+1)
                        IUGL = K
                     END IF
                  ELSE IF (ABS(G).LT.VEPS(X)) THEN
                     GO TO 1070
                  ELSE
                     JM = J - M
                     B = -A(JM-1,M+1)/G
                     IUGL = J - M
                  END IF
                  S = 1.0D0/DSQRT(1.0D0+B*B)
                  C = B*S
                  C2 = C*C
                  S2 = S*S
                  CS = C*S
                  U = C2*A(J-1,1) - 2.0D0*CS*A(J-1,2) + S2*A(J,1)
                  U1 = S2*A(J-1,1) + 2.0D0*CS*A(J-1,2) + C2*A(J,1)
                  A(J-1,2) = CS* (A(J-1,1)-A(J,1)) + (C2-S2)*A(J-1,2)
                  A(J-1,1) = U
                  A(J,1) = U1
                  J2 = J - 2
                  DO 1040 L = IUGL,J2
                     JL = J - L
                     U = C*A(L,JL) - S*A(L,JL+1)
                     A(L,JL+1) = S*A(L,JL) + C*A(L,JL+1)
                     A(L,JL) = U
 1040             CONTINUE
                  JM = J - M
                  IF (J.NE.KR) A(JM-1,M+1) = C*A(JM-1,M+1) - S*G
                  MAXL = M - 1
                  IF (N-J.LT.M-1) MAXL = N - J
                  IF (MAXL.GT.0) THEN
                     DO 1050 L = 1,MAXL
                        U = C*A(J-1,L+2) - S*A(J,L+1)
                        A(J,L+1) = S*A(J-1,L+2) + C*A(J,L+1)
                        A(J-1,L+2) = U
 1050                CONTINUE
                  END IF
                  IF (J+M.LE.N) THEN
                     G = -S*A(J,M+1)
                     A(J,M+1) = C*A(J,M+1)
                  END IF
 1060          CONTINUE
 1070       CONTINUE
 1080    CONTINUE
      END IF
C
      SUB(1) = 0.0D0
      DO 1090 I = 1,N
         DIAG(I) = A(I,1)
 1090 CONTINUE
C
      IF (2.LE.N) THEN
         DO 1100 I = 2,N
            SUB(I) = A(I-1,2)
 1100    CONTINUE
      END IF
*     $st$ unreachable comments ...
C
C
      END
