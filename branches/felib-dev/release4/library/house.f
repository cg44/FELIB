C
      SUBROUTINE HOUSE(A,IA,JA,T,IT,JT,DIAG,IDIAG,SUB,ISUB,N,TOL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      HOUSE uses householder's method to reduce a real symmetric
C      matrix to tridiagonal form for use with QLVAL, QLVEC
C      A is stored as A full matrix
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
C      Commented    19 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array containing the elements of the symmetric
C              matrix
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. N)
C      IT      first dimension of array T (.GE. N)
C      JT      second dimension of T (.GE. N)
C      IDIAG   dimension of vector DIAG (.GE. N)
C      ISUB    dimension of vector SUB (.GE. N)
C      N       order of matrix A
C      TOL     value of RMIN/EPS, where RMIN is the smallest
C              positive number exactly representable on the
C              computer, and EPS is the smallest positive
C              number such that 1.+EPS>1.
C      ITEST   error checking option
C
C ARGUMENTS out
C      T       contains the orthogonal matrix 'q', the product
C              of the householder transformation matrices
C      DIAG    contains diagonal elements of tridiagonal matrix
C      SUB     contains the N-1 off-diagonal elements of the
C              tridiagonal matrix stored in SUB(2) to SUB(N).
C              SUB(1)=0.
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE HOUSE(A,IA,JA,T,IT,JT,DIAG,IDIAG,SUB,ISUB,N,TOL,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IDIAG,IERROR,II,ISUB,IT,ITEST,J,J1,JA,JT,K,L,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,DIAG,F,G,H,HH,SUB,T,TOL,VEPS,X
      DIMENSION A(IA,JA),DIAG(IDIAG),SUB(ISUB),T(IT,JT)
C
      INTRINSIC ABS
      EXTERNAL ERRMES,VEPS
C
      DATA SRNAME/'HOUSE'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IDIAG.LT.N .OR. ISUB.LT.N) IERROR = 2
         IF (IA.LT.N .OR. JA.LT.N) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1010 I = 1,N
         DO 1000 J = 1,I
            T(I,J) = A(I,J)
 1000    CONTINUE
 1010 CONTINUE
C
      IF (N.NE.1) THEN
         DO 1080 II = 2,N
            I = N - II + 2
            L = I - 2
            F = T(I,I-1)
            G = 0.0D0
            IF (L.NE.0) THEN
               DO 1020 K = 1,L
                  G = G + T(I,K)*T(I,K)
 1020          CONTINUE
C
C     If G is too small for orthogonality to be guaranteed the
C     transformation is skipped
C
            END IF
            H = G + F*F
            IF (G.GT.TOL) THEN
               L = L + 1
               G = DSQRT(H)
               IF (F.GE.0.0D0) G = -G
               SUB(I) = G
               H = H - F*G
               T(I,I-1) = F - G
               F = 0.0D0
               DO 1050 J = 1,L
C
C     Form element of A*U
C
                  T(J,I) = T(I,J)/H
                  G = 0.0D0
                  DO 1030 K = 1,J
                     G = G + T(J,K)*T(I,K)
 1030             CONTINUE
                  J1 = J + 1
                  IF (J1.LE.L) THEN
                     DO 1040 K = J1,L
C
C     Form element of p
C
                        G = G + T(K,J)*T(I,K)
 1040                CONTINUE
                  END IF
                  SUB(J) = G/H
C
C     Form K
C
                  F = F + G*T(J,I)
 1050          CONTINUE
C
C     Form reduced A
C
               HH = F/ (H+H)
               DO 1070 J = 1,L
                  F = T(I,J)
                  G = SUB(J) - HH*F
                  SUB(J) = G
                  DO 1060 K = 1,J
                     T(J,K) = T(J,K) - F*SUB(K) - G*T(I,K)
 1060             CONTINUE
 1070          CONTINUE
            ELSE
               SUB(I) = F
               H = 0.0D0
            END IF
            DIAG(I) = H
 1080    CONTINUE
      END IF
C
C     Accumulation of transformation matrices
C
      SUB(1) = 0.0D0
      DIAG(1) = 0.0D0
      DO 1130 I = 1,N
         L = I - 1
         IF (ABS(DIAG(I)).GT.VEPS(X)) THEN
            DO 1110 J = 1,L
               G = 0.0D0
               DO 1090 K = 1,L
                  G = G + T(I,K)*T(K,J)
 1090          CONTINUE
               DO 1100 K = 1,L
                  T(K,J) = T(K,J) - G*T(K,I)
 1100          CONTINUE
 1110       CONTINUE
         END IF
         DIAG(I) = T(I,I)
         T(I,I) = 1.0D0
         IF (L.NE.0) THEN
            DO 1120 J = 1,L
               T(I,J) = 0.0D0
               T(J,I) = 0.0D0
 1120       CONTINUE
         END IF
 1130 CONTINUE
C
      END
