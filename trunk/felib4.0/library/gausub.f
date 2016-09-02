C
      SUBROUTINE GAUSUB(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,R,IR,
     *                  ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      GAUSUB calculates the solution of a set of unsymmetric real
C      banded linear equations with a single rhs.  The banded
C      matrix has previously been decomposed into triangular
C      matrices using GAURDN
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
C      Release 1.1  29 Oct 1979 (CG)
C      Commented    11 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of dimensions (IA, JA).  On entry, contains
C              the elements of the band matrix in LU form,
C              after processing by GAURDN
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. MIN(2*HBAND-1, N))
C      IAL     first dimension of AL (.GE. N)
C      JAL     second dimension of AL (.GE. HBAND-1)
C      N       order of banded matrix A
C      HBAND   semi-bandwidth of A
C      ROPIV   vector of dimension IROPIV.  Contains details
C              of row interchanges performed by GAURDN
C      IROPIV  dimension of ROPIV (.GE. N)
C      R       on entry, contains the vector of the rhs,
C              length IR
C      IR      dimension of R (.GE. N)
C      ITEST   error checking option
C
C ARGUMENTS out
C      R       on exit, contains solution vector
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE GAUSUB(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,R,IR,
C    *                  ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IAL,IERROR,II,IK,IR,IROPIV,ITEST,IW,J,
     *        JA,JAL,K,KK,M,N,ROPIV
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,AL,R,X,Y
      DIMENSION A(IA,JA),AL(IAL,JAL),R(IR),ROPIV(IROPIV)
C
      INTRINSIC MIN
      EXTERNAL ERRMES
C
      DATA SRNAME/'GAUSUB'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IR.LT.N) IERROR = 5
         IF (IA.LT.N .OR. JA.LT.2*HBAND-1) IERROR = 4
         IF (IAL.LT.N .OR. JAL.LT.HBAND) IERROR = 3
         IF (IROPIV.LT.N) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      IW = MIN(N,2*HBAND-1)
      M = HBAND - 1
      DO 1010 K = 1,N
         M = MIN(M+1,N)
         J = ROPIV(K)
         IF (J.NE.K) THEN
            X = R(K)
C
C     Rows K and J interchanged
C
            R(K) = R(J)
            R(J) = X
         END IF
         IK = K + 1
         IF (IK.GT.M) THEN
            GO TO 1020
         ELSE
            X = R(K)
            DO 1000 I = IK,M
               II = I - K
               R(I) = R(I) - X*AL(K,II)
 1000       CONTINUE
C
C     Forward substitution complete
C
         END IF
 1010 CONTINUE
C
 1020 CONTINUE
      DO 1040 K = 1,N
         M = MIN(K,IW)
         I = N + 1 - K
         II = I - 1
         Y = A(I,1)
         X = R(I)
         IF (M.NE.1) THEN
            DO 1030 J = 2,M
               KK = J + II
               X = X - A(I,J)*R(KK)
 1030       CONTINUE
         END IF
C
C     Backward substitution complete
C
         R(I) = X*Y
 1040 CONTINUE
C
      END
