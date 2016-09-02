C
      SUBROUTINE GAURDN(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      GAURDN decomposes A real unsymmetric matrix of order N into
C      triangular matrices using gaussian elimination with
C      partial pivoting
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
C      Commented    10 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of dimensions (IA, JA).  on entry, contains
C              the elements of the band matrix
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. MIN(2*HBAND-1,N))
C      IAL     first dimension of array AL (.GE. N)
C      JAL     second dimension of array AL (.GE. HBAND-1)
C      N       order of matrix A
C      HBAND   semi-bandwidth of matrix A
C      IROPIV  dimension of vector ROPIV
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       array of dimensions (IA, JA).  on exit, contains
C              the upper triangular matrix U, with the
C              diagonal elements of U stored as their
C              reciprocals.  the i'th row of U is stored in the
C              i'th row of A, starting with the diagonal
C              element of U in A(I, 1)
C      ROPIV   vector of length IROPIV, containing details of
C              row interchanges.  if no interchange occurs at
C              the r'th major step then ROPIV(R)=R; if the R
C              and J rows are interchanged then ROPIV(R)=J
C
C ROUTINES called
C      VEPS    ERRMES
C
C     SUBROUTINE GAURDN(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IAL,IERROR,IK,IR,IROPIV,ITEST,IW,J,JA,
     *        JAL,JJ,JR,JTEST,K,L,M,N,ROPIV
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,AL,EPS,ONE,VEPS,X,Y,ZERO
      DIMENSION A(IA,JA),AL(IAL,JAL),ROPIV(IROPIV)
C
*      INTRINSIC MIN
*      EXTERNAL ERRMES, VEPS
C
      DATA ONE/1.0D0/,SRNAME/'GAURDN'/,ZERO/0.0D0/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IA.LT.N .OR. JA.LT.2*HBAND-1) IERROR = 4
         IF (IAL.LT.N .OR. JAL.LT.HBAND) IERROR = 3
         IF (IROPIV.LT.N) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      IERROR = 5
      EPS = VEPS(X)
      IW = MIN(N,2*HBAND-1)
      M = HBAND
      K = IW - HBAND
      IF (K.GT.0) THEN
         DO 1020 I = 1,K
            L = IW - M
            DO 1000 J = 1,M
               JJ = J + L
               A(I,J) = A(I,JJ)
 1000       CONTINUE
            M = M + 1
            DO 1010 J = M,IW
               A(I,J) = ZERO
 1010       CONTINUE
 1020    CONTINUE
      END IF
      M = N - IW + HBAND + 1
      J = IW + 1
      IF (M.LE.N) THEN
         DO 1040 I = M,N
            J = J - 1
            DO 1030 K = J,IW
               A(I,K) = ZERO
 1030       CONTINUE
 1040    CONTINUE
C
C     Zeros inserted
C
      END IF
      DO 1060 I = 1,N
         X = ZERO
         DO 1050 J = 1,IW
            X = X + DABS(A(I,J))
 1050    CONTINUE
         IF (X.GT.ZERO) THEN
C
C     Row norms of A calculated and their reciprocals stored in
C     first column of AL
C
            AL(I,1) = ONE/X
         ELSE
            GO TO 1120
         END IF
 1060 CONTINUE
C
      IERROR = 6
      DO 1110 IR = 1,N
         X = ZERO
         M = MIN0(IR+HBAND-1,N)
         DO 1070 I = IR,M
            Y = DABS(A(I,1))*AL(I,1)
            IF (Y.GT.X) THEN
               X = Y
               J = I
            END IF
 1070    CONTINUE
C
C     IR'TH pivot element selected.
C
         IF (X.LT.EPS) THEN
            GO TO 1130
         ELSE
            ROPIV(IR) = J
            IF (J.NE.IR) THEN
               DO 1080 I = 1,IW
                  X = A(IR,I)
                  A(IR,I) = A(J,I)
                  A(J,I) = X
 1080          CONTINUE
C
C     Row IR and J interchanged.
C
C
               AL(J,1) = AL(IR,1)
            END IF
            JR = IR + 1
            Y = ONE/A(IR,1)
            IF (JR.LE.M) THEN
               DO 1100 I = JR,M
                  X = A(I,1)*Y
                  IF (IW.GE.2) THEN
                     DO 1090 J = 2,IW
                        A(I,J-1) = A(I,J) - X*A(IR,J)
 1090                CONTINUE
                  END IF
                  IK = I - IR
                  AL(IR,IK) = X
                  A(I,IW) = ZERO
 1100          CONTINUE
            END IF
C
C     Elimination completed
C
            A(IR,1) = Y
         END IF
 1110 CONTINUE
C
      RETURN
 1120 CONTINUE
      IR = I
C
 1130 CONTINUE
      A(IR,1) = ZERO
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
C
      END
