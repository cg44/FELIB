      PROGRAM P410
C
C      PROGRAM 4.10 CHARACTERISTIC POLYNOMIAL METHOD FOR
C                   EIGENVALUES OF SYMMETRIC TRIDIAGONAL
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL ALPHA(IN),BETA(IN),DET(0:IN)
C
      READ (5,*) N,(ALPHA(I),I=1,N)
      READ (5,*) (BETA(I),I=1,N-1)
      READ (5,*) J,AL,ALMAX,TOL,ITS
      WRITE(6,*) ('**** CHARACTERISTIC POLYNOMIAL METHOD FOR *******')
      WRITE(6,*) ('* EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX *')
      WRITE(6,*)
      WRITE(6,*) ('MAIN DIAGONAL')
      CALL PRINTV(ALPHA,N,6)
      WRITE(6,*)
      WRITE(6,*) ('OFF-DIAGONAL')
      CALL PRINTV(BETA,N-1,6)
      WRITE(6,*)
      WRITE(6,*) ('EIGENVALUE REQUIRED, 1=LARGEST, 2=NEXT,.. ETC.')
      WRITE(6,1000) J
      WRITE(6,*)
      DET(0) = 1.0
      AOLD = ALMAX
      WRITE(6,*) (' EIGENVALUE   DETERMINANT  NUMBER OF ROOTS SMALLER')
      WRITE(6,*) ('                             THAN CURRENT VALUE   ')  
      ITERS = 0
   10 ITERS = ITERS + 1
          DET(1) = ALPHA(1) - AL
          DO 1 I = 2,N
    1     DET(I) = (ALPHA(I)-AL)*DET(I-1) - BETA(I-1)*BETA(I-1)*DET(I-2)
          NUMBER = 0
          DO 2 I = 1,N
              IF (DET(I).EQ.0.0) GO TO 2
              IF (DET(I-1).EQ.0.0) THEN
                  SIGN = DET(I)*DET(I-2)
              ELSE
                  SIGN = DET(I)*DET(I-1)
              END IF
              IF (SIGN.LT.0.0) NUMBER = NUMBER + 1
    2     CONTINUE
          IF (NUMBER.LE.N-J) THEN
              OLDL = AL
              AL = 0.5* (AL+ALMAX)
          ELSE
              ALMAX = AL
              AL = 0.5* (OLDL+AL)
          END IF
          WRITE (6,100) AL,DET(N),NUMBER
          CALL CHECK(AL,AOLD,TOL,ICON)
          IF(ITERS.NE.ITS.AND.ICON.EQ.0)GOTO 10
      WRITE(6,*)
      WRITE(6,*) ('ITERATIONS TO CONVERGENCE')
      WRITE(6,1000) ITERS
  100 FORMAT (2E12.4,I15)
 1000 FORMAT (I5)
      STOP
      END
