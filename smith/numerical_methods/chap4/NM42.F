      PROGRAM P42
C
C        SHIFTED ITERATION USING THE POWER METHOD
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),X0(IN),X1(IN)
C
      READ (5,*) N,((A(I,J),J=1,N),I=1,N),SHIFT
      WRITE(6,*) ('*************** SHIFTED ITERATION ***************')
      WRITE(6,*)
      WRITE(6,*) ('MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('SHIFT')
      WRITE(6,1000) SHIFT
      WRITE(6,*)
      DO 10 I = 1,N
   10 A(I,I) = A(I,I) - SHIFT
      READ (5,*) (X0(I),I=1,N)
      READ (5,*) TOL,ITS
      ITERS = 0
    1 CALL MVMULT(A,IN,X0,N,N,X1)
      BIG = 0.0
      DO 2 I = 1,N
    2 IF (ABS(X1(I)).GT.ABS(BIG)) BIG = X1(I)
      DO 3 I = 1,N
    3 X1(I) = X1(I)/BIG
      CALL CHECON(X1,X0,N,TOL,ICON)
      ITERS = ITERS + 1
      IF (ICON.EQ.0 .AND. ITERS.LT.ITS) GO TO 1
      SUM = 0.0
      DO 4 I = 1,N
    4 SUM = SUM + X1(I)**2
      SUM = SQRT(SUM)
      DO 5 I = 1,N
    5 X1(I) = X1(I)/SUM
      WRITE(6,*) ('EIGENVALUE')
      WRITE (6,1000) BIG + SHIFT
      WRITE(6,*)
      WRITE(6,*) ('CORRESPONDING EIGENVECTOR')
      CALL PRINTV(X1,N,6)
      WRITE(6,*)
      WRITE(6,*) ('ITERATIONS TO CONVERGENCE')
      WRITE (6,2000) ITERS
 1000 FORMAT (E12.4)
 2000 FORMAT (I5)
      STOP
      END
