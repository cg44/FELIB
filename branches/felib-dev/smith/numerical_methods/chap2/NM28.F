      PROGRAM P28
C
C      PROGRAM 2.8 JACOBI ITERATION
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),LPU(IN,IN),B(IN),X(IN),XNEW(IN),TEMP(IN)
C
      READ (5,*) N,((A(I,J),J=1,N),I=1,N)
      READ (5,*) (B(I),I=1,N)
      READ (5,*) (X(I),I=1,N)
      READ (5,*) TOL,ITS
      WRITE(6,*) ('*************** JACOBI ITERATION *****************')
      WRITE(6,*)
      WRITE(6,*) ('COEFFICIENT MATRIX')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('RIGHT HAND SIDE VECTOR')
      CALL PRINTV(B,N,6)
      WRITE(6,*)
      DO 1 I = 1,N
          DIAG = A(I,I)
          DO 2 J = 1,N
    2     A(I,J) = A(I,J)/DIAG
    1 B(I) = B(I)/DIAG
      CALL NULL(LPU,IN,N,N)
      DO 3 I = 1,N
          DO 3 J = 1,N
    3 IF (I.NE.J) LPU(I,J) = -A(I,J)
      WRITE(6,*) ('FIRST FEW ITERATIONS')
      ITERS = 0
    4 ITERS = ITERS + 1
      CALL MVMULT(LPU,IN,X,N,N,TEMP)
      CALL VECADD(B,TEMP,XNEW,N)
      CALL CHECON(XNEW,X,N,TOL,ICON)
      IF (ITERS.LE.5) CALL PRINTV(X,N,6)
      IF (ICON.EQ.0.AND.ITERS.LT.ITS) GO TO 4 
      WRITE(6,*)
      WRITE(6,*) ('ITERATIONS TO CONVERGENCE')
      WRITE(6,*) ITERS
      WRITE(6,*)
      WRITE(6,*) ('SOLUTION VECTOR')
      CALL PRINTV(X,N,6)
      STOP
      END
