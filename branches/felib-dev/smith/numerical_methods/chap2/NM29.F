      PROGRAM P29
C
C      PROGRAM 2.9 GAUSS-SEIDEL ITERATION
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),U(IN,IN),B(IN),X(IN),XNEW(IN)
C
      READ (5,*) N,((A(I,J),J=1,N),I=1,N)
      READ (5,*) (B(I),I=1,N)
      READ (5,*) (X(I),I=1,N)
      READ (5,*) TOL,ITS
      WRITE(6,*) ('*********** GAUSS-SEIDEL ITERATION ***************')
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
      CALL NULL(U,IN,N,N)
      DO 3 I = 1,N
          DO 3 J = I + 1,N
              U(I,J) = -A(I,J)
    3 A(I,J) = 0.0
      WRITE(6,*) ('FIRST FEW ITERATIONS')
      ITERS = 0
    4 ITERS = ITERS + 1
      CALL MVMULT(U,IN,X,N,N,XNEW)
      CALL VECADD(B,XNEW,XNEW,N)
      CALL SUBFOR(A,IN,XNEW,N)
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
