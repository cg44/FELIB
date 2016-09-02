      PROGRAM P211
C
C      PROGRAM 2.11 STEEPEST DESCENT
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),B(IN),X(IN),R(IN),U(IN),TEMP(IN),XNEW(IN)
C
      READ (5,*) N,((A(I,J),J=1,N),I=1,N)
      READ (5,*) (B(I),I=1,N)
      READ (5,*) (X(I),I=1,N)
      READ (5,*) TOL,ITS
      WRITE(6,*) ('*************** STEEPEST DESCENT *****************')
      WRITE(6,*)
      WRITE(6,*) ('COEFFICIENT MATRIX')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('RIGHT HAND SIDE VECTOR')
      CALL PRINTV(B,N,6)
      WRITE(6,*)
      WRITE(6,*) ('STARTING VECTOR')
      CALL PRINTV(X,N,6)
      WRITE(6,*)
      CALL MVMULT(A,IN,X,N,N,TEMP)
      CALL VECSUB(B,TEMP,R,N)
      WRITE(6,*) ('FIRST FEW ITERATIONS')
      ITERS = 0
    1 ITERS = ITERS + 1
      CALL MVMULT(A,IN,R,N,N,U)
      CALL VDOTV(R,R,UP,N)
      CALL VDOTV(R,U,DOWN,N)
      ALPHA = UP/DOWN
      CALL VECCOP(R,TEMP,N)
      CALL VSMULT(TEMP,ALPHA,N)
      CALL VECADD(X,TEMP,XNEW,N)
      CALL VSMULT(U,ALPHA,N)
      CALL VECSUB(R,U,R,N)
      CALL CHECON(XNEW,X,N,TOL,ICON)
      IF (ITERS.LE.5) CALL PRINTV(X,N,6)
      IF (ICON.EQ.0.AND.ITERS.LT.ITS) GO TO 1
      WRITE(6,*)
      WRITE(6,*) ('ITERATIONS TO CONVERGENCE')
      WRITE(6,*)ITERS
      WRITE(6,*)
      WRITE(6,*) ('SOLUTION VECTOR')
      CALL PRINTV(X,N,6)
      STOP
      END
