      PROGRAM P43

C
C      PROGRAM 4.3 SHIFTED INVERSE ITERATION FOR EIGENVALUE
C                  CLOSEST TO A GIVEN VALUE
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),UPTRI(IN,IN),LOWTRI(IN,IN),X0(IN),X1(IN)
C
      READ (5,*) N,((A(I,J),J=1,N),I=1,N),SHIFT
      WRITE(6,*) ('********** SHIFTED INVERSE ITERATION ************')
      WRITE(6,*)
      WRITE(6,*) ('MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('SHIFT')
      WRITE(6,1000) SHIFT
      WRITE(6,*)
      DO 10 I=1,N
   10 A(I,I)=A(I,I)-SHIFT
      CALL LUFAC(A,UPTRI,LOWTRI,IN,N)
      READ (5,*) (X0(I),I=1,N)
      READ (5,*) TOL,ITS
      ITERS = 0
      CALL VECCOP(X0,X1,N)
    1 CALL SUBFOR(LOWTRI,IN,X1,N)
      CALL SUBBAC(UPTRI,IN,X1,N)
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
      WRITE(6,*) ('EIGENVALUE CLOSEST TO SHIFT')
      WRITE (6,1000) 1./BIG + SHIFT
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
