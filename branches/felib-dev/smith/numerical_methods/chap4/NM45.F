      PROGRAM P45
C
C      PROGRAM 4.5 CONVERSION OF A*X=LAMBDA*B*X
C                  TO STANDARD FORM
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),D(IN),B(IN,IN),C(IN,IN),E(IN),DD(IN,IN)
C
      READ (5,*) N,((A(I,J),J=1,N),I=1,N)
      READ (5,*) ((B(I,J),J=1,N),I=1,N)
      WRITE(6,*) ('******** CONVERSION OF  A*X=LAMBDA*B*X **********')
      WRITE(6,*) ('******** TO STANDARD SYMMETRICAL FORM ***********')
      WRITE(6,*)
      WRITE(6,*) ('MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('MATRIX B')
      CALL PRINTA(B,IN,N,N,6)
      WRITE(6,*)
      CALL LDLT(B,IN,D,N)
      DO 1 I = 1,N
    1 D(I) = SQRT(D(I))
      DO 2 J = 1,N
          DO 2 I = J,N
    2 B(I,J) = B(I,J)/D(J)
      DO 3 I = 1,N - 1
          DO 3 J = I + 1,N
    3 B(I,J) = B(J,I)
      DO 4 J = 1,N
          DO 5 I = 1,N
    5     E(I) = A(I,J)
          CALL SUBFOR(B,IN,E,N)
          DO 6 I = 1,N
    6     C(I,J) = E(I)
    4 CONTINUE
      WRITE(6,*) ('MATRIX C')
      CALL PRINTA(C,IN,N,N,6)
      WRITE(6,*)
      CALL MATRAN(DD,IN,C,IN,N,N)
      DO 7 J = 1,N
          DO 8 I = 1,N
    8     E(I) = DD(I,J)
          CALL SUBFOR(B,IN,E,N)
          DO 9 I = 1,N
    9     DD(I,J) = E(I)
    7 CONTINUE
      WRITE(6,*) ('FINAL SYMMETRICAL MATRIX DD')
      CALL PRINTA(DD,IN,N,N,6)
      STOP
      END
