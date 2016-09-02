      PROGRAM P23
C
C      PROGRAM 2.3 GAUSSIAN ELIMINATION
C      USING A=L*D*LT FACTORISATION LT OVERWRITES A
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),D(IN),B(IN)
C
      READ (5,*) N
      READ (5,*) ((A(I,J),J=1,N),I=1,N)
      READ (5,*) (B(I),I=1,N)
      WRITE(6,*) ('** GAUSSIAN ELIMINATION USING LDLT FACTORISATION **')
      WRITE(6,*) ('************ SYMMETRICAL EQUATIONS ****************')
      WRITE(6,*)
      WRITE(6,*) ('COEFFICIENT MATRIX')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('RIGHT HAND SIDE VECTOR')
      CALL PRINTV(B,N,6)
      WRITE(6,*)
      CALL LDLT(A,IN,D,N)
      WRITE(6,*) ('UT MATRIX')
      DO 1 I=1,N
    1 WRITE(6,100) (A(I,J),J=1,I)
      WRITE(6,*)
      WRITE(6,*) ('DIAGONAL TERMS')
      CALL PRINTV(D,N,6)
      CALL LDLFOR(A,IN,B,N)
      DO 2 I = 1,N
          DO 2 J = 1,N
    2 A(I,J) = A(I,J)/D(I)
      CALL SUBBAC(A,IN,B,N)
      WRITE(6,*)
      WRITE(6,*) ('SOLUTION VECTOR')
      CALL PRINTV(B,N,6)
  100 FORMAT(5E12.4)
      STOP
      END
