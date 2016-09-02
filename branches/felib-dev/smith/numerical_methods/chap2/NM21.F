      PROGRAM P21
C
C      PROGRAM 2.1 GAUSSIAN ELIMINATION
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),B(IN)
C
      READ (5,*) N
      READ (5,*) ((A(I,J),J=1,N),I=1,N)
      READ (5,*) (B(I),I=1,N)
      WRITE(6,*) ('************** GAUSSIAN ELIMINATION ************')
      WRITE(6,*)
      WRITE(6,*) ('COEFFICIENT MATRIX')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('RIGHT HAND SIDE VECTOR')
      CALL PRINTV(B,N,6)
      WRITE(6,*)
C      CONVERT TO UPPER TRIANGULAR FORM
      DO 1 K = 1,N - 1
          IF (ABS(A(K,K)).GT.1.E-6) THEN
              DO 2 I = K + 1,N
                  X = A(I,K)/A(K,K)
                  DO 3 J = K + 1,N
                      A(I,J) = A(I,J) - A(K,J)*X
    3             CONTINUE
                  B(I) = B(I) - B(K)*X
    2         CONTINUE
          ELSE
              WRITE (6,*) ('ZERO PIVOT FOUND IN LINE:')
              WRITE (6,*) K
              STOP
          END IF
    1 CONTINUE
      WRITE(6,*) ('MODIFIED MATRIX')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('MODIFIED RIGHT HAND SIDE VECTOR')
      CALL PRINTV(B,N,6)
      WRITE(6,*)
C          BACK SUBSTITUTION
      DO 5 I = N,1,-1
          SUM = B(I)
          IF (I.LT.N) THEN
              DO 6 J = I + 1,N
                  SUM = SUM - A(I,J)*B(J)
    6         CONTINUE
          END IF
          B(I) = SUM/A(I,I)
    5 CONTINUE
C          PRINT THE RESULTS
      WRITE(6,*) ('SOLUTION VECTOR')
      CALL PRINTV(B,N,6)
      STOP
      END

