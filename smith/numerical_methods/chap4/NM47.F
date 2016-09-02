      PROGRAM P47
C
C      PROGRAM 4.7 HOUSEHOLDER REDUCTION OF A SYMMETRIC
C                  MATRIX TO TRIDIAGONAL FORM
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),A1(IN,IN),P(IN,IN),V(IN)
C
      READ (5,*) N,((A(I,J),J=I,N),I=1,N)
      DO 11 I = 1,N
          DO 11 J = I,N
   11 A(J,I) = A(I,J)
      WRITE(6,*) ('** HOUSEHOLDER REDUCTION OF A SYMMETRIC MATRIX **')
      WRITE(6,*) ('************* TO TRIDIAGONAL FORM ***************')
      WRITE(6,*)
      WRITE(6,*) ('MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      DO 1 K = 1,N - 2
          R = 0.0
          DO 2 L = K,N - 1
    2     R = R + A(K,L+1)*A(K,L+1)
          R = SQRT(R)
          IF (R*A(K,K+1).GT.0.0) R = -R
          H = -1.0/ (R*R-R*A(K,K+1))
          CALL NULVEC(V,N)
          V(K+1) = A(K,K+1) - R
          DO 3 L = K + 2,N
    3     V(L) = A(K,L)
          CALL VVMULT(V,V,P,IN,N,N)
          CALL MSMULT(P,IN,H,N,N)
          DO 4 L = 1,N
    4     P(L,L) = P(L,L) + 1.0
          CALL MATMUL(A,IN,P,IN,A1,IN,N,N,N)
          CALL MATMUL(P,IN,A1,IN,A,IN,N,N,N)
    1 CONTINUE
      WRITE(6,*) ('TRANSFORMED MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      STOP
      END
