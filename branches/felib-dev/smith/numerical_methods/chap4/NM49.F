      PROGRAM P49
C
C      PROGRAM 4.9 LANCZOS REDUCTION OF A SYMMETRIC
C                  MATRIX TO TRIDIAGONAL FORM
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL A(IN,IN),ALPHA(IN),BETA(0:IN-1),V(IN),Y0(IN),Y1(IN),Z(IN)
C
      READ (5,*) N,((A(I,J),J=I,N),I=1,N)
      DO 10 I = 1,N
          DO 10 J = I,N
   10 A(J,I) = A(I,J)
      READ (5,*) (Y1(I),I=1,N)
      WRITE(6,*) ('**** LANCZOS REDUCTION OF A SYMMETRIC MATRIX ****')
      WRITE(6,*) ('************* TO TRIDIAGONAL FORM ***************')
      WRITE(6,*)
      WRITE(6,*) ('MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      CALL NULVEC(Y0,N)
      BETA(0) = 0.0
      DO 1 J = 1,N
          CALL MVMULT(A,IN,Y1,N,N,V)
          DO 2 I = 1,N
    2     V(I) = V(I) - BETA(J-1)*Y0(I)
          CALL VECCOP(Y1,Y0,N)
          CALL VDOTV(Y1,V,ALPHA(J),N)
          IF (J.EQ.N) GO TO 1
          DO 3 I = 1,N
    3     Z(I) = V(I) - ALPHA(J)*Y1(I)
          CALL VDOTV(Z,Z,BETA(J),N)
          BETA(J) = SQRT(BETA(J))
          DO 4 I = 1,N
    4     Y1(I) = Z(I)/BETA(J)
    1 CONTINUE
      WRITE(6,*) ('TRANSFORMED MAIN DIAGONAL OF MATRIX A')
      CALL PRINTV(ALPHA,N,6)
      WRITE(6,*)
      WRITE(6,*) ('TRANSFORMED OFF-DIAGONAL OF MATRIX A')
      WRITE(6,100) (BETA(I),I=1,N-1)
  100 FORMAT(1X,6E12.4)
      STOP
      END
