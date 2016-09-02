      SUBROUTINE MATINV(A,IA,N)
C
C      THIS SUBROUTINE FORMS THE INVERSE OF A MATRIX
C      USING GAUSS-JORDAN TRANSFORMATION
C
      REAL A(IA,*)
      DO 1 K=1,N
      CON=A(K,K)
      A(K,K)=1.
      DO 2 J=1,N
    2 A(K,J)=A(K,J)/CON
      DO 1 I=1,N
      IF(I.EQ.K)GOTO 1
      CON=A(I,K)
      A(I,K)=0.
      DO 3 J=1,N
    3 A(I,J)=A(I,J)-A(K,J)*CON
    1 CONTINUE
      RETURN
      END
