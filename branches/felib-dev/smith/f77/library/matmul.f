      SUBROUTINE MATMUL(A,IA,B,IB,C,IC,L,M,N)
C
C      THIS SUBROUTINE FORMS THE PRODUCT OF TWO MATRICES
C
      REAL A(IA,*),B(IB,*),C(IC,*)
      DO 1 I=1,L
      DO 1 J=1,N
      X=0.0
      DO 2 K=1,M
    2 X=X+A(I,K)*B(K,J)
      C(I,J)=X
    1 CONTINUE
      RETURN
      END
