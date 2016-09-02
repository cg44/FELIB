      SUBROUTINE SPABAC(A,B,N,KDIAG)
C
C      THIS SUBROUTINE PERFORMS THE CHOLESKI BACK-SUBSTITUTION
C      ON THE VARIABLE BANDWIDTH STIFFNESS MATRIX
C
      REAL A(*),B(*)
      INTEGER KDIAG(*)
      B(1)=B(1)/A(1)
      DO 1 I=2,N
      KI=KDIAG(I)-I
      L=KDIAG(I-1)-KI+1
      X=B(I)
      IF(L.EQ.I)GOTO 1
      M=I-1
      DO 2 J=L,M
    2 X=X-A(KI+J)*B(J)
    1 B(I)=X/A(KI+I)
      DO 3 IT=2,N
      I=N+2-IT
      KI=KDIAG(I)-I
      X=B(I)/A(KI+I)
      B(I)=X
      L=KDIAG(I-1)-KI+1
      IF(L.EQ.I)GOTO 3
      M=I-1
      DO 4 K=L,M
    4 B(K)=B(K)-X*A(KI+K)
    3 CONTINUE
      B(1)=B(1)/A(1)
      RETURN
      END
