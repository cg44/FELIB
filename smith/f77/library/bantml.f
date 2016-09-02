      SUBROUTINE BANTML(KB,IKB,LOADS,ANS,N,IW)
C
C      THIS SUBROUTINE MULTIPLIES AN UNSYMMETRIC BANDED MATRIX
C      'PB' BY THE VECTOR 'LOADS'.
C
      REAL KB(IKB,*),LOADS(*),ANS(*)
      DO 1 I=1,N
      X=0.
      K=IW+2
      L=IW+IW+1
      DO 2 J=1,L
      K=K-1
      IF(I-K+1.GT.N)GOTO 2
      IF(I-K+1.LT.1)GOTO 2
      X=X+KB(I,J)*LOADS(I-K+1)
    2 CONTINUE
      ANS(I)=X
    1 CONTINUE
      RETURN
      END
