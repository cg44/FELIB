      SUBROUTINE CHOBAC(KB,IKB,LOADS,N,IW)
C
C      THIS SUBROUTINE PERFORMS THE CHOLESKI BACK-SUBSTITUTION
C
      REAL KB(IKB,*),LOADS(*)
      LOADS(1)=LOADS(1)/KB(1,IW+1)
      DO 1 I=2,N
      X=0.0
      K=1
      IF(I.LE.IW+1)K=IW-I+2
      DO 2 J=K,IW
    2 X=X+KB(I,J)*LOADS(I+J-IW-1)
    1 LOADS(I)=(LOADS(I)-X)/KB(I,IW+1)
      LOADS(N)=LOADS(N)/KB(N,IW+1)
      I=N-1
    3 X=0.0
      L=I+IW
      IF(I.GT.N-IW)L=N
      M=I+1
      DO 4 J=M,L
    4 X=X+KB(J,IW+I-J+1)*LOADS(J)
      LOADS(I)=(LOADS(I)-X)/KB(I,IW+1)
      I=I-1
      IF(I)5,5,3
    5 CONTINUE
      RETURN
      END
