      SUBROUTINE BACSUB(BK,LOADS,N,IW)
C
C      THIS SUBROUTINE PERFORMS THE GAUSSIAN BACK-SUBSTITUTION
C
      REAL BK(*),LOADS(*)
      LOADS(1)=LOADS(1)/BK(1)
      DO 1 I=2,N
      SUM=LOADS(I)
      I1=I-1
      NKB=I-IW
      IF(NKB)2,2,3
    2 NKB=1
    3 DO 4 K=NKB,I1
      JN=(I-K)*N+K
      SUM=SUM-BK(JN)*LOADS(K)
    4 CONTINUE
      LOADS(I)=SUM/BK(I)
    1 CONTINUE
      DO 5 JJ=2,N
      I=N-JJ+1
      SUM=0.
      I1=I+1
      NKB=I+IW
      IF(NKB-N)7,7,6
    6 NKB=N
    7 DO 8 K=I1,NKB
      JN=(K-I)*N+I
    8 SUM=SUM+BK(JN)*LOADS(K)
      LOADS(I)=LOADS(I)-SUM/BK(I)
    5 CONTINUE
      RETURN
      END
