      SUBROUTINE COMBAC(BK,R,L,IW)
C
C      THIS SUBROUTINE PERFORMS THE BACK-SUBSTITUTION OF THE
C      COMPLEX STIFFNESS EQUATIONS
C
      COMPLEX BK(*),R(*),SUM
      KB=IW+1
      R(1)=R(1)/BK(1)
      DO 1 I=2,L
      SUM=R(I)
      I1=I-1
      NKB=I-KB+1
      IF(NKB)2,2,3
    2 NKB=1
    3 DO 4 K=NKB,I1
      JN=(I-K)*L+K
      SUM=SUM-BK(JN)*R(K)
    4 CONTINUE
      R(I)=SUM/BK(I)
    1 CONTINUE
      DO 5 JJ=2,L
      I=L-JJ+1
      SUM=0.0
      I1=I+1
      NKB=I-1+KB
      IF(NKB-L)7,7,6
    6 NKB=L
    7 DO 8 K=I1,NKB
      JN=(K-I)*L+I
    8 SUM=SUM+BK(JN)*R(K)
      R(I)=R(I)-SUM/BK(I)
    5 CONTINUE
      RETURN
      END
