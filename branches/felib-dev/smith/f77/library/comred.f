      SUBROUTINE COMRED(BK,L,IW)
C
C      THIS SUBROUTINE REDUCES THE COMPLEX STIFFNESS MATRIX
C
      COMPLEX BK(*),SUM
      KB=IW+1
      DO 1 I=2,L
      IL1=I-1
      KBL=IL1+KB
      IF(KBL-L)3,3,2
    2 KBL=L
    3 DO 1 J=I,KBL
      IJ=(J-I)*L+I
      SUM=BK(IJ)
      NKB=J-KB+1
      IF(NKB)4,4,5
    4 NKB=1
    5 IF(NKB-IL1)6,6,8
    6 DO 7 N=NKB,IL1
      NI=(I-N)*L+N
      NJ=(J-N)*L+N
    7 SUM=SUM-BK(NI)*BK(NJ)/BK(N)
    8 BK(IJ)=SUM
    1 CONTINUE
      RETURN
      END
