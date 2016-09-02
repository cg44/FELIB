      SUBROUTINE LBBT(L,BT,IBT,IW,N)
C
C      THIS SUBROUTINE SOLVES LA=BT
C
      REAL L(IBT,*),BT(IBT,*)
      DO 1 J=1,N
      BT(J,1)=BT(1,J)/L(1,IW+1)
      IF(J.EQ.1)GO TO 1
      DO 2 I=2,J
      C=0.0
      IP=I-1
      DO 3 M=1,IP
      IF(M.GT.IW)GO TO 3
    4 C=C-BT(J,I-M)*L(I,IW-M+1)
    3 CONTINUE
      C=C+BT(I,J)
      C=C/L(I,IW+1)
    2 BT(J,I)=C
    1 CONTINUE
      RETURN
      END
