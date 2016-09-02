      SUBROUTINE LBKBAN(L,B,K,IK,IW,N)
C
C      THIS SUBROUTINE SOLVES LB=K
C
      REAL L(IK,*),B(IK,*),K(IK,*)
      IP=IW+1
      DO 1 J=1,IP
    1 B(1,J)=K(J,IW+2-J)/L(1,IW+1)
      DO 2 I=2,IP
      IQ=IW+I
      DO 2 J=1,IQ
      IF(J.GT.N)GO TO 2
      X=0.0
      IR=I-1
      DO 3 M=1,IR
    3 X=X+L(I,M-I+IW+1)*B(M,J)
      IF(J.LE.I)Y=K(I,J-I+IW+1)
      IF(J.GT.I)Y=K(J,IW+1-J+I)
      B(I,J)=(Y-X)/L(I,IW+1)
    2 CONTINUE
      IP=IW+2
      IF(IP.GT.N)GO TO 5
      DO 4 I=IP,N
      DO 4 J=1,N
      IF(IW+1-J+I.LT.1)GO TO 4
      X=.0
      DO 6 M=1,IW
    6 X=X+L(I,M)*B(I-IW+M-1,J)
      IF(I-J.LE.IW)GO TO 8
      Y=.0
      GO TO 7
    8 IF(J.GT.I)GO TO 10
      Y=K(I,J-I+IW+1)
      GO TO 7
   10 Y=K(J,IW+1-J+I)
    7 B(I,J)=(Y-X)/L(I,IW+1)
    4 CONTINUE
    5 CONTINUE
      RETURN
      END
