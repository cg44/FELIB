      SUBROUTINE GAUSBA(PB,IPB,WORK,IWORK,N,IW)
C
C      THIS SUBROUTINE PERFORMS GAUSSIAN REDUCTION OF AN
C      UNSYMMETRIC BANDED MATRIX 'PB' .  ARRAY 'WORK'
C      USED AS WORKING SPACE.
C
      REAL PB(IPB,*),WORK(IWORK,*)
      IWP1=IW+1
      IQ=2*IWP1-1
      IQP=IWP1
      IWP11=IWP1-1
      DO 1 I=1,IWP11
      DO 1 J=1,IQ
      IF(J.GE.IWP1+I)GOTO 2
      PB(I,J)=PB(I,J+IWP1-I)
      GOTO 1
    2 PB(I,J)=0.
      PB(N-I+1,J)=0.
    1 CONTINUE
      DO 3 K=1,N
      L=K+IWP1-1
      IF(L.GT.N)L=N
      IP=0
      S=1.E-10
      DO 4 I=K,L
      IF(ABS(PB(I,1)).LE.S)GOTO 4
      S=ABS(PB(I,1))
      IP=I
    4 CONTINUE
      IF(IP.EQ.0)GOTO 5
      IF(K.EQ.N)GOTO 11
      WORK(IWP1,K)=IP
      IQP=IQP-1
      J=IWP1+IP-K
      IF(IQP.LT.J)IQP=J
      IF(J.EQ.IWP1)GOTO 6
      DO 7 J=1,IQP
      S=PB(K,J)
      PB(K,J)=PB(IP,J)
      PB(IP,J)=S
    7 CONTINUE
    6 K1=K+1
      DO 8 I=K1,L
      S=PB(I,1)/PB(K,1)
      DO 9 J=2,IQ
      IF(J.GT.IQP)GOTO 10
      PB(I,J-1)=PB(I,J)-S*PB(K,J)
      GOTO 9
   10 PB(I,J-1)=PB(I,J)
    9 CONTINUE
      PB(I,IQ)=0.
      WORK(I-K,K)=S
    8 CONTINUE
    3 CONTINUE
    5 WRITE(6,*)' SINGULAR'
   11 RETURN
      END
