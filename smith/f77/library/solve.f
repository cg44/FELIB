      SUBROUTINE SOLVE(K,IK,U,F,N)
C
C      THIS SUBROUTINE PERFORMS GAUSSIAN ELIMINATION WITH
C      PARTIAL PIVOTING ON A FULL N*N MATRIX
C
      REAL K(IK,*),F(*),U(*)
C
C      PIVOTING STAGE
C
      DO 1,I=1,N-1
      BIG=ABS(K(I,I))
      IHOLD=I
      DO 10 J=I+1,N
      IF(ABS(K(J,I)).GT.BIG)THEN
      BIG=ABS(K(J,I))
      IHOLD=J
      END IF
   10 CONTINUE
      IF(IHOLD.NE.I)THEN
      DO 12 J=I,N
      HOLD=K(I,J)
      K(I,J)=K(IHOLD,J)
      K(IHOLD,J)=HOLD
   12 CONTINUE
      HOLD=F(I)
      F(I)=F(IHOLD)
      F(IHOLD)=HOLD
      END IF
C
C      ELIMINATION STAGE
C
      DO 3 J=I+1,N
      FAC=K(J,I)/K(I,I)
      DO 4 L=I,N
    4 K(J,L)=K(J,L)-K(I,L)*FAC
      F(J)=F(J)-F(I)*FAC
    3 CONTINUE
    1 CONTINUE
C
C      BACK-SUBSTITUTION STAGE
C
      DO 9 I=N,1,-1
      SUM=0.
      DO 6 L=I+1,N
    6 SUM=SUM+K(I,L)*U(L)
      U(I)=(F(I)-SUM)/K(I,I)
    9 CONTINUE
      RETURN
      END
