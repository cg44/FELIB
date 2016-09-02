      SUBROUTINE FDIAGV(BK,KM,IKM,G,KDIAG,IDOF)
C
C      THIS SUBROUTINE ASSEMBLES A DIAGONAL ELEMENT MATRIX
C      INTO THE GLOBAL SYSTEM
C
      REAL BK(*),KM(IKM,*)
      INTEGER G(*),KDIAG(*)
      DO 1 I=1,IDOF
      J=G(I)
      IF(J.EQ.0)GO  TO 1
      K=KDIAG(J)
      BK(K)=BK(K)+KM(I,I)
    1 CONTINUE
      RETURN
      END