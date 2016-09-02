      SUBROUTINE EA15GD(J,LALFA,ALFA,BETA,E,MATCH,Z)
C FIND AN APPROXIMATE EIGENVECTOR OF A LANCZOS TRIDIAGONAL MATRIX BY
C     FORWARD AND BACKWARD RECURRENCE.
C J IS THE ORDER OF THE MATRIX AND IS NOT ALTERED.
C ALFA,BETA MUST BE EXACTLY AS LEFT BY LANCZ1. THEY ARE NOT ALTERED.
C E IS THE EIGENVALUE CORRESPONDING TO WHICH AN EIGENVECTOR IS WANTED.
C     IT IS NOT ALTERED.
C MATCH IS THE MATCHING POINT BETWEEN THE RECURRENCES. IT IS NOT
C     ALTERED.
C Z IS SET TO THE REQUIRED EIGENVECTOR.
      REAL             ALFA(LALFA),BETA(LALFA),E,Z(J)
      REAL             PSI,W,W0,W1,W2,Z1,Z2
      REAL             ZERO,ONE
      DATA ZERO/0.0E0/, ONE/1.0E0/
      J1=J-1
C
C FORWARD RECURRENCE.
      W1=ZERO
      W2=ONE
      W=ALFA(1)-E
      IF(J1.LE.0)GO TO 15
      DO 10 K=1,J1
      W0=W1
      W1=W2
      W2=-W/BETA(K+1)
      W=(ALFA(K+1)-E)*W2+BETA(K+1)*W1
      Z(K)=W1
      IF(K.EQ.MATCH)GO TO 20
   10 CONTINUE
   15 K=J
      W0=W1
      W1=W2
      Z(K)=W1
C BACKWARD RECURRENCE
   20 Z1=ONE
      PSI=ALFA(J)-E
      IF(K.GT.J1)GO TO 40
      DO 30 KK=K,J1
      I=K+J1-KK
      Z2=Z1
      Z1=-PSI/BETA(I+1)
      Z(I+1)=Z2
      PSI=(ALFA(I)-E)*Z1+BETA(I+1)*Z2
   30 CONTINUE
C RESCALE THE LAST SET OF COMPONENTS
   40 W1=W1/Z1
      IF(K.GT.J1)GO TO 60
      DO 50 I=K,J1
      Z(I+1)=Z(I+1)*W1
   50 CONTINUE
   60 RETURN
      END