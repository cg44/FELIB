      SUBROUTINE PRTTOP(TOTELS,ELTOP,IELTOP,JELTOP,ITEST)
      INTEGER TOTELS,ELTOP,IELTOP,JELTOP,ITEST
      INTEGER NOUT,I,L,K,J
      INTEGER ERRMES,IERROR
      DOUBLE PRECISION SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP)
      DATA NOUT /6/,SRNAME /8H PRTTOP /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IELTOP.LT.TOTELS) IERROR=2
      IF(TOTELS.LE.0) IERROR=1
      ITEST=ERRMES(IERROR,ITEST,SRNAME)
      IF(ITEST.NE.0) RETURN
999   WRITE(NOUT,9999)
      WRITE(NOUT,9998) TOTELS
      WRITE(NOUT,9997)
      DO 100 I=1,TOTELS
      L=ELTOP(I,2)
      K=L+2
      IF(ITEST.EQ.-1) GO TO 998
      IERROR=0
      IF(JELTOP.LT.K) IERROR=3
      ITEST=ERRMES(IERROR,ITEST,SRNAME)
      IF(ITEST.NE.0) RETURN
998   WRITE(NOUT,9996) I, (ELTOP(I,J),J=1,K)
100   CONTINUE
      RETURN
9999  FORMAT(1H ,////,27H **** ELEMENT TOPOLOGY ****//)
9998  FORMAT(1H ,20HNUMBER OF ELEMENTS =,I3)
9997   FORMAT(1H0,2X,4HELEM,4X,5HELTYP,4X,5HNODEL,4X,5HNODES)
9996  FORMAT(1H0,2X,I3,5X,I3,7X,I3,4X,15(I3,2X))
      END