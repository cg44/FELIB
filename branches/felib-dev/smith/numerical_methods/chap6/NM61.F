      PROGRAM P61
C
C      PROGRAM 6.1  REPEATED NEWTON-COTES RULE
C
      REAL WEIGHT(5)
C
      READ (5,*) A,B,NNCP,NR
      CALL NEWCOT(WEIGHT,NNCP)
      WR = (B-A)/NR
      H = WR
      IF (NNCP.GE.2) H = WR/ (NNCP-1)
      AREA = 0.
      DO 10 I = 1,NR
          BR = A + (I-1)*WR
          DO 10 J = 1,NNCP
              XS = BR + (J-1)*H
              W = WEIGHT(J)
   10 AREA = AREA + W*H*F(XS)
      WRITE (6,*) ('****** REPEATED NEWTON-COTES RULE *****')
      WRITE (6,100) A,B
      WRITE (6,101) NNCP
      WRITE (6,102) NR
      WRITE (6,103) AREA
  100 FORMAT (/,'  LIMITS OF INTEGRATION                    ',2F12.4)
  101 FORMAT ('  NUMBER OF SAMPLING POINTS PER STRIP      ',I7)
  102 FORMAT ('  NUMBER OF REPETITIONS                    ',I7)
  103 FORMAT ('  COMPUTED RESULT                          ',1F12.4)
      STOP
      END
C
      FUNCTION F(X)
C
C      THIS FUNCTION PROVIDES THE VALUE OF F(X)
C      AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      F = SIN(X)*SIN(X) + X
      RETURN
      END
