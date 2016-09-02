      PROGRAM P62
C
C      PROGRAM 6.2  REPEATED GAUSS-LEGENDRE RULE
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (ISAMP=7)
C
      REAL SAMP(ISAMP,2)
C
      READ (5,*) A,B,NGP,NR
      CALL GAULEG(SAMP,ISAMP,NGP)
      WR = (B-A)/NR
      HR = 0.5*WR
      AREA = 0.
      DO 10 I = 1,NR
          CR = A + (I-1)*WR + HR
          DO 10 J = 1,NGP
              X = SAMP(J,1)
              W = SAMP(J,2)
              XS = CR + X*HR
   10 AREA = AREA + W*HR*F(XS)
      WRITE (6,*) ('****** REPEATED GAUSS-LEGENDRE RULE *****')
      WRITE (6,100) A,B
      WRITE (6,101) NGP
      WRITE (6,102) NR
      WRITE (6,103) AREA
  100 FORMAT (/,'  LIMITS OF INTEGRATION                 ',2F12.4)
  101 FORMAT ('  NUMBER OF GAUSS POINTS PER STRIP      ',I7)
  102 FORMAT ('  NUMBER OF REPETITIONS                 ',I7)
  103 FORMAT ('  COMPUTED RESULT                       ',F12.4)
      STOP
      END
C
      FUNCTION F(X)
C
C      THIS FUNCTION PROVIDES THE VALUE OF F(X)
C      AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      F = COS(X)*X
      RETURN
      END
