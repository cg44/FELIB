      PROGRAM P63
C
C      PROGRAM 6.3  GAUSS-LAGUERRE RULE
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (ISAMP=5)
C
      REAL SAMP(ISAMP,2)
C
      READ (5,*) NGL
      CALL GAULAG(SAMP,ISAMP,NGL)
      AREA = 0.
      DO 10 I = 1,NGL
          X = SAMP(I,1)
          W = SAMP(I,2)
   10 AREA = AREA + W*F(X)
      WRITE (6,*) ('****** GAUSS-LAGUERRE RULE *****')
      WRITE (6,101) NGL
      WRITE (6,103) AREA
  101 FORMAT (/,'  NUMBER OF GAUSS POINTS                ',I7)
  103 FORMAT ('  COMPUTED RESULT                       ',F12.4)
      STOP
      END
C
      FUNCTION F(X)
C
C      THIS FUNCTION PROVIDES THE VALUE OF F(X)
C      AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      F = SIN(X)
      RETURN
      END
