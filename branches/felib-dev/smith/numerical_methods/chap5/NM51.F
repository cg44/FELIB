      PROGRAM P51
C
C      PROGRAM 5.1  INTERPOLATION USING LAGRANGIAN POLYNOMIALS
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (MAX=10)
C
      REAL X(MAX),Y(MAX)
C
      WRITE (6,*) ('******* LAGRANGIAN INTERPOLATION ******')
      WRITE (6,100)
      READ (5,*) NP
      DO 1 I = 1,NP
          READ (5,*) X(I),Y(I)
    1 WRITE (6,102) X(I),Y(I)
      READ (5,*) XI
      YI = 0.
      DO 2 I = 1,NP
          FAC = 1.
          DO 3 J = 1,NP
              IF (J.EQ.I) GO TO 3
              FAC = FAC* (XI-X(J))/ (X(I)-X(J))
    3     CONTINUE
    2 YI = YI + FAC*Y(I)
      WRITE (6,101) NP - 1
      WRITE (6,102) XI,YI
  100 FORMAT (/,'DATA POINTS')
  101 FORMAT (/,'INTERPOLATED POINT USING POLYNOMIAL (ORDER',I2,')')
  102 FORMAT (2E12.4)
      STOP
      END
