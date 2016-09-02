      PROGRAM P64
C
C      PROGRAM 6.4 AREA INTEGRAL OVER A QUADRILATERAL USING
C      GAUSS-LEGENDRE QUADRATURE
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (ISAMP=7)
C
      REAL SAMP(ISAMP,2),COORD(4,2),FUN(4),JAC(2,2),DER(2,4)
C
      DO 10 I = 1,4
   10 READ (5,*) COORD(I,1),COORD(I,2)
      READ (5,*) NGP
      CALL GAULEG(SAMP,ISAMP,NGP)
      VOL = 0.
      DO 20 I = 1,NGP
          DO 20 J = 1,NGP
              CALL BILIN(DER,FUN,SAMP,ISAMP,I,J)
              CALL MATMUL(DER,2,COORD,4,JAC,2,2,4,2)
              DET = JAC(1,1)*JAC(2,2) - JAC(1,2)*JAC(2,1)
              X = 0.
              Y = 0.
              DO 30 K = 1,4
                  X = X + FUN(K)*COORD(K,1)
   30         Y = Y + FUN(K)*COORD(K,2)
   20 VOL = VOL + DET*SAMP(I,2)*SAMP(J,2)*F(X,Y)
      WRITE (6,*) ('****** AREA INTEGRAL OVER A QUADRILATERAL ******')
      WRITE (6,*) ('************ GAUSS-LEGENDRE IN 2-D *************')
      WRITE (6,*)
      DO 40 I = 1,4
   40 WRITE (6,100) COORD(I,1),COORD(I,2)
      WRITE (6,101) NGP
      WRITE (6,103) VOL
  100 FORMAT ('  X,Y COORDINATES                       ',2F12.4)
  101 FORMAT (/,'  NO. OF GAUSS POINTS IN EACH DIRECTION ',I7)
  103 FORMAT (/,'  COMPUTED RESULT                       ',F12.4)
      STOP
      END
C
      FUNCTION F(X,Y)
C
C      THIS FUNCTION PROVIDES THE VALUE OF F(X,Y)
C      AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      F = X**2*Y**2
      RETURN
      END
