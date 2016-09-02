      PROGRAM P52
C
C      PROGRAM 5.2  INTERPOLATION USING CUBIC SPLINES
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (MAX=10)
C
      REAL X(0:MAX),Y(0:MAX),U(MAX),F(MAX),D2(0:MAX)
      REAL DX(MAX),DY(MAX),K(MAX,MAX)
C
      WRITE (6,*) ('***** CUBIC SPLINE INTERPOLATION ******')
      WRITE (6,100)
      READ (5,*) NP
      N = NP - 1
      DO 1 I = 0,N
          READ (5,*) X(I),Y(I)
    1 WRITE (6,102) X(I),Y(I)
      READ (5,*) XI
      DO 2 I = 1,N
          DX(I) = X(I) - X(I-1)
    2 DY(I) = Y(I) - Y(I-1)
      CALL NULL(K,MAX,N-1,N-1)
      DO 3 I = 1,N - 1
          K(I,I) = 2.* (DX(I)+DX(I+1))
    3 F(I) = 6.* (DY(I+1)/DX(I+1)-DY(I)/DX(I))
      DO 4 I = 2,N - 1
          K(I-1,I) = DX(I)
    4 K(I,I-1) = DX(I)
      CALL SOLVE(K,MAX,U,F,N-1)
      DO 5 I = 1,N - 1
    5 D2(I) = U(I)
      D2(0) = 0.
      D2(N) = 0.
      I = 1
    6 IF (XI.LE.X(I)) GO TO 7
      I = I + 1
      GO TO 6
    7 YI = (D2(I-1)* (X(I)-XI)**3+D2(I)* (XI-X(I-1))**3)/ (6.*DX(I)) +
     +     (Y(I-1)/DX(I)-D2(I-1)*DX(I)/6.)* (X(I)-XI) +
     +     (Y(I)/DX(I)-D2(I)*DX(I)/6.)* (XI-X(I-1))
      WRITE (6,101)
      WRITE (6,102) XI,YI
  100 FORMAT (/,'DATA POINTS')
  101 FORMAT (/,'INTERPOLATED POINT USING CUBIC SPLINES')
  102 FORMAT (2E12.4)
      STOP
      END
