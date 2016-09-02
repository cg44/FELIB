      PROGRAM P53
C
C      PROGRAM 5.3  CURVE FITTING USING 'LEAST SQUARES'
 
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (MAX=10,IC=5)
C
      REAL X(MAX),Y(MAX),F(IC),C(IC),LS(IC,IC),R(IC)
C
      WRITE (6,*) ('******* LEAST SQUARES CURVE FITTING *******')
      WRITE (6,100)
      READ (5,*) NP,NC
      DO 1 I = 1,NP
          READ (5,*) X(I),Y(I)
    1 WRITE (6,101) X(I),Y(I)
      READ (5,*) XI
      CALL NULL(LS,IC,NC,NC)
      CALL NULVEC(R,NC)
      DO 2 I = 1,NP
          CALL FUNC(F,X(I))
          DO 3 J = 1,NC
              R(J) = R(J) + F(J)*Y(I)
              DO 3 K = 1,NC
    3     LS(J,K) = LS(J,K) + F(J)*F(K)
    2 CONTINUE
      CALL SOLVE(LS,IC,C,R,NC)
      CALL FUNC(F,XI)
      YI = 0.
      DO 4 I = 1,NC
    4 YI = YI + C(I)*F(I)
      WRITE (6,*)
      WRITE (6,*) ('OPTIMISED FUNCTION COEFFICIENTS')
      WRITE (6,101) (C(I),I=1,NC)
      WRITE (6,*)
      WRITE (6,*) ('INTERPOLATED POINT')
      WRITE (6,101) XI,YI
  100 FORMAT (/,'DATA POINTS')
  101 FORMAT (5E12.4)
      STOP
      END
C
      SUBROUTINE FUNC(F,X)
C
C      THIS SUBROUTINE PROVIDES THE VALUES OF THE REQUIRED
C      FUNCTIONS AT A PARTICULAR VALUE OF X
C
      REAL F(*)
C
      F(1) = 1.
      F(2) = ALOG(X)
      RETURN
      END
