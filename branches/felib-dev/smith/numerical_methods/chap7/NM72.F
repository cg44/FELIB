      PROGRAM P72
C
C      PROGRAM 7.2  ONE-STEP METHODS FOR SYSTEMS OF EQUATIONS 
C
C      ITYPE= 1 (EULER'S METHOD)   ITYPE= 2 (MODIFIED EULER'S METHOD)
C      ITYPE= 3 (MID-POINT RULE)   ITYPE= 4 (4TH ORDER RUNGE-KUTTA)
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER(MEQ=10)
C
      REAL Y(MEQ),Y0(MEQ),K0(MEQ),K1(MEQ),K2(MEQ),K3(MEQ)
C
      READ (5,*) ITYPE,N,NSTEPS,H,X
      READ (5,*) (Y(I),I=1,N)
      WRITE (6,*) ('******* SYSTEMS OF EQUATIONS **********')
      WRITE (6,*)
      GO TO (1,2,3,4) ITYPE
    1 WRITE (6,*) ('************ EULERS METHODS ***********')
      WRITE (6,*)
      WRITE (6,200) N
      DO 10 J = 0,NSTEPS
          WRITE (6,100) X, (Y(I),I=1,N)
          CALL FUNC(K0,Y,X)
          DO 11 I = 1,N
   11         Y(I) = Y(I) + K0(I)*H
   10     X = X + H
      GO TO 50
    2 WRITE (6,*) ('******* MODIFIED EULERS METHOD ********')
      WRITE (6,*)
      WRITE (6,200) N
      DO 20 J = 0,NSTEPS
          WRITE (6,100) X, (Y(I),I=1,N)
          CALL FUNC(K0,Y,X)
          DO 21 I = 1,N
              Y0(I) = Y(I)
   21         Y(I) = Y0(I) + K0(I)*H
          X = X + H
          CALL FUNC(K1,Y,X)
          DO 22 I = 1,N
   22         Y(I) = Y0(I) + .5*H* (K0(I)+K1(I))
   20     CONTINUE
      GO TO 50
    3 WRITE (6,*) ('*********** MID-POINT RULE ************')
      WRITE (6,*)
      WRITE (6,200) N
      DO 30 J = 0,NSTEPS
          WRITE (6,100) X, (Y(I),I=1,N)
          CALL FUNC(K0,Y,X)
          DO 31 I = 1,N
              Y0(I) = Y(I)
   31         Y(I) = Y0(I) + .5*H*K0(I)
          X = X + .5*H
          CALL FUNC(K1,Y,X)
          DO 32 I = 1,N
   32         Y(I) = Y0(I) + H*K1(I)
   30     X = X + .5*H
      GO TO 50
    4 WRITE (6,*) ('******** 4TH ORDER RUNGE-KUTTA ********')
      WRITE (6,*)
      WRITE (6,200) N
      DO 40 J = 0,NSTEPS
          WRITE (6,100) X, (Y(I),I=1,N)
          CALL FUNC(K0,Y,X)
          DO 42 I = 1,N
              Y0(I) = Y(I)
   42         Y(I) = Y0(I) + K0(I)*.5*H
          X = X + H*.5
          CALL FUNC(K1,Y,X)
          DO 43 I = 1,N
   43         Y(I) = Y0(I) + K1(I)*.5*H
          CALL FUNC(K2,Y,X)
          DO 44 I = 1,N
   44         Y(I) = Y0(I) + K2(I)*H
          X = X + H*.5
          CALL FUNC(K3,Y,X)
          DO 45 I = 1,N
   45         Y(I) = Y0(I) + (K0(I)+2.* (K1(I)+K2(I))+K3(I))/6.*H
   40     CONTINUE
   50 CONTINUE
  100 FORMAT (10E13.5)
  200 FORMAT ('      X            Y(I) , I = 1,',I2)
      STOP
      END
C
      SUBROUTINE FUNC(F,Y,X)
C
C      THIS SUBROUTINE PROVIDES THE VALUES OF F(X,Y) FOR EACH
C      EQUATION AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      REAL F(*),Y(*)
C
      F(1) = 3.*x*y(2)+4.
      F(2) = x*y(1)-y(2)-exp(x)
      RETURN
      END
