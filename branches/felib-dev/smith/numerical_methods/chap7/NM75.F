      PROGRAM P75
C
C      PROGRAM 7.5  SHOOTING METHODS FOR NONLINEAR 2ND ORDER
C                   BOUNDARY VALUE PROBLEMS
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (MAX=25)
C
      REAL Y(2),Y0(MAX),Y1(MAX),YS(MAX)
C
      READ (5,*) NSTEPS,H,XA,YA,XB,YB,A0,A1,TOL,ITS
      WRITE (6,*) ('******** SHOOTING METHODS FOR *********')
      WRITE (6,*) ('******** NONLINEAR BVP"S USING ********')
      WRITE (6,*) ('******** LINEAR INTERPOLATION *********')
      WRITE (6,*)
      X = XA
      Y(1) = YA
      Y(2) = A0
      Y0(1) = YA
      DO 10 I = 2,NSTEPS + 1
          CALL RUNG4(Y,2,X,H)
   10 Y0(I) = Y(1)
      B0 = Y0(NSTEPS+1)
      X = XA
      Y(1) = YA
      Y(2) = A1
      Y1(1) = YA
      DO 20 I = 2,NSTEPS + 1
          CALL RUNG4(Y,2,X,H)
   20 Y1(I) = Y(1)
      B1 = Y1(NSTEPS+1)
      IF ((B0-YB)* (B1-YB).GT.0.) THEN
          WRITE (6,*) ('TRY NEW GRADIENTS')
          STOP
      END IF
      ITERS = 0
   30 ITERS = ITERS + 1
      AS = A0 + (YB-B0)* (A1-A0)/ (B1-B0)
      X = XA
      Y(1) = YA
      Y(2) = AS
      YS(1) = YA
      DO 40 I = 2,NSTEPS + 1
          CALL RUNG4(Y,2,X,H)
   40 YS(I) = Y(1)
      BS = YS(NSTEPS+1)
      IF ((BS-YB)* (B0-YB).GT.0.) THEN
          A0 = AS
          DO 50 I = 2,NSTEPS + 1
   50     Y0(I) = YS(I)
          B0 = BS
      ELSE
          A1 = AS
          DO 60 I = 2,NSTEPS + 1
   60     Y1(I) = YS(I)
          B1 = BS
      END IF
      IF (ABS((BS-YB)/YB).GT.TOL .AND. ITERS.LT.ITS) GO TO 30
      WRITE (6,*) ('     X            Y')
      X = XA
      WRITE (6,100) X,YS(1)
      DO 70 I = 2,NSTEPS + 1
          X = X + H
   70 WRITE (6,100) X,YS(I)
      WRITE (6,*)
      WRITE (6,*) 'AFTER',ITERS,'ITERATIONS'
  100 FORMAT (10E13.5)
      STOP
      END
C
      SUBROUTINE RUNG4(Y,N,X,H)
C
C      ADVANCES A SYSTEM OF 1ST ORDER
C      EQUATIONS USING THE 4TH ORDER RUNGE-KUTTA METHOD
C
      REAL Y(10),K0(10),K1(10),K2(10),K3(10),Y0(10)
C
      CALL FUNC(K0,Y,X)
      DO 2 I = 1,N
          Y0(I) = Y(I)
    2 Y(I) = Y0(I) + K0(I)*.5*H
      X = X + H*.5
      CALL FUNC(K1,Y,X)
      DO 3 I = 1,N
    3 Y(I) = Y0(I) + K1(I)*.5*H
      CALL FUNC(K2,Y,X)
      DO 4 I = 1,N
    4 Y(I) = Y0(I) + K2(I)*H
      X = X + H*.5
      CALL FUNC(K3,Y,X)
      DO 5 I = 1,N
    5 Y(I) = Y0(I) + (K0(I)+2.* (K1(I)+K2(I))+K3(I))/6.*H
      RETURN
      END
C
      SUBROUTINE FUNC(F,Y,X)
C
C      THIS SUBROUTINE PROVIDES THE VALUES OF F(X,Y) FOR EACH
C      EQUATION AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      REAL F(2),Y(2)
 
      F(1) = Y(2)
      F(2) = -2*Y(2)**2/Y(1)
      RETURN
      END
