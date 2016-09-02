      PROGRAM P74
C
C      PROGRAM 7.4  SHOOTING METHODS FOR LINEAR 2ND ORDER
C                   BOUNDARY VALUE PROBLEMS
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (MAX=25)
C
      REAL Y(2),Y0(MAX),Y1(MAX),YS(MAX)
C
      READ (5,*) NSTEPS,H,XA,YA,XB,YB,A0,A1
      WRITE (6,*) ('******** SHOOTING METHOD FOR **********')
      WRITE (6,*) ('******** LINEAR BVP"S USING ***********')
      WRITE (6,*) ('******* 2 TRIAL INITIAL GRADIENTS *****')
      WRITE (6,*)
      X = XA
      Y(1) = YA
      Y(2) = AO
      Y0(1) = YA
      DO 10 I = 2,NSTEPS + 1
          CALL RUNG4(Y,X,H)
   10 Y0(I) = Y(1)
      B0 = Y0(NSTEPS+1)
      X = XA
      Y(1) = YA
      Y(2) = A1
      Y1(1) = YA
      DO 20 I = 2,NSTEPS + 1
          CALL RUNG4(Y,X,H)
   20 Y1(I) = Y(1)
      B1 = Y1(NSTEPS+1)
      WRITE (6,*) ('     X            Y')
      X = XA
      DO 30 I = 1,NSTEPS + 1
          YS(I) = ((YB-B1)*Y0(I)+ (B0-YB)*Y1(I))/ (B0-B1)
          WRITE (6,100) X,YS(I)
   30 X = X + H
  100 FORMAT (10E13.5)
      STOP
      END
C
      SUBROUTINE RUNG4(Y,X,H)
C
C      ADVANCES TWO 1ST ORDER EQUATIONS
C      USING THE 4TH ORDER RUNGE-KUTTA METHOD
C
      REAL Y(2),K0(2),K1(2),K2(2),K3(2),Y0(2)
C
      CALL FUNC(K0,Y,X)
      DO 2 I = 1,2
          Y0(I) = Y(I)
    2 Y(I) = Y0(I) + K0(I)*.5*H
      X = X + H*.5
      CALL FUNC(K1,Y,X)
      DO 3 I = 1,2
    3 Y(I) = Y0(I) + K1(I)*.5*H
      CALL FUNC(K2,Y,X)
      DO 4 I = 1,2
    4 Y(I) = Y0(I) + K2(I)*H
      X = X + H*.5
      CALL FUNC(K3,Y,X)
      DO 5 I = 1,2
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
C
      F(1) = Y(2)
      F(2) = 3.*X + 4.*Y(1)
      RETURN
      END
