      PROGRAM P71
C
C      PROGRAM 7.1  ONE-STEP METHODS (SINGLE EQUATION)
C
C      ITYPE= 1 (EULER'S METHOD)   ITYPE= 2 (MODIFIED EULER'S METHOD)
C      ITYPE= 3 (MID-POINT RULE)   ITYPE= 4 (4TH ORDER RUNGE-KUTTA)
C
      READ (5,*) ITYPE,NSTEPS,H
      READ (5,*) X,Y
      WRITE (6,*) ('*********** ONE-STEP METHODS **********')
      WRITE (6,*)
      GO TO (1,2,3,4) ITYPE
    1 WRITE (6,*) ('************ EULERS METHODS ***********')
      WRITE (6,*)
      WRITE (6,*) ('     X            Y')
      DO 10 I = 0,NSTEPS
          WRITE (6,100) X,Y
          RK0 = H*F(X,Y)
          X = X + H
   10 Y = Y + RK0
      GO TO 50
    2 WRITE (6,*) ('******* MODIFIED EULERS METHOD ********')
      WRITE (6,*)
      WRITE (6,*) ('     X            Y')
      DO 20 I = 0,NSTEPS
          WRITE (6,100) X,Y
          RK0 = H*F(X,Y)
          RK1 = H*F(X+H,Y+RK0)
          X = X + H
   20 Y = Y + 0.5* (RK0+RK1)
      GO TO 50
    3 WRITE (6,*) ('*********** MID-POINT RULE ************')
      WRITE (6,*)
      WRITE (6,*) ('     X            Y')
      DO 30 I = 0,NSTEPS
          WRITE (6,100) X,Y
          RK0 = H*F(X,Y)
          RK1 = H*F(X+0.5*H,Y+0.5*RK0)
          X = X + H
   30 Y = Y + RK1
      GO TO 50
    4 WRITE (6,*) ('******** 4TH ORDER RUNGE-KUTTA ********')
      WRITE (6,*)
      WRITE (6,*) ('     X            Y')
      DO 40 I = 0,NSTEPS
          WRITE (6,100) X,Y
          RK0 = H*F(X,Y)
          RK1 = H*F(X+0.5*H,Y+0.5*RK0)
          RK2 = H*F(X+0.5*H,Y+0.5*RK1)
          RK3 = H*F(X+H,Y+RK2)
          X = X + H
   40 Y = Y + (RK0+2.*RK1+2.*RK2+RK3)/6.
   50 CONTINUE
  100 FORMAT (2E13.5)
      STOP
      END
C
      FUNCTION F(X,Y)
C
C      THIS FUNCTION PROVIDES THE VALUE OF F(X,Y)
C      AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      F = (Y+X)**2
      RETURN
      END
