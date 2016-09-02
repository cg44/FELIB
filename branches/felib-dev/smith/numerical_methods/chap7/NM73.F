      PROGRAM P73
C
C      PROGRAM 7.3  PREDICTOR CORRECTOR METHOD USING
C                   4TH ORDER ADAMS/BASHFORTH/MOULTON
C
C      FOUR INITIAL VALUES OF X AND Y PROVIDED (NIC=4)
C
C      ONE INITIAL VALUE PROVIDED AND THE REMAINING THREE
C      GENERATED USING 4TH ORDER RUNGE-KUTTA (NIC=1)
C
      REAL X(5),Y(5)
C
      READ (5,*) NIC,NSTEPS,H,TOL
      WRITE (6,*) ('********* PREDICTOR CORRECTOR METHOD *********')
      WRITE (6,*)
      WRITE (6,*) ('********** ADAMS/BASHFORTH/MOULTON ***********')
      WRITE (6,*)
      IF (NIC.EQ.4) THEN
          WRITE (6,*) ('********** 4 INITIAL VALUES PROVIDED *********')
          WRITE (6,*)
          WRITE (6,*) (
     +      '     X            Y           ERR      ITERATIONS')
          DO 10 I = 1,4
              READ (5,*) X(I),Y(I)
   10     WRITE (6,100) X(I),Y(I)
      ELSE
          WRITE (6,*) ('********** 1 INITIAL VALUE PROVIDED **********')
          WRITE (6,*) ('*** REMAINING 3 FROM 4TH ORDER RUNGE-KUTTA ***')
          WRITE (6,*)
          WRITE (6,*) (
     +      '     X            Y           ERR      ITERATIONS')
          READ (5,*) X(1),Y(1)
          WRITE (6,100) X(1),Y(1)
          DO 20 I = 1,3
              RK0 = H*F(X(I),Y(I))
              RK1 = H*F(X(I)+0.5*H,Y(I)+0.5*RK0)
              RK2 = H*F(X(I)+0.5*H,Y(I)+0.5*RK1)
              RK3 = H*F(X(I)+H,Y(I)+RK2)
              X(I+1) = X(I) + H
              Y(I+1) = Y(I) + (RK0+2.*RK1+2.*RK2+RK3)/6.
   20     WRITE (6,100) X(I+1),Y(I+1)
      END IF
      DO 30 I = 1,NSTEPS
          Y(5) = Y(4) + H/24.* (55.*F(X(4),Y(4))-59.*F(X(3),Y(3))+
     +           37.*F(X(2),Y(2))-9.*F(X(1),Y(1)))
          YHOLD = Y(5)
          X(5) = X(4) + H
          ITERS = 0
   40     ITERS = ITERS + 1
          Y(5) = Y(4) + H/24.* (9.*F(X(5),Y(5))+19.*F(X(4),Y(4))-
     +           5.*F(X(3),Y(3))+F(X(2),Y(2)))
          IF (ITERS.EQ.1) ERR = -1./14.* (Y(5)-YHOLD)
          IF (ABS((Y(5)-YHOLD)/YHOLD).LE.TOL) GO TO 50
          YHOLD = Y(5)
          GO TO 40
   50     WRITE (6,100) X(5),Y(5),ERR,ITERS
          DO 60 J = 1,4
              Y(J) = Y(J+1)
   60     X(J) = X(J+1)
   30 CONTINUE
  100 FORMAT (3E13.5,I8)
      STOP
      END
C
      FUNCTION F(X,Y)
C
C      THIS FUNCTION PROVIDES THE VALUE OF F(X,Y))
C      AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
      F = X*X + Y*Y
      RETURN
      END
