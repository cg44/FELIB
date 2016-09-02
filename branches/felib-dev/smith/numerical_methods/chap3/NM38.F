      PROGRAM P38
C
C      PROGRAM 3.8 MODIFIED  NEWTON-RAPHSON FOR MANY EQUATIONS
C
C      ALETR NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20)
C
      REAL X(IN),X1(IN),DELX(IN),F(IN),DF(IN,IN)
C
      READ (5,*) N,(X(I),I=1,N)
      READ (5,*) TOL,ITS
      WRITE(6,*) ('************ MODIFIED NEWTON-RAPHSON *************')
      WRITE(6,*) ('*********** FOR A SYSTEM OF EQUATIONS ************')
      WRITE(6,*)
      WRITE(6,*) ('INITIAL VALUES')
      CALL PRINTV(X,N,6)
      WRITE(6,*)
      ITERS = 0
      ICON = 0
      CALL FANDDF(X,F,DF,IN,N)
      CALL MATINV(DF,IN,N)
    1 CALL FUN(X,F)
      CALL MVMULT(DF,IN,F,N,N,DELX)
      CALL VECSUB(X,DELX,X1,N)
      CALL CHECON(X1,X,N,TOL,ICON)
      ITERS = ITERS + 1
      IF (ICON.EQ.0 .AND. ITERS.LT.ITS) GO TO 1
      WRITE(6,*) ('SOLUTIONS AND ITERATIONS TO CONVERGENCE')
      CALL PRINTV(X1,N,6)
      WRITE (6,100) ITERS
  100 FORMAT (I10)
      STOP
      END
C
      SUBROUTINE FANDDF(X,F,DF,IN,N)
C
C      THIS SUBROUTINE PROVIDES THE VALUES OF THE FUNCTIONS
C      OF X AND THEIR DERIVATIVES, AND WILL VARY FROM ONE
C      PROBLEM TO THE NEXT
C
C
      REAL X(*),F(*),DF(IN,*)
C
      F(1) = 2.0*X(1)*X(1) + X(2)*X(2) - 4.32
      F(2) = X(1)*X(1) - X(2)*X(2)
      DF(1,1) = 4.0*X(1)
      DF(2,1) = 2.0*X(1)
      DF(1,2) = 2.0*X(2)
      DF(2,2) = -2.0*X(2)
      RETURN
      END
C
      SUBROUTINE FUN(X,F)
C
C      THIS SUBROUTINE PROVIDES THE VALUES OF THE FUNCTIONS
C      OF X AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
C
      REAL X(*),F(*)
C
      F(1) = 2.0*X(1)*X(1) + X(2)*X(2) - 4.32
      F(2) = X(1)*X(1) - X(2)*X(2)
      RETURN
      END
