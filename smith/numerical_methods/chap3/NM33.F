      PROGRAM P33
C
C      PROGRAM 3.3 FALSE POSITION METHOD FOR A SINGLE ROOT
C
      READ (5,*) XN,XN1,TOL,ITS
      WRITE(6,*) ('********** FALSE POSITION METHOD ***************')
      WRITE(6,*)
      WRITE(6,*) ('INITIAL VALUES')
      WRITE (6,100) XN,XN1
      ITERS = 0
      ICON = 0
      XOLD = XN
    1 IF (ICON.EQ.0 .AND. ITERS.LE.ITS) THEN
          XNEW = XN - F(XN)* (XN1-XN)/ (F(XN1)-F(XN))
          PROD = F(XN)*F(XNEW)
          IF (PROD.LT..0) THEN
              XN1 = XNEW
          ELSE
              XN = XNEW
          END IF
          CALL CHECK(XNEW,XOLD,TOL,ICON)
          ITERS = ITERS + 1
          XOLD = XNEW
          GO TO 1
      END IF
      WRITE(6,*)
      WRITE(6,*) ('SOLUTION AND ITERATIONS TO CONVERGENCE')
      WRITE (6,200) XNEW,ITERS
  100 FORMAT (2E12.4)
  200 FORMAT (E12.4,I10)
      STOP
      END
C
      FUNCTION F(X)
C
C      THIS FUNCTION PROVIDES THE VALUE OF F(X)
C      AND WILL VARY FROM ONE PROBLEM TO THE NEXT
C
C
      F = X**3 - X - 1.0
      RETURN
      END
