      PROGRAM P32
C
C      PROGRAM 3.2 BISECTION METHOD FOR A SINGLE ROOT
C
      READ (5,*) XN,XN1,TOL,ITS
      WRITE(6,*) ('*************** BISECTION METHOD ***************')
      WRITE(6,*)
      WRITE(6,*) ('INITIAL VALUES')
      WRITE (6,100) XN,XN1
      ITERS = 0
      ICON = 0
      XOLD = XN
    1 IF (ICON.EQ.0 .AND. ITERS.LE.ITS) THEN
          XMID = 0.5* (XN+XN1)
          PROD = F(XN)*F(XMID)
          IF (PROD.LT..0) THEN
              XN1 = XMID
          ELSE
              XN = XMID
          END IF
          CALL CHECK(XMID,XOLD,TOL,ICON)
          ITERS = ITERS + 1
          GO TO 1
      END IF
      WRITE(6,*)
      WRITE(6,*) ('SOLUTION AND ITERATIONS TO CONVERGENCE')
      WRITE (6,200) XMID,ITERS
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
