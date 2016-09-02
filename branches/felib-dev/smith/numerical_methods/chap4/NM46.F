      PROGRAM P46
C
C      PROGRAM 4.6 JACOBI TRANSFORMATION FOR EIGENVALUES
C                  OF SYMMETRIC MATRICES
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE)
C
      PARAMETER (IN=10)
C
      REAL A(IN,IN),A1(IN,IN),P(IN,IN),PT(IN,IN),ENEW(IN),EOLD(IN)
C
      PI = 4.*ATAN(1.)
      READ(5,*) N,((A(I,J),J=I,N),I=1,N)
      DO 2 I = 1,N
          DO 2 J = I,N
    2 A(J,I) = A(I,J)
      READ (5,*) TOL,ITS
      WRITE(6,*) ('************ JACOBI ITERATION *******************')
      WRITE(6,*) ('******** FOR SYMMETRIC MATRICES ***************')
      WRITE(6,*)
      WRITE(6,*) ('MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      CALL NULVEC(EOLD,N)
      ITERS = 0
   99 ITERS = ITERS + 1
      BIG = 0.
      DO 3 I = 1,N
          DO 3 J = I + 1,N
              IF (ABS(A(I,J)).GT.BIG) THEN
                  BIG = ABS(A(I,J))
                  HOLD = A(I,J)
                  IROW = I
                  ICOL = J
              END IF
    3 CONTINUE
      DEN = A(IROW,IROW) - A(ICOL,ICOL)
      IF (DEN.EQ.0.) THEN
          ALPHA = PI*.25
          IF (HOLD.LT.0.) ALPHA = -ALPHA
      ELSE
          ALPHA = .5*ATAN(2.*HOLD/DEN)
      END IF
      CT = COS(ALPHA)
      ST = SIN(ALPHA)
      CALL NULL(P,IN,N,N)
      DO 4 I = 1,N
    4 P(I,I) = 1.
      P(IROW,IROW) = CT
      P(ICOL,ICOL) = CT
      P(IROW,ICOL) = -ST
      P(ICOL,IROW) = ST
      ALPHA = ALPHA*180./PI
      CALL MATRAN(PT,IN,P,IN,N,N)
      CALL MATMUL(PT,IN,A,IN,A1,IN,N,N,N)
      CALL MATMUL(A1,IN,P,IN,A,IN,N,N,N)
      DO 6 I = 1,N
    6 ENEW(I) = A(I,I)
      CALL CHECON(ENEW,EOLD,N,TOL,ICON)
      IF (ITERS.LT.ITS .AND. ICON.EQ.0) GO TO 99
      WRITE(6,*) ('FINAL TRANSFORMATION MATRIX P')
      CALL PRINTA(P,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('FINAL TRANSFORMED MATRIX A')
      CALL PRINTA(A,IN,N,N,6)
      WRITE(6,*)
      WRITE(6,*) ('ITERATIONS TO CONVERGENCE')
      WRITE(6,1000) ITERS
 1000 FORMAT (I5)
      STOP
      END
