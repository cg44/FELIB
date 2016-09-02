      PROGRAM P24
C
C      PROGRAM 2.4 CHOLESKY FACTORISATION
C                  WITH BANDED STORAGE
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IN=20,IIW=10)
C
      REAL LB(IN,IIW),B(IN)
C
      READ (5,*) N,IW
      IWP1 = IW + 1
      READ (5,*) ((LB(I,J),J=1,IWP1),I=1,N)
      READ (5,*) (B(I),I=1,N)
      WRITE(6,*) ('*** CHOLESKY FACTORISATION WITH BANDED STORAGE ****')
      WRITE(6,*) ('************ SYMMETRICAL EQUATIONS ****************')
      WRITE(6,*)
      WRITE(6,*) ('BANDED COEFFICIENT MATRIX')
      DO 1 I=1,N
    1 WRITE(6,100) (LB(I,J),J=1,IWP1)
      WRITE(6,*)
      WRITE(6,*) ('RIGHT HAND SIDE VECTOR')
      CALL PRINTV(B,N,6)
      WRITE(6,*)
      CALL CHOFAC(LB,IN,N,IW)
      WRITE(6,*) ('LT IN BAND FORM')
      CALL PRINTA(LB,IN,N,IWP1,6)
      WRITE(6,*)
      CALL CHOSUB(LB,IN,B,N,IW)
      WRITE(6,*) ('SOLUTION VECTOR')
      CALL PRINTV(B,N,6)
  100 FORMAT(5E12.4)
      STOP
      END
