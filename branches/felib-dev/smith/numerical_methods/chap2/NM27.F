      PROGRAM P27
C
C      PROGRAM 2.7 CHOLESKY FACTORISATION
C      WITH SKYLINE STORAGE
C      PRESCRIBED VARIABLES BY 'BIG SPRINGS'
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER (IA=50,IB=20,INO=5)
C
      REAL A(IA),B(IB),VAL(INO)
      INTEGER KDIAG(IB),NO(INO)
C
      READ (5,*) N
      READ (5,*) (KDIAG(I),I=1,N)
      IR = KDIAG(N)
      READ (5,*) (A(I),I=1,IR)
      READ (5,*) (B(I),I=1,N)
      READ (5,*) NFIX
      READ (5,*) (NO(I),VAL(I),I=1,NFIX)
      WRITE(6,*) ('** CHOLESKY FACTORISATION WITH SKYLINE STORAGE ***')
      WRITE(6,*) ('*********** SYMMETRICAL EQUATIONS ****************')
      WRITE(6,*) ('*********** PRESCRIBED VARIABLES *****************')
      WRITE(6,*)
      WRITE(6,*) ('COEFFICIENT VECTOR')
      CALL PRINTV(A,IR,6)
      WRITE(6,*)
      WRITE(6,*) ('DIAGONAL LOCATIONS')
      WRITE(6,100) (KDIAG(I),I=1,N)
      WRITE(6,*)
      WRITE(6,*) ('RIGHT HAND SIDE VECTOR')
      CALL PRINTV(B,N,6)
      WRITE(6,*)
      WRITE(6,*) ('FREEDOM NUMBER AND PRESCRIBED VALUE')
      DO 1 I=1,NFIX
    1 WRITE(6,200)NO(I),VAL(I)
      WRITE(6,*)
      DO 2 I=1,NFIX
      A(KDIAG(NO(I))) = A(KDIAG(NO(I))) + 1.E10
    2 B(NO(I)) = A(KDIAG(NO(I)))*VAL(I)
      CALL SKYFAC(A,N,KDIAG)
      CALL SKYSUB(A,B,N,KDIAG)
      WRITE(6,*) ('SOLUTION VECTOR')
      CALL PRINTV(B,N,6)
  100 FORMAT(10I6)
  200 FORMAT(I5,E12.4)
      STOP
      END
