C
      SUBROUTINE PRTVAL(VAL,IVAL,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRTVAL prints out the nodal values of the solution in a
C      standard format
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 1.1  29 Oct 1979 (CG)
C      Commented    14 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      VAL     vector of dimension IVAL containing solution
C              values,and prescribed boundary values
C      IVAL    dimension of VAL (.GE. TOTAL number of freedoms
C              in system)
C      NF      integer array of dimension (INF, JNF) containing
C              freedom numbers associated with each NODE
C      INF     first dimension of NF (.GE. TOTNOD)
C      JNF     second dimension of NF (.GE. DOFNOD)
C      DOFNOD  number of degrees of freedom at each NODE
C      TOTNOD  total number of nodes in mesh
C      NOUT    fortran unit number
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE PRTVAL(VAL,IVAL,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C***********************************************************************
C
      INTEGER DOFNOD,ERRMES,I,IERROR,INF,ITEST,IVAL,J,JNF,JTEST,K,N,NF,
     *        NODE,NOUT,TOTNOD
      CHARACTER*6 SRNAME
      DOUBLE PRECISION VAL,WORK
      DIMENSION VAL(IVAL),WORK(5),NF(INF,JNF),NODE(5)
      DATA SRNAME/'PRTVAL'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (INF.LT.TOTNOD .OR. JNF.LT.DOFNOD) IERROR = 2
         IF (DOFNOD.LE.0 .OR. TOTNOD.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      IF (DOFNOD.EQ.1) THEN
C
         WRITE (NOUT,FMT=9950)
         N = 0
         DO 1000 I = 1,TOTNOD
            N = N + 1
            NODE(N) = I
            K = NF(I,1)
            WORK(N) = 0.0D0
            IF (K.NE.0) THEN
               IF (K.LE.0) K = IABS(K)
               WORK(N) = VAL(K)
            END IF
            IF (N.EQ.4) THEN
               WRITE (NOUT,FMT=9930) (NODE(J),WORK(J),J=1,4)
               N = 0
            END IF
 1000    CONTINUE
C
         IF (N.NE.0) WRITE (NOUT,FMT=9930) (NODE(J),WORK(J),J=1,N)
      ELSE
         IF (DOFNOD.EQ.2) WRITE (NOUT,FMT=9990)
         IF (DOFNOD.EQ.3) WRITE (NOUT,FMT=9980)
         IF (DOFNOD.EQ.4) WRITE (NOUT,FMT=9970)
         IF (DOFNOD.GE.5) WRITE (NOUT,FMT=9960)
         N = 0
         DO 1020 I = 1,TOTNOD
            DO 1010 J = 1,DOFNOD
               N = N + 1
               NODE(N) = I
               K = NF(I,J)
C
C     Range checking on K
C
               IF (JTEST.NE.-1) THEN
                  IERROR = 0
                  IF (IVAL.LT.K) IERROR = 3
                  ITEST = ERRMES(JTEST,IERROR,SRNAME)
                  IF (ITEST.NE.0) RETURN
               END IF
C
               WORK(N) = 0.0D0
               IF (K.NE.0) THEN
                  IF (K.LE.0) K = IABS(K)
                  WORK(N) = VAL(K)
               END IF
               IF (((DOFNOD.NE.2).OR. (N.EQ.4)) .AND.
     *             (DOFNOD.EQ.2)) THEN
                  WRITE (NOUT,FMT=9920) NODE(1),WORK(1),WORK(2),
     *               NODE(3),WORK(3),WORK(4)
                  N = 0
               ELSE IF ((DOFNOD.NE.2) .AND. (N.EQ.DOFNOD)) THEN
                  WRITE (NOUT,FMT=9940) I, (WORK(K),K=1,DOFNOD)
                  N = 0
               END IF
 1010       CONTINUE
 1020    CONTINUE
C
         IF (N.NE.0) THEN
            IF (DOFNOD.EQ.2) WRITE (NOUT,FMT=9920) NODE(1),WORK(1),
     *          WORK(2)
            IF (DOFNOD.NE.2) WRITE (NOUT,FMT=9940) I, (WORK(K),K=1,N)
         END IF
      END IF
C
 9990 FORMAT (/' ',2 ('NODE',6X,2 ('VALUE',9X)),/' ')
 9980 FORMAT (/' NODE',6X,3 ('VALUE',9X),/' ')
 9970 FORMAT (/' NODE',6X,4 ('VALUE',9X),/' ')
 9960 FORMAT (/' NODE',6X,5 ('VALUE',9X),/' ')
 9950 FORMAT (/' ',4 ('NODE',5X,'VALUE',5X),/' ')
 9940 FORMAT (' ',I4,1X,5 (D12.5,2X),/' ',5X,5 (D12.5,2X),/' ',5X,
     *       5 (D12.5,2X))
 9930 FORMAT (' ',4 (I4,1X,D12.5,2X))
 9920 FORMAT (' ',2 (I4,1X,2 (D12.5,2X),5X))
      END
