C
      SUBROUTINE CPRTVL(VAL,IVAL,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CPRTVL prints out the nodal values of the solution in a
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
C      Release 2.0  20 Jan 1984 (CRIE)
C      Commented     1 Nov 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      VAL     complex vector of dimension (2*IVAL) containing
C              solution values,and prescribed boundary values
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
C     SUBROUTINE CPRTVL(VAL,IVAL,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C***********************************************************************
C
      INTEGER DOFNOD,ERRMES,I,IERROR,INF,ITEST,IVAL,J,JNF,K,L,NODE,NOUT,
     *        TOTNOD,N,NF
      CHARACTER*6 SRNAME
      DOUBLE PRECISION VAL,WORK
      DIMENSION VAL(2,IVAL),WORK(2,10),NF(INF,JNF),NODE(2)
      DATA SRNAME/'CPRTVL'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (INF.LT.TOTNOD .OR. JNF.LT.DOFNOD) IERROR = 2
         IF (DOFNOD.LE.0 .OR. TOTNOD.LE.0 .OR. NOUT.LE.0) IERROR = 1
         IF (DOFNOD.GT.10) IERROR = 1
         ITEST = ERRMES(IERROR,ITEST,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      IF (DOFNOD.EQ.1) THEN
C
C     For DOFNOD=1
C
         WRITE (NOUT,FMT=9920)
         N = 0
         DO 1000 I = 1,TOTNOD
            N = N + 1
            NODE(N) = I
            K = NF(I,1)
            WORK(1,N) = 0.0D0
            WORK(2,N) = 0.0D0
            IF (K.NE.0) THEN
               WORK(1,N) = VAL(1,K)
               WORK(2,N) = VAL(2,K)
            END IF
            IF (N.EQ.2) THEN
               WRITE (NOUT,FMT=9910) (NODE(J),WORK(1,J),WORK(2,J),J=1,
     *            2)
               N = 0
            END IF
 1000    CONTINUE
         IF (N.NE.0) WRITE (NOUT,FMT=9910) NODE(1),WORK(1,1),WORK(2,1)
      ELSE
         WRITE (NOUT,FMT=9990)
         DO 1080 I = 1,TOTNOD
            DO 1010 J = 1,DOFNOD
               K = NF(I,J)
               WORK(1,J) = 0.0D0
               WORK(2,J) = 0.0D0
               IF (K.NE.0) THEN
                  WORK(1,J) = VAL(1,K)
                  WORK(2,J) = VAL(2,K)
               END IF
 1010       CONTINUE
C
C     SELECT number of degrees of freedom
C
            GO TO (1080,1020,1060,1030,1060,1040,1060,1050,
     *             1060,1070) DOFNOD
            GO TO 1060
 1020       CONTINUE
            WRITE (NOUT,FMT=9970) I, (WORK(1,L),WORK(2,L),L=1,DOFNOD)
            GO TO 1080
 1030       CONTINUE
            WRITE (NOUT,FMT=9960) I, (WORK(1,L),WORK(2,L),L=1,DOFNOD)
            GO TO 1080
 1040       CONTINUE
            WRITE (NOUT,FMT=9950) I, (WORK(1,L),WORK(2,L),L=1,DOFNOD)
            GO TO 1080
 1050       CONTINUE
            WRITE (NOUT,FMT=9940) I, (WORK(1,L),WORK(2,L),L=1,DOFNOD)
            GO TO 1080
 1060       CONTINUE
            WRITE (NOUT,FMT=9980) I, (WORK(1,L),WORK(2,L),L=1,DOFNOD)
            GO TO 1080
 1070       CONTINUE
            WRITE (NOUT,FMT=9930) I, (WORK(1,L),WORK(2,L),L=1,DOFNOD)
 1080    CONTINUE
      END IF
C
 9990 FORMAT (/' ','NODE',6X,2 ('REAL',7X,'IMAGINARY',16X),/' ')
 9980 FORMAT (' ',I4,2X,2 (D12.5,2X),8X,2 (D12.5,2X),
     *       8 (/' ',6X,2 (D12.5,2X),8X,2 (D12.5,2X)))
 9970 FORMAT (' ',I4,2X,2 (D12.5,2X),8X,2 (D12.5,2X))
 9960 FORMAT (' ',I4,2X,2 (D12.5,2X),8X,2 (D12.5,2X),
     *       1 (/' ',6X,2 (D12.5,2X),8X,2 (D12.5,2X)))
 9950 FORMAT (' ',I4,2X,2 (D12.5,2X),8X,2 (D12.5,2X),
     *       2 (/' ',6X,2 (D12.5,2X),8X,2 (D12.5,2X)))
 9940 FORMAT (' ',I4,2X,2 (D12.5,2X),8X,2 (D12.5,2X),
     *       3 (/' ',6X,2 (D12.5,2X),8X,2 (D12.5,2X)))
 9930 FORMAT (' ',I4,2X,2 (D12.5,2X),8X,2 (D12.5,2X),
     *       4 (/' ',6X,2 (D12.5,2X),8X,2 (D12.5,2X)))
 9920 FORMAT (/' ',2 ('NODE',6X,'REAL',7X,'IMAGINARY',6X),/' ')
 9910 FORMAT (2 (' ',I4,2X,2 (D12.5,2X),1X))
      END
