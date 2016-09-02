C***********************************************************************
C$SPLIT$DCSQUA$*********************************************************
C***********************************************************************
      SUBROUTINE DCSQUA(JACIN, IJACIN, JJACIN, SIDNUM, COSIN,
     *     ICOSIN, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES THE DIRECTION COSINES OF THE OUTWARD NORMAL
C      TO THE SPECIFIED SIDE OF A QUADRILATERAL ELEMENT GIVEN THE
C      JACOBIAN INVERSE
C
C HISTORY
C      RELEASE 2.0  1 FEB 1981 (CG) --- SERC COPYRIGHT
C      COMMENTED    1 FEB 1981 (CG)
C
C ARGUMENTS IN
C      JACIN   ARRAY OF DIMENSION (IJACIN,JJACIN). CONTAINS THE
C              INVERSE OF THE TRANSFORMATION JACOBIAN.
C      IJACIN  FIRST DIMENSION OF JACIN (.GE.2)
C      JJACIN  SECOND DIMENSION OF JACIN (.GE.2)
C      SIDNUM  THE SIDE NUMBER FOR WHICH THE OUTWARD NORMAL
C              IS REQUIRED (.LE.4)
C      ICOSIN  DIMENSION OF VECTOR COSIN (.GE.2)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      COSIN   VECTOR OF DIMENSION ICOSIN. CONTAINS THE DIR-
C              ECTION COSINES OF THE OUTWARD NORMAL
C
C ROUTINES CALLED
C      ERRMES  MATVEC
C
C
C     SUBROUTINE DCSQUA(JACIN,IJACIN,JJACIN,SIDNUM,COSIN,ICOSIN,ITEST)
C***********************************************************************
C
      INTEGER DIMEN, ERRMES, ICOSIN, IERROR, IJACIN, ITEST,
     *     IWORK, JJACIN, SIDNUM
      DOUBLE PRECISION AMOD, COSIN, JACIN, SRNAME, WORK
      DIMENSION COSIN(ICOSIN), JACIN(IJACIN,JJACIN), WORK(2)
      DATA DIMEN /2/, IWORK /2/, SRNAME /8H DCSQUA /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IJACIN.LT.2 .OR. JJACIN.LT.2)
     *                      IERROR = 3
                        IF (ICOSIN.LT.2) IERROR = 2
                        IF (SIDNUM.LE.0 .OR. SIDNUM.GT.4)
     *                      IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 GO TO (1020, 1030, 1040, 1050), SIDNUM
 1020 WORK(1) = -1.0D0
      WORK(2) = 0.0D0
      GO TO 1060
 1030 WORK(1) = 0.0D0
      WORK(2) = 1.0D0
      GO TO 1060
 1040 WORK(1) = 1.0D0
      WORK(2) = 0.0D0
      GO TO 1060
 1050 WORK(1) = 0.0D0
      WORK(2) = -1.0D0
 1060 CALL MATVEC(JACIN, IJACIN, JJACIN, WORK, IWORK, DIMEN,
     *     DIMEN, COSIN, ICOSIN, ITEST)
      AMOD = DSQRT(COSIN(1)*COSIN(1)+COSIN(2)*COSIN(2))
      COSIN(1) = COSIN(1)/AMOD
      COSIN(2) = COSIN(2)/AMOD
      RETURN
      END
