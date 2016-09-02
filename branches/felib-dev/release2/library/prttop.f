C***********************************************************************
C$SPLIT$PRTTOP$*********************************************************
C***********************************************************************
      SUBROUTINE PRTTOP(TOTELS, ELTOP, IELTOP, JELTOP, NOUT, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRINTS ELEMENT TOPOLOGIES IN A STANDARD FORMAT
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    14 OCT 1980 (KR)
C
C ARGUMENTS IN
C      TOTELS  TOTAL NUMBER OF ELEMENTS IN THE MESH
C      ELTOP   INTEGER ARRAY OF DIMENSION (IELTOP,JELTOP)
C              CONTAINING ELEMENT TOPOLOGIES, ELEMENT TYPE, AND
C              NUMBER OF NODES ON THE ELEMENT
C      IELTOP  FIRST DIMENSION OF ELTOP (.GE.TOTELS)
C      JELTOP  SECOND DIMENSION OF ELTOP (.GE.NUMBER OF NODES
C              ON ELEMENT + 2)
C      NOUT    FORTRAN UNIT NUMBER
C      ITEST   ERROR CHECKING OPTION
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE PRTTOP(TOTELS,ELTOP,IELTOP,JELTOP,ITEST)
C***********************************************************************
C
      INTEGER ELTOP, ERRMES, I, IELTOP, IERROR, ITEST, J, JELTOP,
     *     K, L, NOUT, TOTELS
      DOUBLE PRECISION SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP)
      DATA SRNAME /8H PRTTOP /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IELTOP.LT.TOTELS) IERROR = 2
                        IF (TOTELS.LE.0) IERROR = 1
                        ITEST = ERRMES(IERROR,ITEST,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 WRITE (NOUT,9010)
      WRITE (NOUT,9020) TOTELS
      WRITE (NOUT,9030)
      DO 1030 I=1,TOTELS
      L = ELTOP(I,2)
      K = L + 2
                        IF (ITEST.EQ.-1) GO TO 1020
                        IERROR = 0
                        IF (JELTOP.LT.K) IERROR = 3
                        ITEST = ERRMES(IERROR,ITEST,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1020 WRITE (NOUT,9040) I, (ELTOP(I,J),J=1,K)
 1030 CONTINUE
      RETURN
 9010 FORMAT (1H ////27H **** ELEMENT TOPOLOGY ****//1H )
 9020 FORMAT (22H NUMBER OF ELEMENTS = , I3)
 9030 FORMAT (/1H , 2X, 4HELEM, 4X, 5HELTYP, 4X, 5HNODEL, 4X,
     *     5HNODES/1H )
 9040 FORMAT (1H , 2X, I4, 4X, I3, 7X, I3, 4X, 9(I4, 2X)/1H , 27X,
     *     9(I4, 2X))
      END
