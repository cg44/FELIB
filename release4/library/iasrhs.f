C
      SUBROUTINE IASRHS(RHS,IRHS,VALUE,IVALUE,STEER,ISTEER,DOFEL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      IASRHS the routine adds into the right-hand side of a system
C      the values contianed in an element vector, thus
C      forming the right-hand side VALUE is real, RHS is complex
C      (ordered pairs), with VALUE added into imaginary part of RHS.
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
C      Release 2.0   1 Jul 1984 (CRIE)
C      Commented     1 Jul 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      RHS     the right-hand side of the system
C              complex ( ordered pairs )
C      IRHS    dimension of array rsh
C      VALUE   the element vector of the current element to
C              be added into the right-hand side, imaginary component
C      IVALUE  dimension of array VALUE
C      STEER   the steering vector containing the freedom
C              numbers of the freedoms associated with the
C              current element in the local order
C      ISTEER  dimension of array STEER
C      DOFEL   the maximum number of degrees of freedom on
C              an element of the current type
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE IASRHS(RHS,IRHS,VALUE,IVALUE,STEER,ISTEER,DOFEL,ITEST)
C***********************************************************************
C
      INTEGER DOFEL,ERRMES,IERROR,IRHS,ISTEER,ITEST,IVALUE,K,L,STEER,
     *        JTEST
      CHARACTER*6 SRNAME
      DOUBLE PRECISION RHS,VALUE
      DIMENSION RHS(2,IRHS),STEER(ISTEER),VALUE(IVALUE)
C
      EXTERNAL ERRMES
C
      DATA SRNAME/'IASRHS'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (ISTEER.LT.DOFEL) IERROR = 3
         IF (IVALUE.LT.DOFEL) IERROR = 2
         IF (DOFEL.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      DO 1000 K = 1,DOFEL
         IF (STEER(K).NE.0) THEN
            L = STEER(K)
            IF (JTEST.NE.-1) THEN
               IERROR = 0
               IF (L.GT.IRHS) IERROR = 4
               ITEST = ERRMES(ITEST,IERROR,SRNAME)
               IF (ITEST.NE.0) RETURN
            END IF
            RHS(2,L) = RHS(2,L) + VALUE(K)
         END IF
 1000 CONTINUE
C
      END
