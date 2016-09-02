C
      SUBROUTINE UPDATE(PHI,IPHI,RHS,IRHS,TOTNOD,DOFNOD,TOTDOF,NF,INF,
     *                  JNF,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      UPDATE takes a full solution vector and a set of
C      updates and updates the solution vector
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
CC
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307

C      Release 2.0  29 Jun 1986 (CJH, CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      PHI     solution vector
C      IPHI    dimension of vector PHI (IPHI .GE. TOTNOD*DOFNOD)
C      RHS     vector of updates
C      IRHS    dimension of vector RHS (IRHS .GE. TOTDOF)
C      TOTNOD  the number of nodes in the problem
C      DOFNOD  the maximum number of nodes per node
C      TOTDOF  the number of freedoms in RHS
C      NF      the nodal freedom array
C      INF     first dimension of NF (INF .GE. TOTNOD)
C      JNF     second dimension of NF (JNF .GE. DOFNOD)
C      ITEST   error checking option
C
C ARGUMENTS out
C      PHI     the updated solution vector
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE UPDATE(PHI,IPHI,RHS,IRHS,TOTNOD,DOFNOD,TOTDOF,NF,INF,
C    *                  JNF,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,DOFNOD,I,IERROR,INF,IPHI,IRHS,ITEST,J,JNF,K,L,NF,
     *        TOTDOF,TOTNOD
      CHARACTER*6 SRNAME
      DOUBLE PRECISION PHI,RHS
      DIMENSION NF(INF,JNF),PHI(IPHI),RHS(IRHS)
      DATA SRNAME/'UPDATE'/
C
C     Check array bounds
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IPHI.LT.TOTNOD*DOFNOD) IERROR = 1
         IF (IRHS.LT.TOTDOF) IERROR = 2
         IF (INF.LT.TOTNOD) IERROR = 3
         IF (JNF.LT.DOFNOD) IERROR = 4
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      DO 1010 I = 1,TOTNOD
C
         DO 1000 J = 1,DOFNOD
C
            K = NF(I,J)
C
C     Don't UPDATE restrained variables
C
            IF (K.NE.0) THEN
C
               L = I + J - 1
C
C     UPDATE solution vector
C
               PHI(L) = PHI(L) + RHS(K)
            END IF
 1000    CONTINUE
 1010 CONTINUE
C
C
      END
