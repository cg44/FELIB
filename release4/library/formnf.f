C
      SUBROUTINE FORMNF(REST,IREST,JREST,RESNOD,TOTNOD,DOFNOD,NF,INF,
     *                  JNF,TOTDOF,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMNF constructs the nodal freedom array from the restrained
C      freedom data
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
C      Release 1.1  29 Oct 1979 (IMS)
C      Commented    18 Feb 1980 (KR)
C      Recoded      01 Nov 1981 (NB)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      REST    integer array; REST(I ,J) contains the I'TH set
C              of restraint information - REST(I, 1) contains
C              the node number, REST(I, J+1) for J=1(1)DOFNOD
C              contains the local freedom numbers of the
C              freedoms that are restrained
C      IREST   first dimension of array REST (.GE. RESNOD)
C      JREST   second dimension of REST (.GE. DOFNOD)
C      RESNOD  number of nodes at which freedoms are restrained
C      TOTNOD  total number of nodes in mesh
C      INF     first dimension of array INF (.GE. TOTNOD)
C      JNF     second dimension of INF (.GE. DOFNOD)
C      ITEST   error checking option
C
C ARGUMENTS out
C      NF      NF(I, J), J=1(1)DOFNOD, contains the freedom
C              numbers associated with the I'TH node
C      TOTDOF  total number of freedoms in problem under
C              consideration
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE FORMNF(REST,IREST,JREST,RESNOD,TOTNOD,DOFNOD,NF,INF,
C    *                  JNF,TOTDOF,ITEST)
C***********************************************************************
C
      INTEGER DOFNOD,ERRMES,I,IERROR,INF,IREST,ITEST,J,JNF,JREST,
     *        JTEST,K,L,M,NF,RESNOD,REST,TOTDOF,TOTNOD
      CHARACTER*6 SRNAME
      LOGICAL SWITCH
      DIMENSION NF(INF,JNF),REST(IREST,JREST)
C
      INTRINSIC ABS
      EXTERNAL ERRMES
C
      DATA SRNAME/'FORMNF'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (INF.LT.TOTNOD .OR. JNF.LT.DOFNOD) IERROR = 3
         IF (IREST.LT.RESNOD .OR. JREST.LT.DOFNOD+1) IERROR = 2
         IF (RESNOD.LT.0 .OR. TOTNOD.LE.0 .OR. DOFNOD.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      SWITCH = .TRUE.
C
C     Initialise nodel freedom array
C
      DO 1010 I = 1,TOTNOD
         DO 1000 J = 1,DOFNOD
            NF(I,J) = 1
 1000    CONTINUE
 1010 CONTINUE
C
C     If no restrained nodes branch
C
      IF (RESNOD.NE.0) THEN
         DO 1030 I = 1,RESNOD
            K = REST(I,1)
            DO 1020 J = 1,DOFNOD
               L = REST(I,J+1)
               M = ABS(L)
C
C     Range checking on K and M
C
               IF (JTEST.NE.-1) THEN
                  IERROR = 0
                  IF (K.GT.TOTNOD .OR. M.GT.DOFNOD) IERROR = 4
                  ITEST = ERRMES(JTEST,IERROR,SRNAME)
                  IF (ITEST.NE.0) RETURN
               END IF
C
               IF (L.GT.0) NF(K,L) = 0
               IF (L.LT.0) THEN
                  NF(K,M) = L
                  SWITCH = .FALSE.
               END IF
 1020       CONTINUE
 1030    CONTINUE
      END IF
C
C     Renumber nodal freedom array
C
      K = 1
      DO 1050 I = 1,TOTNOD
         DO 1040 J = 1,DOFNOD
            IF (NF(I,J).GT.0) THEN
               NF(I,J) = K
               K = K + 1
            END IF
 1040    CONTINUE
 1050 CONTINUE
C
C     Set up for prescribed values
C
      TOTDOF = K - 1
      IF (.NOT.SWITCH) THEN
         DO 1070 I = 1,TOTNOD
            DO 1060 J = 1,DOFNOD
               IF (NF(I,J).LT.0) THEN
                  NF(I,J) = -K
                  K = K + 1
               END IF
 1060       CONTINUE
 1070    CONTINUE
      END IF
      END
