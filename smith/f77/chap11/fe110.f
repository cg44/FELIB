      PROGRAM P110
C
C      PROGRAM 11.0 FORCED VIBRATION OF A RECTANGULAR SOLID IN
C      PLANE STRAIN USING 8-NODE QUADRILATERALS,LUMPED MASS,
C      MODAL SUPERPOSITION
C
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER(IBIGK=103,INF=85,IMOD=30)
C
      REAL DEE(3,3),SAMP(3,2),COORD(8,2),FUN(8),JAC(2,2),JAC1(2,2),
     +DER(2,8),DERIV(2,8),BEE(3,16),DBEE(3,16),BTDB(16,16),
     +BT(16,3),KM(16,16),BIGK(IBIGK,IBIGK),XMOD(IMOD),
     +LOADS(IBIGK),DIAG(IBIGK),UDIAG(IBIGK),EMM(16,16)
      INTEGER NF(INF,2),G(16)
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/
      DATA IH,ISAMP,IDEE,IBEE,IDBEE/5*3/
      DATA ICOORD,NOD/2*8/,IBTDB,IKM,IBT,IEMM,IDOF/5*16/
C
C      INPUT AND INITIALISATION
C
      READ(5,*)NXE,NYE,N,NN,NR,NGP,AA,BB,RHO,E,V,DR,NMODES,
     +ISTEP,OMEGA,NPRI
      CALL READNF(NF,INF,NN,NODOF,NR)
      PI=ACOS(-1.)
      PERIOD=2.*PI/OMEGA
      DTIM=PERIOD/20.
      CALL NULL(BIGK,IBIGK,N,N)
      CALL NULL(DEE,IDEE,IH,IH)
      CALL FMDEPS(DEE,IDEE,E,V)
      CALL GAUSS(SAMP,ISAMP,NGP)
      CALL NULL(EMM,IEMM,IDOF,IDOF)
      CALL NULVEC(DIAG,N)
C
C      FORM LUMPED MASS MATRIX
C
      DO 30 I=1,IDOF
   30 EMM(I,I)=AA*BB*RHO*.2
      DO 40 I=1,13,4
   40 EMM(I,I)=EMM(3,3)*.25
      DO 50 I=2,14,4
   50 EMM(I,I)=EMM(3,3)*.25
C
C      ELEMENT STIFFNESS AND MASS INTEGRATION AND ASSEMBLY
C
      DO 10 IP=1,NXE
      DO 10 IQ=1,NYE
      CALL GEOM8Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)
      CALL NULL(KM,IKM,IDOF,IDOF)
      DO 20 I=1,NGP
      DO 20 J=1,NGP
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)
      CALL NULL(BEE,IBEE,IH,IDOF)
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)
      QUOT=DET*SAMP(I,2)*SAMP(J,2)
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)
   20 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)
      CALL FMLUMP(DIAG,EMM,IEMM,G,IDOF)
   10 CALL FMBIGK(BIGK,IBIGK,KM,IKM,G,IDOF)
C
C      REDUCE TO STANDARD EIGENVALUE PROBLEM
C
      DO 60 I=1,N
      DIAG(I)=1./SQRT(DIAG(I))
   60 LOADS(I)=DIAG(I)
      DO 70 I=1,N
      DO 70 J=1,N
   70 BIGK(I,J)=BIGK(I,J)*DIAG(I)*DIAG(J)
C
C      EXTRACT EIGENVALUES
C
      CALL TRIDIA(N,1.E-20,BIGK,DIAG,UDIAG,BIGK,IBIGK)
      IFAIL=1
      CALL EVECTS(N,1.E-20,DIAG,UDIAG,BIGK,IBIGK,IFAIL)
C
C      TIME STEPPING LOOP
C
      TIM=0.
      DO 80 J=1,ISTEP
      TIM=TIM+DTIM
      DO 90 M=1,NMODES
      DO 100 I=1,N
  100 UDIAG(I)=BIGK(I,M)*LOADS(I)
      F=UDIAG(N)
      X1=DIAG(M)-OMEGA*OMEGA
      X2=X1*X1+4.*OMEGA*OMEGA*DR*DR*DIAG(M)
      X3=F*X1/X2
      X4=F*2.*OMEGA*DR*SQRT(DIAG(M))/X2
      XMOD(M)=X3*COS(OMEGA*TIM)+X4*SIN(OMEGA*TIM)
   90 CONTINUE
C
C      SUPERIMPOSE THE MODES
C
      DO 110 I=1,N
      SUM=0.
      DO 120 M=1,NMODES
  120 SUM=SUM+BIGK(I,M)*LOADS(I)*XMOD(M)
  110 UDIAG(I)=SUM
      IF(J/NPRI*NPRI.EQ.J)WRITE(6,1000)TIM,COS(OMEGA*TIM),UDIAG(N)
   80 CONTINUE
 1000 FORMAT(3E12.4)
      END
