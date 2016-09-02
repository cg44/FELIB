      PROGRAM P113
C
C      PROGRAM 11.3 FORCED VIBRATION OF A RECTANGULAR SOLID IN
C      PLANE STRAIN USING 8-NODE QUADRILATERALS,
C      LUMPED MASS,COMPLEX RESPONSE METHOD
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER(IKC=1000,ILOADS=103,INF=85)
C
      REAL DEE(3,3),SAMP(3,2),COORD(8,2),FUN(8),JAC(2,2),JAC1(2,2),
     +DER(2,8),DERIV(2,8),BEE(3,16),DBEE(3,16),BTDB(16,16),
     +BT(16,3),KM(16,16),CM(16,16),EMM(16,16)
      COMPLEX KC(IKC),LOADS(ILOADS)
      INTEGER NF(INF,2),G(16)
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/
      DATA IH,ISAMP,IDEE,IBEE,IDBEE/5*3/
      DATA ICOORD,NOD/2*8/,IBTDB,IKM,IBT,IEMM,IDOF/5*16/
C
C      INPUT AND INITIALISATION
C
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,AA,BB,RHO,E,V,DR,OMEGA,
     +ISTEP,NPRI
      CALL READNF(NF,INF,NN,NODOF,NR)
      IR=N*(IW+1)
      PI=ACOS(-1.)
      PERIOD=2.*PI/OMEGA
      DTIM=PERIOD/20.
      CALL FMDEPS(DEE,IDEE,E,V)
      CALL GAUSS(SAMP,ISAMP,NGP)
C
C      SOLVE FOR SINGLE INPUT LOADING HARMONIC
C
      CALL NULL(EMM,IEMM,IDOF,IDOF)
      DO 10 I=1,IR
   10 KC(I)=(0.,0.)
      DO 20 I=1,N
   20 LOADS(I)=(0.,0.)
C
C      FORM ELEMENT LUMPED MASS METRIX
C
      DO 30 I=1,IDOF
   30 EMM(I,I)=AA*BB*RHO*.2*OMEGA**2
      DO 31 I=1,13,4
   31 EMM(I,I)=EMM(3,3)*.25
      DO 32 I=2,14,4
   32 EMM(I,I)=EMM(3,3)*.25
C
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY
C
      DO 40 IP=1,NXE
      DO 40 IQ=1,NYE
      CALL GEOM8Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)
      CALL NULL(KM,IKM,IDOF,IDOF)
      DO 50 I=1,NGP
      DO 50 J=1,NGP
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)
      CALL NULL(BEE,IBEE,IH,IDOF)
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)
      QUOT=DET*SAMP(I,2)*SAMP(J,2)*(1.-2.*DR*DR)
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)
   50 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)
C
C      COMPLEX MATRIX ASSEMBLY
C
      SNP=2.*DR*SQRT(1.-DR*DR)/(1.-2.*DR*DR)
      DO 60 I=1,IDOF
      DO 60 J=1,IDOF
      CM(I,J)=KM(I,J)*SNP
   60 KM(I,J)=KM(I,J)-EMM(I,J)
   40 CALL FORMKC(KC,KM,IKM,CM,IKM,G,N,IDOF)
C
C      COMPLEX EQUATION SOLUTION
C
      LOADS(N)=(1.,0.)
      CALL COMRED(KC,N,IW)
      CALL COMBAC(KC,LOADS,N,IW)
      A=REAL(LOADS(N))
      B=AIMAG(LOADS(N))
      TIM=0.
      DO 70 J=1,ISTEP
      TIM=TIM+DTIM
      IF(J/NPRI*NPRI.EQ.J)
     +WRITE(6,1000)TIM,A*COS(OMEGA*TIM)-B*SIN(OMEGA*TIM)
   70 CONTINUE
 1000 FORMAT(3E12.4)
      STOP
      END
