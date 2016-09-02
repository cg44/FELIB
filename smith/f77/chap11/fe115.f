      PROGRAM P115
C
C      PROGRAM 11.5 FORCED VIBRATION OF A RECTANGULAR
C      SOLID IN PLANE STRAIN USING 4-NODE QUADRILATERALS
C      IMPLICIT AND EXPLICIT INTEGRATION
C      LUMPED OR CONSISTENT MASS
C
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE
C
      PARAMETER(IKV=1000,ILOADS=100,INF=100,INX=25,INY=25)
C
      REAL DEE(3,3),SAMP(3,2),COORD(4,2),JAC(2,2),JAC1(2,2),BTDB(8,8),
     +DER(2,4),DERIV(2,4),BEE(3,8),DBEE(3,8),KM(8,8),EMM(8,8),
     +ECM(8,8),FUN(4),BT(8,3),TN(8,8),NT(8,2),KV(IKV),MM(IKV),
     +LOADS(ILOADS),X0(ILOADS),D1X0(ILOADS),D2X0(ILOADS),
     +X1(ILOADS),D1X1(ILOADS),D2X1(ILOADS)
      INTEGER NF(INF,2),KDIAG(ILOADS),TYPE(INX,INY),G(8)
      DATA IBTDB,IKM,IEMM,IECM,IBT,ITN,INT,IDOF/8*8/,ICOORD,NOD/2*4/
      DATA ISAMP,IDEE,IBEE,IDBEE,IH/5*3/
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/
C
C      INPUT AND INITIALISATION
C
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,RHO,
     +         ISTEP,DTIM,AA,BB,V,E,GAMMA,BETA
      READ(5,*)((TYPE(I,J),J=1,NYE),I=1,NXE)
      CALL READNF(NF,INF,NN,NODOF,NR)
      IR=N*(IW+1)
      C1=1./DTIM/DTIM/BETA
      C2=GAMMA/DTIM/BETA
      CALL NULVEC(MM,IR)
      CALL NULL(DEE,IDEE,IH,IH)
      CALL FMDEPS(DEE,IDEE,E,V)
      CALL GAUSS(SAMP,ISAMP,NGP)
C
C      VARIABLE BANDWIDTH STORE
C
      DO 10 I=1,N
   10 KDIAG(I)=0
      DO 20 IP=1,NXE
      DO 20 IQ=1,NYE
      CALL GEOM4Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)
      IF(TYPE(IP,IQ).NE.1)CALL FKDIAG(KDIAG,G,IDOF)
   20 CONTINUE
      DO 30 I=1,N
   30 IF(KDIAG(I).EQ.0)KDIAG(I)=1
      KDIAG(1)=1
      DO 40 I=2,N
   40 KDIAG(I)=KDIAG(I)+KDIAG(I-1)
      IR=KDIAG(N)
      CALL NULVEC(KV,IR)
C
C      ELEMENT STIFFNESS AND MASS INTEGRATION AND ASSEMBLY
C
      DO 50 IP=1,NXE
      DO 50 IQ=1,NYE
      AREA=0.
      CALL GEOM4Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)
      CALL NULL(KM,IKM,IDOF,IDOF)
      CALL NULL(EMM,IEMM,IDOF,IDOF)
      IF(TYPE(IP,IQ).EQ.1)THEN
      DO 60 I=1,IDOF
   60 EMM(I,I)=1.
      END IF
      DO 70 I=1,NGP
      DO 70 J=1,NGP
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)
      CALL NULL(BEE,IBEE,IH,IDOF)
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)
      IF(TYPE(IP,IQ).NE.1)
     +CALL ECMAT(ECM,IECM,TN,ITN,NT,INT,FUN,NOD,NODOF)
      QUOT=DET*SAMP(I,2)*SAMP(J,2)
      AREA=AREA+QUOT
      DO 80 K=1,IDOF
      DO 80 L=1,IDOF
      BTDB(K,L)=BTDB(K,L)*QUOT
      IF(TYPE(IP,IQ).NE.1)ECM(K,L)=ECM(K,L)*QUOT*RHO*C1
   80 CONTINUE
      CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)
      IF(TYPE(IP,IQ).NE.1)CALL MATADD(EMM,IEMM,ECM,IECM,IDOF,IDOF)
   70 CONTINUE
      AREA=AREA/NOD*RHO
      IF(TYPE(IP,IQ).EQ.1)THEN
      DO 90 K=1,IDOF
   90 EMM(K,K)=EMM(K,K)*AREA*C1
      CALL FDIAGV(KV,EMM,IEMM,G,KDIAG,IDOF)
      DO 100 I=1,IDOF
      DO 100 J=1,IDOF
  100 KM(I,J)=-KM(I,J)
      CALL FORMKV(MM,KM,IKM,G,N,IDOF)
      CALL FORMKV(MM,EMM,IEMM,G,N,IDOF)
      ELSE
      CALL FSPARV(KV,KM,IKM,G,KDIAG,IDOF)
      CALL FSPARV(KV,EMM,IEMM,G,KDIAG,IDOF)
      CALL FORMKV(MM,EMM,IEMM,G,N,IDOF)
      END IF
   50 CONTINUE
C
C      TIME INTEGRATION BY SIMPLE PREDICTOR-CORRECTOR
C
      TIM=0.
      DO 110 JR=1,N
      X0(JR)=0.
      D1X0(JR)=1.
  110 D2X0(JR)=0.
      CALL SPARIN(KV,N,KDIAG)
      DO 120 J=1,ISTEP
      TIM=TIM+DTIM
      DO 130 JR=1,N
      LOADS(JR)=.0
  130 D1X1(JR)=X0(JR)+D1X0(JR)*DTIM+D2X0(JR)*.5*DTIM*DTIM*
     +(1.-2.*BETA)
      CALL LINMUL(MM,D1X1,X1,N,IW)
      CALL VECADD(X1,LOADS,X1,N)
      CALL SPABAC(KV,X1,N,KDIAG)
      DO 140 JR=1,N
      D2X1(JR)=(X1(JR)-D1X1(JR))/DTIM/DTIM/BETA
  140 D1X1(JR)=D1X0(JR)+D2X0(JR)*DTIM*(1.-GAMMA)+D2X1(JR)*DTIM*
     +GAMMA
      WRITE(6,1000)TIM,X1(N),D1X1(N),D2X1(N)
      CALL VECCOP(X1,X0,N)
      CALL VECCOP(D1X1,D1X0,N)
      CALL VECCOP(D2X1,D2X0,N)
  120 CONTINUE
 1000 FORMAT(F8.5,4E12.4)
      STOP
      END
