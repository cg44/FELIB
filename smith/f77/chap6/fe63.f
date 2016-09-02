      PROGRAM P63                                                               
C                                                                               
C      PROGRAM 6.3 AXISYMMETRIC STRAIN OF AN UNDRAINED                          
C      ELASTIC-PLASTIC (MOHR-COULOMB) SOLID USING                               
C      8-NODED QUADRILATERAL ELEMENTS,VISCOPLASTIC STRAIN METHOD                
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=3500,ILOADS=100,INF=100,                                    
     +          IEVPT=500,INR=20,IND=20,INO=10)                                 
C                                                                               
      REAL DEE(4,4),SAMP(4,2),COORD(8,2),JAC(2,2),JAC1(2,2),                    
     +DER(2,8),DERIV(2,8),BEE(4,16),DBEE(4,16),RAD(INR),DEP(IND),               
     +BTDB(16,16),KM(16,16),ELD(16),EPS(4),SIGMA(4),BLOAD(16),                  
     +BT(16,4),FUN(8),KV(IKV),LOADS(ILOADS),ELOAD(16),                          
     +TOTD(ILOADS),BDYLDS(ILOADS),EVPT(IEVPT),OLDIS(ILOADS),                    
     +SX(INR,IND,4),SY(INR,IND,4),TXY(INR,IND,4),SZ(INR,IND,4),                 
     +EX(INR,IND,4),EY(INR,IND,4),GXY(INR,IND,4),EZ(INR,IND,4),                 
     +STORKV(INO),ERATE(4),EVP(4),DEVP(4),M1(4,4),M2(4,4),M3(4,4),              
     +FLOW(4,4),STRESS(4),PORE(INR,IND,4)                                       
      INTEGER NF(INF,2),G(16),NO(INO)                                           
      DATA IDEE,IBEE,IDBEE,IH,IFLOW/5*4/,IDOF,IBTDB,IBT,IKM/4*16/               
      DATA IJAC,IJAC1,NODOF,IT,IDER,IDERIV/6*2/,ICOORD,NOD/2*8/                 
      DATA ISAMP/4/                                                             
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)PHI,C,PSI,E,V,BW,CONS,NRE,NDE,N,IW,NN,NR,NGP                     
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      READ(5,*)(RAD(I),I=1,NRE+1)                                               
      READ(5,*)(DEP(I),I=1,NDE+1)                                               
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL NULVEC(OLDIS,N)                                                      
      CALL NULVEC(TOTD,N)                                                       
      CALL FMDRAD(DEE,IDEE,E,V)                                                 
C                                                                               
C      ADD FLUID BULK MODULUS                                                   
C                                                                               
      DO 10 I=1,IH                                                              
      DO 10 J=1,IH                                                              
      IF(I.EQ.3.OR.J.EQ.3)GOTO 10                                               
      DEE(I,J)=DEE(I,J)+BW                                                      
   10 CONTINUE                                                                  
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      SNPH=SIN(PHI*4.*ATAN(1.)/180.)                                            
      DT=4.*(1.+V)*(1.-2.*V)/(E*(1.-2.*V+SNPH*SNPH))                            
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 20 IP=1,NRE                                                            
      DO 20 IQ=1,NDE                                                            
      CALL GEOV8Y(IP,IQ,NDE,RAD,DEP,COORD,ICOORD,G,NF,INF)                      
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      IG=0                                                                      
      DO 30 I=1,NGP                                                             
      DO 30 J=1,NGP                                                             
      IG=IG+1                                                                   
      SX(IP,IQ,IG)=CONS                                                         
      SY(IP,IQ,IG)=CONS                                                         
      TXY(IP,IQ,IG)=0.                                                          
      SZ(IP,IQ,IG)=CONS                                                         
      EX(IP,IQ,IG)=0.                                                           
      EY(IP,IQ,IG)=0.                                                           
      GXY(IP,IQ,IG)=0.                                                          
      EZ(IP,IQ,IG)=0.                                                           
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FMBRAD(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)               
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=SUM*DET*SAMP(I,2)*SAMP(J,2)                                          
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
   30 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   20 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      READ PRESCRIBED DISPLACEMENTS                                            
C      AUGMENT AND REDUCE STIFFNESS MATRIX                                      
C                                                                               
      READ(5,*)NL,(NO(I),I=1,NL),PRESC,INCS,ITS                                 
      DO 40 I=1,NL                                                              
      KV(NO(I))=KV(NO(I))+1.E20                                                 
   40 STORKV(I)=KV(NO(I))                                                       
      CALL BANRED(KV,N,IW)                                                      
C                                                                               
C      DISPLACEMENT INCREMENT LOOP                                              
C                                                                               
      CALL FMDRAD(DEE,IDEE,E,V)                                                 
      DO 50 IY=1,INCS                                                           
      PTOT=PRESC*IY                                                             
      ITERS=0                                                                   
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(EVPT,NRE*NDE*IH*NGP*NGP)                                      
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
   60 ITERS=ITERS+1                                                             
      CALL NULVEC(LOADS,N)                                                      
      DO 70 I=1,NL                                                              
   70 LOADS(NO(I))=STORKV(I)*PRESC                                              
      CALL VECADD(LOADS,BDYLDS,LOADS,N)                                         
      CALL BACSUB(KV,LOADS,N,IW)                                                
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(LOADS,OLDIS,N,0.0001,ICON)                                    
      IF(ITERS.EQ.1)ICON=0                                                      
C                                                                               
C      INSPECT ALL GAUSS POINTS                                                 
C                                                                               
      NM=0                                                                      
      DO 80 IP=1,NRE                                                            
      DO 80 IQ=1,NDE                                                            
      NM=NM+1                                                                   
      CALL NULVEC(BLOAD,IDOF)                                                   
      CALL GEOV8Y(IP,IQ,NDE,RAD,DEP,COORD,ICOORD,G,NF,INF)                      
      DO 90 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   90 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      IG=0                                                                      
      DO 100 I=1,NGP                                                            
      DO 100 J=1,NGP                                                            
      IG=IG+1                                                                   
      IN=NGP*NGP*IH*(NM-1)+IH*(IG-1)                                            
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FMBRAD(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)               
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      DO 110 K=1,IH                                                             
  110 EPS(K)=EPS(K)-EVPT(IN+K)                                                  
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      STRESS(1)=SIGMA(1)+SX(IP,IQ,IG)                                           
      STRESS(2)=SIGMA(2)+SY(IP,IQ,IG)                                           
      STRESS(3)=SIGMA(3)+TXY(IP,IQ,IG)                                          
      STRESS(4)=SIGMA(4)+SZ(IP,IQ,IG)                                           
      CALL INVAR(STRESS,SIGM,DSBAR,THETA)                                       
C                                                                               
C      CHECK WHETHER YIELD IS VIOLATED                                          
C                                                                               
      CALL MOCOUF(PHI,C,SIGM,DSBAR,THETA,F)                                     
      IF(F.LT.0.)GOTO 120                                                       
      CALL MOCOUQ(PSI,DSBAR,THETA,DQ1,DQ2,DQ3)                                  
      CALL FORMM(STRESS,M1,M2,M3)                                               
      DO 130 L=1,IH                                                             
      DO 130 M=1,IH                                                             
  130 FLOW(L,M)=F*(M1(L,M)*DQ1+M2(L,M)*DQ2+M3(L,M)*DQ3)                         
      CALL MVMULT(FLOW,IFLOW,STRESS,IH,IH,ERATE)                                
      DO 140 K=1,IH                                                             
      EVP(K)=ERATE(K)*DT                                                        
  140 EVPT(IN+K)=EVPT(IN+K)+EVP(K)                                              
      CALL MVMULT(DEE,IDEE,EVP,IH,IH,DEVP)                                      
      CALL MVMULT(BT,IBT,DEVP,IDOF,IH,ELOAD)                                    
      QUOT=SUM*DET*SAMP(I,2)*SAMP(J,2)                                          
      DO 150 K=1,IDOF                                                           
  150 BLOAD(K)=BLOAD(K)+ELOAD(K)*QUOT                                           
  120 IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 100                                    
C                                                                               
C      UPDATE GAUSS POINT STRESSES AND STRAINS                                  
C      CALCULATE PORE PRESSURES                                                 
C                                                                               
      SX(IP,IQ,IG)=STRESS(1)                                                    
      SY(IP,IQ,IG)=STRESS(2)                                                    
      TXY(IP,IQ,IG)=STRESS(3)                                                   
      SZ(IP,IQ,IG)=STRESS(4)                                                    
      EX(IP,IQ,IG)=EX(IP,IQ,IG)+EPS(1)+EVPT(IN+1)                               
      EY(IP,IQ,IG)=EY(IP,IQ,IG)+EPS(2)+EVPT(IN+2)                               
      GXY(IP,IQ,IG)=GXY(IP,IQ,IG)+EPS(3)+EVPT(IN+3)                             
      EZ(IP,IQ,IG)=EZ(IP,IQ,IG)+EPS(4)+EVPT(IN+4)                               
      PORE(IP,IQ,IG)=(EX(IP,IQ,IG)+EY(IP,IQ,IG)+EZ(IP,IQ,IG))*BW                
  100 CONTINUE                                                                  
C                                                                               
C      COMPUTE TOTAL BODYLOADS VECTOR                                           
C                                                                               
      DO 160 M=1,IDOF                                                           
      IF(G(M).EQ.0)GOTO 160                                                     
      BDYLDS(G(M))=BDYLDS(G(M))+BLOAD(M)                                        
  160 CONTINUE                                                                  
   80 CONTINUE                                                                  
      IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 60                                     
      CALL VECADD(TOTD,LOADS,TOTD,N)                                            
      WRITE(6,1000)PTOT                                                         
      WRITE(6,1000)SX(1,1,1),SY(1,1,1),SZ(1,1,1)                                
      WRITE(6,1000)DSBAR,PORE(1,1,1)                                            
      WRITE(6,2000)ITERS                                                        
      IF(ITERS.EQ.ITS)GOTO 170                                                  
   50 CONTINUE                                                                  
  170 CONTINUE                                                                  
 1000 FORMAT(10E12.4)                                                           
 2000 FORMAT(10I12)                                                             
      STOP                                                                      
      END                                                                       
