      PROGRAM P64                                                               
C                                                                               
C      PROGRAM 6.4 THREE-DIMENSIONAL ELASTO-PLASTIC                             
C      ANALYSIS USING 20-NODE BRICK ELEMENTS                                    
C      MOHR-COULOMB'S FAILURE CRITERION                                         
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=1000,ILOADS=200,INF=100,IEL=10,                             
     +          IEVPT=500,INO=20)                                               
C                                                                               
      REAL DEE(6,6),SAMP(3,2),COORD(20,3),JAC(3,3),JAC1(3,3),                   
     +DER(3,20),DERIV(3,20),BEE(6,60),DBEE(6,60),BTDB(60,60),KM(60,60),         
     +ELD(60),EPS(6),SIGMA(6),BT(60,6),FUN(20),KV(IKV),LOADS(ILOADS),           
     +SX(IEL,27),SY(IEL,27),SZ(IEL,27),TXY(IEL,27),TYZ(IEL,27),                 
     +TZX(IEL,27),BLOAD(60),ELOAD(60),TOTD(ILOADS),OLDIS(ILOADS),               
     +VAL(INO),ERATE(6),EVP(6),EVPT(IEVPT),BDYLDS(ILOADS),DEVP(6),              
     +M1(6,6),M2(6,6),M3(6,6),FLOW(6,6),STRESS(6),STORKV(INO)                   
      INTEGER NO(INO),G(60),NF(INF,3),KDIAG(ILOADS)                             
      DATA IBTDB,IKM,IBT,IDOF/4*60/,IDEE,IBEE,IDBEE,IH,IFLOW/5*6/               
      DATA ICOORD,NOD/2*20/,IJAC,IJAC1,IDER,IDERIV,NODOF,IT,ISAMP/7*3/          
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)PHI,C,PSI,E,V,CONS,NXE,NYE,NZE,N,NN,NR,NGP,AA,BB,CC              
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      DO 10 I=1,N                                                               
   10 KDIAG(I)=0                                                                
      DO 20 IP=1,NXE                                                            
      DO 20 IQ=1,NYE                                                            
      DO 20 IS=1,NZE                                                            
      CALL GE203D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)              
   20 CALL FKDIAG(KDIAG,G,IDOF)                                                 
      KDIAG(1)=1                                                                
      DO 30 I=2,N                                                               
   30 KDIAG(I)=KDIAG(I)+KDIAG(I-1)                                              
      IR=KDIAG(N)                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL NULVEC(TOTD,N)                                                       
      CALL NULVEC(OLDIS,N)                                                      
      CALL FORMD3(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      SNPH=SIN(PHI*ACOS(-1.)/180.)                                              
      DT=4.*(1.+V)*(1.-2.*V)/(E*(1.-2.*V+SNPH*SNPH))                            
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      NM=0                                                                      
      DO 40 IP=1,NXE                                                            
      DO 40 IQ=1,NYE                                                            
      DO 40 IS=1,NZE                                                            
      NM=NM+1                                                                   
      CALL GE203D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)              
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      IG=0                                                                      
      DO 50 I=1,NGP                                                             
      DO 50 J=1,NGP                                                             
      DO 50 K=1,NGP                                                             
      IG=IG+1                                                                   
      SX(NM,IG)=CONS                                                            
      SY(NM,IG)=CONS                                                            
      SZ(NM,IG)=CONS                                                            
      TXY(NM,IG)=0.                                                             
      TYZ(NM,IG)=0.                                                             
      TZX(NM,IG)=0.                                                             
      CALL FMQUA3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)                                
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TREEX3(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)                                    
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)*SAMP(K,2)                                    
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
   50 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   40 CALL FSPARV(KV,KM,IKM,G,KDIAG,IDOF)                                       
C                                                                               
C      READ PRESCRIBED FREEDOMS                                                 
C      AUGMENT AND REDUCE STIFFNESS MATRIX                                      
C                                                                               
      READ(5,*)NL,(NO(I),I=1,NL),PRESC,INCS,ITS                                 
      DO 60 I=1,NL                                                              
      KV(KDIAG(NO(I)))=KV(KDIAG(NO(I)))+1.E20                                   
   60 STORKV(I)=KV(KDIAG(NO(I)))                                                
      CALL SPARIN(KV,N,KDIAG)                                                   
C                                                                               
C      DISPLACEMENT INCREMENT LOOP                                              
C                                                                               
      DO 70 IY=1,INCS                                                           
      PTOT=PRESC*IY                                                             
      ITERS=0                                                                   
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(EVPT,NXE*NYE*NZE*IH*NGP*NGP*NGP)                              
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
   80 ITERS=ITERS+1                                                             
      CALL NULVEC(LOADS,N)                                                      
      DO 90 I=1,NL                                                              
   90 LOADS(NO(I))=STORKV(I)*PRESC                                              
      CALL VECADD(LOADS,BDYLDS,LOADS,N)                                         
      CALL SPABAC(KV,LOADS,N,KDIAG)                                             
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(LOADS,OLDIS,N,.0001,ICON)                                     
      IF(ITERS.EQ.1)ICON=0                                                      
      IF(ICON.EQ.1.OR.ITERS.EQ.ITS)CALL NULVEC(BDYLDS,N)                        
C                                                                               
C      INSPECT ALL GAUSS POINTS                                                 
C                                                                               
      NM=0                                                                      
      DO 100 IP=1,NXE                                                           
      DO 100 IQ=1,NYE                                                           
      DO 100 IS=1,NZE                                                           
      NM=NM+1                                                                   
      CALL NULVEC(BLOAD,IDOF)                                                   
      CALL GE203D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)              
      DO 110 M=1,IDOF                                                           
      IF(G(M).EQ.0)ELD(M)=0.                                                    
  110 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      IG=0                                                                      
      DO 120 I=1,NGP                                                            
      DO 120 J=1,NGP                                                            
      DO 120 K=1,NGP                                                            
      IG=IG+1                                                                   
      IN=NGP*NGP*NGP*IH*(NM-1)+IH*(IG-1)                                        
      CALL FMQUA3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)                                
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TREEX3(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)                                    
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      DO 130 L=1,IH                                                             
  130 EPS(L)=EPS(L)-EVPT(IN+L)                                                  
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      STRESS(1)=SIGMA(1)+SX(NM,IG)                                              
      STRESS(2)=SIGMA(2)+SY(NM,IG)                                              
      STRESS(3)=SIGMA(3)+SZ(NM,IG)                                              
      STRESS(4)=SIGMA(4)+TXY(NM,IG)                                             
      STRESS(5)=SIGMA(5)+TYZ(NM,IG)                                             
      STRESS(6)=SIGMA(6)+TZX(NM,IG)                                             
      CALL INVAR3(STRESS,SIGM,DSBAR,THETA)                                      
C                                                                               
C      CHECK WHETHER YIELD IS VIOLATED                                          
C                                                                               
      CALL MOCOUF(PHI,C,SIGM,DSBAR,THETA,F)                                     
      IF(ICON.EQ.1.OR.ITERS.EQ.ITS)GOTO 140                                     
      IF(F.LT.0.)GOTO 150                                                       
      CALL MOCOUQ(PSI,DSBAR,THETA,DQ1,DQ2,DQ3)                                  
      CALL FORMM3(STRESS,M1,M2,M3)                                              
      DO 160 L=1,IH                                                             
      DO 160 M=1,IH                                                             
  160 FLOW(L,M)=F*(M1(L,M)*DQ1+M2(L,M)*DQ2+M3(L,M)*DQ3)                         
      CALL MVMULT(FLOW,IFLOW,STRESS,IH,IH,ERATE)                                
      DO 170 L=1,IH                                                             
      EVP(L)=ERATE(L)*DT                                                        
  170 EVPT(IN+L)=EVPT(IN+L)+EVP(L)                                              
      CALL MVMULT(DEE,IDEE,EVP,IH,IH,DEVP)                                      
      GOTO 180                                                                  
  140 CALL VECCOP(STRESS,DEVP,IH)                                               
  180 CALL MVMULT(BT,IBT,DEVP,IDOF,IH,ELOAD)                                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)*SAMP(K,2)                                    
      DO 190 L=1,IDOF                                                           
  190 BLOAD(L)=BLOAD(L)+ELOAD(L)*QUOT                                           
  150 IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 120                                    
C                                                                               
C      UPDATE GAUSS POINT STRESSES                                              
C                                                                               
      SX(NM,IG)=STRESS(1)                                                       
      SY(NM,IG)=STRESS(2)                                                       
      SZ(NM,IG)=STRESS(3)                                                       
      TXY(NM,IG)=STRESS(4)                                                      
      TYZ(NM,IG)=STRESS(5)                                                      
      TZX(NM,IG)=STRESS(6)                                                      
  120 CONTINUE                                                                  
C                                                                               
C      COMPUTE TOTAL BODYLOADS VECTOR                                           
C                                                                               
      DO 200 M=1,IDOF                                                           
      IF(G(M).EQ.0)GOTO 200                                                     
      BDYLDS(G(M))=BDYLDS(G(M))+BLOAD(M)                                        
  200 CONTINUE                                                                  
  100 CONTINUE                                                                  
      IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 80                                     
      CALL VECADD(TOTD,LOADS,TOTD,N)                                            
      WRITE(6,1000)PTOT                                                         
      WRITE(6,1000)SZ(1,1),SX(1,1),SY(1,1)                                      
      WRITE(6,2000)ITERS                                                        
      IF(ITERS.EQ.ITS)GOTO 210                                                  
   70 CONTINUE                                                                  
  210 CONTINUE                                                                  
 1000 FORMAT(10E12.4)                                                           
 2000 FORMAT(10I12)                                                             
      STOP                                                                      
      END                                                                       
