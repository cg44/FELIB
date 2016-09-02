      PROGRAM P60                                                               
C                                                                               
C      PROGRAM 6.0 PLANE STRAIN OF AN ELASTIC-PLASTIC                           
C      (VON-MISES) SOLID USING 8-NODE QUADRILATERAL ELEMENTS                    
C      VISCOPLASTIC STRAIN METHOD                                               
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=184,IKB2=30,ILOADS=184,INF=130,                            
     +          IEVPT=512,INX=20,INY=20,INO=10,IQINC=20)                        
C                                                                               
      REAL DEE(4,4),SAMP(4,2),COORD(8,2),JAC(2,2),JAC1(2,2),                    
     +DER(2,8),DERIV(2,8),BEE(4,16),DBEE(4,16),WIDTH(INX),DEPTH(INY),           
     +BTDB(16,16),KM(16,16),ELD(16),EPS(4),SIGMA(4),BLOAD(16),                  
     +BT(16,4),FUN(8),KB(IKB1,IKB2),LOADS(ILOADS),ELOAD(16),                    
     +TOTD(ILOADS),BDYLDS(ILOADS),EVPT(IEVPT),OLDIS(ILOADS),                    
     +SX(INX,INY,4),SY(INX,INY,4),TXY(INX,INY,4),SZ(INX,INY,4),                 
     +VAL(INO),ERATE(4),EVP(4),DEVP(4),M1(4,4),M2(4,4),M3(4,4),                 
     +FLOW(4,4),STRESS(4),QINC(IQINC)                                           
      INTEGER NF(INF,2),G(16),NO(INO)                                           
      DATA IDEE,IBEE,IDBEE,IH,IFLOW/5*4/,IDOF,IBTDB,IBT,IKM/4*16/               
      DATA IJAC,IJAC1,NODOF,IT,IDER,IDERIV/6*2/,ICOORD,NOD/2*8/                 
      DATA ISAMP/4/                                                             
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)CU,E,V,NXE,NYE,N,IW,NN,NR,NGP,ITS                                
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      READ(5,*)(WIDTH(I),I=1,NXE+1)                                             
      READ(5,*)(DEPTH(I),I=1,NYE+1)                                             
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL NULVEC(OLDIS,N)                                                      
      CALL NULVEC(TOTD,N)                                                       
      CALL FMDRAD(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DT=4.*(1.+V)/(3.*E)                                                       
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEOV8Y(IP,IQ,NYE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)                  
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      IG=0                                                                      
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      IG=IG+1                                                                   
      SX(IP,IQ,IG)=0.                                                           
      SY(IP,IQ,IG)=0.                                                           
      TXY(IP,IQ,IG)=0.                                                          
      SZ(IP,IQ,IG)=0.                                                           
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
   10 CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
C                                                                               
C      READ LOAD WEIGHTINGS AND REDUCE EQUATIONS                                
C                                                                               
      READ(5,*)NL,(NO(I),VAL(I),I=1,NL)                                         
      CALL CHOLIN(KB,IKB1,N,IW)                                                 
C                                                                               
C      LOAD INCREMENT LOOP                                                      
C                                                                               
      READ(5,*)INCS,(QINC(I),I=1,INCS)                                          
      PTOT=0.                                                                   
      DO 30 IY=1,INCS                                                           
      PTOT=PTOT+QINC(IY)                                                        
      ITERS=0                                                                   
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(EVPT,NXE*NYE*IH*NGP*NGP)                                      
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
   40 ITERS=ITERS+1                                                             
      CALL NULVEC(LOADS,N)                                                      
      DO 50 I=1,NL                                                              
   50 LOADS(NO(I))=VAL(I)*QINC(IY)                                              
      CALL VECADD(LOADS,BDYLDS,LOADS,N)                                         
      CALL CHOBAC(KB,IKB1,LOADS,N,IW)                                           
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(LOADS,OLDIS,N,0.001,ICON)                                     
      IF(ITERS.EQ.1)ICON=0                                                      
      IF(ICON.EQ.1.OR.ITERS.EQ.ITS)CALL NULVEC(BDYLDS,N)                        
C                                                                               
C      INSPECT ALL GAUSS POINTS                                                 
C                                                                               
      NM=0                                                                      
      DO 60 IP=1,NXE                                                            
      DO 60 IQ=1,NYE                                                            
      NM=NM+1                                                                   
      CALL NULVEC(BLOAD,IDOF)                                                   
      CALL GEOV8Y(IP,IQ,NYE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)                  
      DO 70 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   70 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      IG=0                                                                      
      DO 80 I=1,NGP                                                             
      DO 80 J=1,NGP                                                             
      IG=IG+1                                                                   
      IN=NGP*NGP*IH*(NM-1)+IH*(IG-1)                                            
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      DO 90 K=1,IH                                                              
   90 EPS(K)=EPS(K)-EVPT(IN+K)                                                  
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      STRESS(1)=SIGMA(1)+SX(IP,IQ,IG)                                           
      STRESS(2)=SIGMA(2)+SY(IP,IQ,IG)                                           
      STRESS(3)=SIGMA(3)+TXY(IP,IQ,IG)                                          
      STRESS(4)=SIGMA(4)+SZ(IP,IQ,IG)                                           
      CALL INVAR(STRESS,SIGM,DSBAR,THETA)                                       
C                                                                               
C      CHECK WHETHER YIELD IS VIOLATED                                          
C                                                                               
      F=DSBAR-SQRT(3.)*CU                                                       
      IF(ICON.EQ.1.OR.ITERS.EQ.ITS)GOTO 100                                     
      IF(F.LT.0.)GOTO 110                                                       
      DQ1=0.                                                                    
      DQ2=1.5/DSBAR                                                             
      DQ3=0.                                                                    
      CALL FORMM(STRESS,M1,M2,M3)                                               
      DO 120 L=1,IH                                                             
      DO 120 M=1,IH                                                             
  120 FLOW(L,M)=F*(M1(L,M)*DQ1+M2(L,M)*DQ2+M3(L,M)*DQ3)                         
      CALL MVMULT(FLOW,IFLOW,STRESS,IH,IH,ERATE)                                
      DO 130 K=1,IH                                                             
      EVP(K)=ERATE(K)*DT                                                        
  130 EVPT(IN+K)=EVPT(IN+K)+EVP(K)                                              
      CALL MVMULT(DEE,IDEE,EVP,IH,IH,DEVP)                                      
      GOTO 140                                                                  
  100 CALL VECCOP(STRESS,DEVP,IH)                                               
  140 CALL MVMULT(BT,IBT,DEVP,IDOF,IH,ELOAD)                                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      DO 150 K=1,IDOF                                                           
  150 BLOAD(K)=BLOAD(K)+ELOAD(K)*QUOT                                           
  110 IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 80                                     
C                                                                               
C      UPDATE GAUSS POINT STRESSES                                              
C                                                                               
      SX(IP,IQ,IG)=STRESS(1)                                                    
      SY(IP,IQ,IG)=STRESS(2)                                                    
      TXY(IP,IQ,IG)=STRESS(3)                                                   
      SZ(IP,IQ,IG)=STRESS(4)                                                    
   80 CONTINUE                                                                  
C                                                                               
C      COMPUTE TOTAL BODYLOADS VECTOR                                           
C                                                                               
      DO 160 M=1,IDOF                                                           
      IF(G(M).EQ.0)GOTO 160                                                     
      BDYLDS(G(M))=BDYLDS(G(M))+BLOAD(M)                                        
  160 CONTINUE                                                                  
   60 CONTINUE                                                                  
      IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 40                                     
      CALL VECADD(TOTD,LOADS,TOTD,N)                                            
      WRITE(6,1000)PTOT                                                         
      WRITE(6,1000)(TOTD(NO(I)),I=1,NL)                                         
      WRITE(6,2000)ITERS                                                        
      IF(ITERS.EQ.ITS)GOTO 170                                                  
   30 CONTINUE                                                                  
  170 CONTINUE                                                                  
 1000 FORMAT(10E12.4)                                                           
 2000 FORMAT(10I12)                                                             
      STOP                                                                      
      END                                                                       
