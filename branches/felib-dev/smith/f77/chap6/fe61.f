      PROGRAM P61                                                               
C                                                                               
C      PROGRAM 6.1 PLANE STRAIN OF AN ELASTIC-PLASTIC                           
C      (MOHR-COULOMB) SOLID USING 8-NODE QUADRILATERAL ELEMENTS                 
C      VISCOPLASTIC STRAIN METHOD (GRAVITY LOADING)                             
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=184,IKB2=38,ILOADS=184,INF=96,                             
     +          IEVPT=512,INX=20,INY=20,IFOS=10)                                
C                                                                               
      REAL DEE(4,4),SAMP(4,2),COORD(8,2),JAC(2,2),JAC1(2,2),DER(2,8),           
     +DERIV(2,8),BEE(4,16),DBEE(4,16),DEPTH(INY),GRAVLO(ILOADS),                
     +BTDB(16,16),KM(16,16),ELD(16),EPS(4),SIGMA(4),BLOAD(16),                  
     +BT(16,4),FUN(8),KB(IKB1,IKB2),LOADS(ILOADS),ELOAD(16),                    
     +BDYLDS(ILOADS),EVPT(IEVPT),OLDIS(ILOADS),                                 
     +ERATE(4),EVP(4),DEVP(4),M1(4,4),M2(4,4),M3(4,4),                          
     +FOS(IFOS),FLOW(4,4),STRESS(4),TOP(INX),BOT(INX)                           
      INTEGER NF(INF,2),G(16)                                                   
      DATA IDEE,IBEE,IDBEE,IH,IFLOW/5*4/,IDOF,IBTDB,IBT,IKM/4*16/               
      DATA IJAC,IJAC1,NODOF,IT,IDER,IDERIV/6*2/,ICOORD,NOD/2*8/                 
      DATA ISAMP/4/                                                             
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)PHI,C,PSI,GAMA,E,V,NXE,NYE,N,IW,NN,NR,NGP,ITS                    
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      READ(5,*)(TOP(I),I=1,NXE+1)                                               
      READ(5,*)(BOT(I),I=1,NXE+1)                                               
      READ(5,*)(DEPTH(I),I=1,NYE+1)                                             
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL NULVEC(OLDIS,N)                                                      
      CALL NULVEC(GRAVLO,N)                                                     
      CALL FMDRAD(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      PI=ACOS(-1.)                                                              
      TNPH=TAN(PHI*PI/180.)                                                     
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL SLOGEO(IP,IQ,NYE,TOP,BOT,DEPTH,COORD,ICOORD,G,NF,INF)                
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      CALL NULVEC(ELD,IDOF)                                                     
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
      DO 30 K=2,IDOF,2                                                          
   30 ELD(K)=ELD(K)+FUN(K/2)*QUOT                                               
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
   20 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
      CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
      DO 40 K=1,IDOF                                                            
   40 IF(G(K).NE.0)GRAVLO(G(K))=GRAVLO(G(K))-ELD(K)*GAMA                        
   10 CONTINUE                                                                  
C                                                                               
C      REDUCE EQUATIONS                                                         
C                                                                               
      CALL CHOLIN(KB,IKB1,N,IW)                                                 
C                                                                               
C      TRIAL FACTOR OF SAFETY LOOP                                              
C                                                                               
      READ(5,*)INCS,(FOS(I),I=1,INCS)                                           
      DO 50 IY=1,INCS                                                           
      PHIF=ATAN(TNPH/FOS(IY))*180./PI                                           
      SNPH=SIN(PHIF*PI/180.)                                                    
      DT=4.*(1.+V)*(1.-2.*V)/(E*(1.-2.*V+SNPH**2.))                             
      CF=C/FOS(IY)                                                              
      ITERS=0                                                                   
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(EVPT,NXE*NYE*IH*NGP*NGP)                                      
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
   60 ITERS=ITERS+1                                                             
      CALL VECADD(GRAVLO,BDYLDS,LOADS,N)                                        
      CALL CHOBAC(KB,IKB1,LOADS,N,IW)                                           
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(LOADS,OLDIS,N,0.0001,ICON)                                    
      IF(ITERS.EQ.1)ICON=0                                                      
      IF(ICON.EQ.1.OR.ITERS.EQ.ITS)CALL NULVEC(BDYLDS,N)                        
C                                                                               
C      INSPECT ALL GAUSS POINTS                                                 
C                                                                               
      NM=0                                                                      
      DO 70 IP=1,NXE                                                            
      DO 70 IQ=1,NYE                                                            
      NM=NM+1                                                                   
      CALL NULVEC(BLOAD,IDOF)                                                   
      CALL SLOGEO(IP,IQ,NYE,TOP,BOT,DEPTH,COORD,ICOORD,G,NF,INF)                
      DO 80 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   80 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      IG=0                                                                      
      DO 90 I=1,NGP                                                             
      DO 90 J=1,NGP                                                             
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
      DO 100 K=1,IH                                                             
  100 EPS(K)=EPS(K)-EVPT(IN+K)                                                  
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL INVAR(SIGMA,SIGM,DSBAR,THETA)                                        
C                                                                               
C      CHECK WHETHER YIELD IS VIOLATED                                          
C                                                                               
      CALL MOCOUF(PHIF,CF,SIGM,DSBAR,THETA,F)                                   
      IF(ICON.EQ.1.OR.ITERS.EQ.ITS)GOTO 110                                     
      IF(F.LT.0.)GOTO 90                                                        
      CALL MOCOUQ(PSI,DSBAR,THETA,DQ1,DQ2,DQ3)                                  
      CALL FORMM(SIGMA,M1,M2,M3)                                                
      DO 120 L=1,IH                                                             
      DO 120 M=1,IH                                                             
  120 FLOW(L,M)=F*(M1(L,M)*DQ1+M2(L,M)*DQ2+M3(L,M)*DQ3)                         
      CALL MVMULT(FLOW,IFLOW,SIGMA,IH,IH,ERATE)                                 
      DO 130 K=1,IH                                                             
      EVP(K)=ERATE(K)*DT                                                        
  130 EVPT(IN+K)=EVPT(IN+K)+EVP(K)                                              
      CALL MVMULT(DEE,IDEE,EVP,IH,IH,DEVP)                                      
      GOTO 140                                                                  
  110 CALL VECCOP(SIGMA,DEVP,IH)                                                
  140 CALL MVMULT(BT,IBT,DEVP,IDOF,IH,ELOAD)                                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      DO 150 K=1,IDOF                                                           
  150 BLOAD(K)=BLOAD(K)+ELOAD(K)*QUOT                                           
   90 CONTINUE                                                                  
C                                                                               
C      COMPUTE TOTAL BODYLOADS VECTOR                                           
C                                                                               
      DO 160 M=1,IDOF                                                           
      IF(G(M).EQ.0)GOTO 160                                                     
      BDYLDS(G(M))=BDYLDS(G(M))+BLOAD(M)                                        
  160 CONTINUE                                                                  
   70 CONTINUE                                                                  
      IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 60                                     
      BIG=0.                                                                    
      DO 170 I=1,N                                                              
  170 IF(ABS(LOADS(I)).GT.BIG)BIG=ABS(LOADS(I))                                 
      WRITE(6,1000)FOS(IY),BIG                                                  
      WRITE(6,2000)ITERS                                                        
      IF(ITERS.EQ.ITS)GOTO 180                                                  
   50 CONTINUE                                                                  
  180 CONTINUE                                                                  
 1000 FORMAT(10E12.4)                                                           
 2000 FORMAT(10I12)                                                             
      STOP                                                                      
      END                                                                       
