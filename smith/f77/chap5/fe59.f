      PROGRAM P59                                                               
C                                                                               
C      PROGRAM 5.9 THREE-DIMENSIONAL ELASTIC                                    
C      ANALYSIS USING 8-NODE BRICK ELEMENTS                                     
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=1000,ILOADS=100,INF=50)                                     
C                                                                               
      REAL DEE(6,6),SAMP(4,2),COORD(8,3),JAC(3,3),JAC1(3,3),                    
     +DER(3,8),DERIV(3,8),BEE(6,24),DBEE(6,24),BTDB(24,24),KM(24,24),           
     +ELD(24),EPS(6),SIGMA(6),BT(24,6),FUN(8),KV(IKV),LOADS(ILOADS)             
      INTEGER G(24),NF(INF,3)                                                   
      DATA ISAMP/4/,IBTDB,IKM,IBT,IDOF/4*24/,IDEE,IBEE,IDBEE,IH/4*6/            
      DATA ICOORD,NOD/2*8/,IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*3/                 
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,NZE,N,IW,NN,NR,NGP,AA,BB,CC,E,V                          
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL FORMD3(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      DO 10 IS=1,NZE                                                            
      CALL GEO83D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)              
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      DO 20 K=1,NGP                                                             
      CALL FMLIN3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)                                
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
      CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   20 CONTINUE                                                                  
      CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
   10 CONTINUE                                                                  
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
C                                                                               
C      RECOVER ELEMENT STRAINS AND STRESSES                                     
C      AT BRICK CENTROIDS                                                       
C                                                                               
      NGP=1                                                                     
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 30 IP=1,NXE                                                            
      DO 30 IQ=1,NYE                                                            
      DO 30 IS=1,NZE                                                            
      CALL GEO83D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)              
      DO 40 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.0                                                   
   40 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      DO 30 I=1,NGP                                                             
      DO 30 J=1,NGP                                                             
      DO 30 K=1,NGP                                                             
      CALL FMLIN3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)                                
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TREEX3(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)                                    
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL PRINTV(SIGMA,IH)                                                     
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
