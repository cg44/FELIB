      PROGRAM P510                                                              
C                                                                               
C      PROGRAM 5.10 THREE-DIMENSIONAL ELASTIC                                   
C      ANALYSIS USING 20-NODE BRICK ELEMENTS                                    
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=4388,ILOADS=124,INF=70)                                     
C                                                                               
      REAL DEE(6,6),SAMP(4,2),COORD(20,3),JAC(3,3),JAC1(3,3),                   
     +DER(3,20),DERIV(3,20),BEE(6,60),DBEE(6,60),BTDB(60,60),KM(60,60),         
     +ELD(60),EPS(6),SIGMA(6),BT(60,6),FUN(20),KV(IKV),LOADS(ILOADS)            
      INTEGER G(60),NF(INF,3),KDIAG(ILOADS)                                     
      DATA ISAMP/4/,IBTDB,IKM,IBT,IDOF/4*60/,IDEE,IBEE,IDBEE,IH/4*6/            
      DATA ICOORD,NOD/2*20/,IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*3/                
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,NZE,N,NN,NR,NGP,AA,BB,CC,E,V                             
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
      CALL FORMD3(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 40 IP=1,NXE                                                            
      DO 40 IQ=1,NYE                                                            
      DO 40 IS=1,NZE                                                            
      CALL GE203D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)              
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 50 I=1,NGP                                                             
      DO 50 J=1,NGP                                                             
      DO 50 K=1,NGP                                                             
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
      CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   50 CONTINUE                                                                  
      CALL FSPARV(KV,KM,IKM,G,KDIAG,IDOF)                                       
   40 CONTINUE                                                                  
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL SPARIN(KV,N,KDIAG)                                                   
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL SPABAC(KV,LOADS,N,KDIAG)                                             
      CALL PRINTV(LOADS,34)                                                     
C                                                                               
C      RECOVER ELEMENT STRAINS AND STRESSES                                     
C      AT BRICK CENTROIDS                                                       
C                                                                               
      NGP=1                                                                     
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 60 IP=1,NXE                                                            
      DO 60 IQ=1,NYE                                                            
      DO 60 IS=1,NZE                                                            
      CALL GE203D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)              
      DO 70 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.0                                                   
   70 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      DO 60 I=1,NGP                                                             
      DO 60 J=1,NGP                                                             
      DO 60 K=1,NGP                                                             
      CALL FMQUA3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)                                
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TREEX3(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)                                    
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL PRINTV(SIGMA,IH)                                                     
   60 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
