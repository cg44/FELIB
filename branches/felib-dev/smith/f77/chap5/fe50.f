      PROGRAM P50                                                               
C                                                                               
C      PROGRAM 5.0 PLANE STRESS OF AN ELASTIC                                   
C      SOLID USING 3-NODE TRIANGULAR ELEMENTS                                   
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=2400,ILOADS=200,INF=100)                                    
C                                                                               
      REAL DEE(3,3),SAMP(16,2),COORD(3,2),JAC(2,2),JAC1(2,2),DER(2,3),          
     +DERIV(2,3),BEE(3,6),DBEE(3,6),BTDB(6,6),KM(6,6),ELD(6),                   
     +EPS(3),SIGMA(3),BT(6,3),FUN(3),WT(16),KV(IKV),LOADS(ILOADS)               
      INTEGER NF(INF,2),G(6)                                                    
      DATA ISAMP/16/,IBTDB,IKM,IBT,IDOF/4*6/                                    
      DATA IDEE,ICOORD,IBEE,IDBEE,IH,NOD/6*3/                                   
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/                                 
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NCE,NYE,N,IW,NN,NR,NIP,AA,BB,E,V                                 
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IR=(IW+1)*N                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL FMDSIG(DEE,IDEE,E,V)                                                 
      CALL NUMINT(SAMP,ISAMP,WT,NIP)                                            
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NCE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEOM3X(IP,IQ,NCE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 20 I=1,NIP                                                             
      CALL FMTRI3(DER,IDER,FUN,SAMP,ISAMP,I)                                    
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=.5*DET*WT(I)                                                         
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
   20 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   10 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
C                                                                               
C      RECOVER STRESSES AT TRIANGLE CENTRES                                     
C                                                                               
      DO 30 IP=1,NCE                                                            
      DO 30 IQ=1,NYE                                                            
      CALL GEOM3X(IP,IQ,NCE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      DO 40 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   40 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      DO 30 I=1,NIP                                                             
      CALL FMTRI3(DER,IDER,FUN,SAMP,ISAMP,I)                                    
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL PRINTV(SIGMA,IH)                                                     
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
