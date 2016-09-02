      PROGRAM P52                                                               
C                                                                               
C      PROGRAM 5.2 PLANE STRAIN OF AN ELASTIC                                   
C      SOLID USING 15-NODE TRIANGULAR ELEMENTS                                  
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=200,IKB2=42,ILOADS=200,INF=100,                            
     +INC=20,INY=20)                                                            
C                                                                               
      REAL DEE(3,3),SAMP(16,2),COORD(15,2),JAC(2,2),JAC1(2,2),DER(2,15),        
     +DERIV(2,15),BEE(3,30),DBEE(3,30),BTDB(30,30),KM(30,30),ELD(30),           
     +EPS(3),SIGMA(3),BT(30,3),FUN(15),WT(16),KB(IKB1,IKB2),                    
     +LOADS(ILOADS),WID(INC),DEP(INY)                                           
      INTEGER NF(INF,2),G(30)                                                   
      DATA ISAMP/16/,IBTDB,IKM,IBT,IDOF/4*30/                                   
      DATA IDEE,IBEE,IDBEE,IH/4*3/,ICOORD,NOD/2*15/                             
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/                                 
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NCE,NYE,N,IW,NN,NR,NIP,E,V                                       
      READ(5,*)(WID(I),I=1,NCE+1)                                               
      READ(5,*)(DEP(I),I=1,(NYE+2)/2)                                           
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL FMDEPS(DEE,IDEE,E,V)                                                 
      CALL NUMINT(SAMP,ISAMP,WT,NIP)                                            
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NCE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEO15Y(IP,IQ,NYE,WID,DEP,COORD,ICOORD,G,NF,INF)                      
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 20 I=1,NIP                                                             
      CALL FMTR15(DER,IDER,FUN,SAMP,ISAMP,I)                                    
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
   10 CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL CHOLIN(KB,IKB1,N,IW)                                                 
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL CHOBAC(KB,IKB1,LOADS,N,IW)                                           
      CALL PRINTV(LOADS,N)                                                      
C                                                                               
C      RECOVER STRESSES AT TRIANGLE CENTRES                                     
C                                                                               
      NIP=1                                                                     
      CALL NUMINT(SAMP,ISAMP,WT,NIP)                                            
      DO 30 IP=1,NCE                                                            
      DO 30 IQ=1,NYE                                                            
      CALL GEO15Y(IP,IQ,NYE,WID,DEP,COORD,ICOORD,G,NF,INF)                      
      DO 40 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   40 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      DO 30 I=1,NIP                                                             
      CALL FMTR15(DER,IDER,FUN,SAMP,ISAMP,I)                                    
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
