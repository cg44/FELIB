      PROGRAM P58                                                               
C                                                                               
C      PROGRAM 5.8 THREE-DIMENSIONAL ELASTIC                                    
C      ANALYSIS USING 4-NODE TETRAHEDRON ELEMENTS                               
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=320,IKB2=12,ILOADS=320,ISTOR=20)                           
C                                                                               
      REAL DEE(6,6),SAMP(16,3),COORD(4,3),JAC(3,3),JAC1(3,3),                   
     +DER(3,4),DERIV(3,4),BEE(6,12),DBEE(6,12),BTDB(12,12),KM(12,12),           
     +ELD(12),EPS(6),SIGMA(6),BT(12,6),FUN(4),WT(16),                           
     +KB(IKB1,IKB2),LOADS(ILOADS),STOREC(ISTOR,4,3)                             
      INTEGER G(12),STOREG(ISTOR,12)                                            
      DATA ISAMP/16/,IBTDB,IKM,IBT,IDOF/4*12/,IDEE,IBEE,IDBEE,IH/4*6/           
      DATA ICOORD,NOD/2*4/,IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*3/                 
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NEL,N,IW,NIP,E,V                                                 
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL FORMD3(DEE,IDEE,E,V)                                                 
      CALL NUMIN3(SAMP,ISAMP,WT,NIP)                                            
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NEL                                                            
      READ(5,*)((COORD(I,J),J=1,3),I=1,4)                                       
      DO 20 I=1,NOD                                                             
      DO 20 J=1,IT                                                              
   20 STOREC(IP,I,J)=COORD(I,J)                                                 
      READ(5,*)(G(I),I=1,IDOF)                                                  
      DO 30 I=1,IDOF                                                            
   30 STOREG(IP,I)=G(I)                                                         
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 40 I=1,NIP                                                             
      CALL FMTET4(DER,IDER,FUN,SAMP,ISAMP,I)                                    
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TREEX3(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)                                    
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=DET*WT(I)/6.                                                         
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
      CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   40 CONTINUE                                                                  
      CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
   10 CONTINUE                                                                  
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL CHOLIN(KB,IKB1,N,IW)                                                 
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL CHOBAC(KB,IKB1,LOADS,N,IW)                                           
      CALL PRINTV(LOADS,N)                                                      
C                                                                               
C      RECOVER ELEMENT STRAINS AND STRESSES                                     
C      AT TETRAHEDRON CENTROIDS                                                 
C                                                                               
      DO 50 IP=1,NEL                                                            
      DO 60 I=1,NOD                                                             
      DO 60 J=1,IT                                                              
   60 COORD(I,J)=STOREC(IP,I,J)                                                 
      DO 70 I=1,IDOF                                                            
   70 G(I)=STOREG(IP,I)                                                         
      DO 50 I=1,NIP                                                             
      CALL FMTET4(DER,IDER,FUN,SAMP,ISAMP,I)                                    
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TREEX3(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)                                    
      DO 80 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.0                                                   
   80 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL PRINTV(SIGMA,IH)                                                     
   50 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
