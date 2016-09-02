      PROGRAM P51                                                               
C                                                                               
C      PROGRAM 5.1 PLANE STRESS OF AN ELASTIC                                   
C      SOLID USING 6-NODE TRIANGULAR ELEMENTS                                   
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=200,IKB2=22,ILOADS=200,INF=100)                            
C                                                                               
      REAL DEE(3,3),SAMP(16,2),COORD(6,2),JAC(2,2),JAC1(2,2),DER(2,6),          
     +DERIV(2,6),BEE(3,12),DBEE(3,12),BTDB(12,12),KM(12,12),ELD(12),            
     +EPS(3),SIGMA(3),BT(12,3),FUN(6),WT(16),KB(IKB1,IKB2),LOADS(ILOADS)        
      INTEGER NF(INF,2),G(12)                                                   
      DATA ISAMP/16/,IBTDB,IKM,IBT,IDOF/4*12/                                   
      DATA IDEE,IBEE,IDBEE,IH/4*3/,ICOORD,NOD/2*6/                              
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/                                 
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NCE,NYE,N,IW,NN,NR,NIP,AA,BB,E,V                                 
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL FMDSIG(DEE,IDEE,E,V)                                                 
      CALL NUMINT(SAMP,ISAMP,WT,NIP)                                            
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NCE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEOM6X(IP,IQ,NCE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 20 I=1,NIP                                                             
      CALL FMTRI6(DER,IDER,FUN,SAMP,ISAMP,I)                                    
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
      WRITE(6,'(10E12.4)')(LOADS(I),I=1,81,20)                                  
C                                                                               
C      RECOVER STRESSES AT TRIANGLE CENTRES                                     
C                                                                               
      NIP=1                                                                     
      CALL NUMINT(SAMP,ISAMP,WT,NIP)                                            
      DO 30 IP=1,NCE                                                            
      DO 30 IQ=1,NYE                                                            
      CALL GEOM6X(IP,IQ,NCE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      DO 40 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   40 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      DO 30 I=1,NIP                                                             
      CALL FMTRI6(DER,IDER,FUN,SAMP,ISAMP,I)                                    
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      IF(IP.EQ.1.AND.IQ/2*2.NE.IQ)WRITE(6,'(E12.4)')SIGMA(2)                    
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
