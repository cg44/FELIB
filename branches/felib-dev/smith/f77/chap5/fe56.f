      PROGRAM P56                                                               
C                                                                               
C      PROGRAM 5.6 AXISYMMETRIC STRAIN OF AN ELASTIC                            
C      SOLID USING 4-NODE QUADRILATERAL ELEMENTS                                
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=3500,ILOADS=100,INF=100,INR=20,IND=20,                      
     +          IPROP=20)                                                       
C                                                                               
      REAL DEE(4,4),SAMP(4,2),COORD(4,2),JAC(2,2),JAC1(2,2),                    
     +DER(2,4),DERIV(2,4),BEE(4,8),DBEE(4,8),KV(IKV),LOADS(ILOADS),             
     +BTDB(8,8),KM(8,8),ELD(8),EPS(4),SIGMA(4),                                 
     +BT(8,4),FUN(4),RAD(INR),DEP(IND),PROP(IPROP,2)                            
      INTEGER G(8),NF(INF,2)                                                    
      DATA IDEE,IBEE,IDBEE,IH/4*4/,IDOF,IBTDB,IBT,IKM/4*8/                      
      DATA IJAC,IJAC1,NODOF,IT,IDER,IDERIV/6*2/,ICOORD,NOD/2*4/                 
      DATA ISAMP/4/                                                             
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NRE,NDE,N,IW,NN,NR,NGP                                           
      READ(5,*)(RAD(I),I=1,NRE+1)                                               
      READ(5,*)(DEP(I),I=1,NDE+1)                                               
      READ(5,*)((PROP(I,J),I=1,NDE),J=1,2)                                      
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NRE                                                            
      DO 10 IQ=1,NDE                                                            
      E=PROP(IQ,1)                                                              
      V=PROP(IQ,2)                                                              
      CALL FMDRAD(DEE,IDEE,E,V)                                                 
      CALL GEOV4Y(IP,IQ,NDE,RAD,DEP,COORD,ICOORD,G,NF,INF)                      
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
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
C      RECOVER STRESSES AT ELEMENT CENTRES                                      
C                                                                               
      NGP=1                                                                     
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 30 IP=1,NRE                                                            
      DO 30 IQ=1,NDE                                                            
      E=PROP(IQ,1)                                                              
      V=PROP(IQ,2)                                                              
      CALL FMDRAD(DEE,IDEE,E,V)                                                 
      CALL GEOV4Y(IP,IQ,NDE,RAD,DEP,COORD,ICOORD,G,NF,INF)                      
      DO 40 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   40 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      DO 30 I=1,NGP                                                             
      DO 30 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FMBRAD(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)               
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL PRINTV(EPS,IH)                                                       
      CALL PRINTV(SIGMA,IH)                                                     
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
