      PROGRAM P57                                                               
C                                                                               
C      PROGRAM 5.7 NON-AXISYMMETRIC STRAIN OF AXISYMMETRIC                      
C      SOLIDS USING 8-NODE QUADRILATERAL ELEMENTS                               
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=100,IKB2=35,ILOADS=100,INF=100)                            
C                                                                               
      REAL DEE(6,6),SAMP(4,2),COORD(8,2),JAC(2,2),JAC1(2,2),                    
     +DER(2,8),DERIV(2,8),BEE(6,24),DBEE(6,24),                                 
     +BTDB(24,24),KM(24,24),ELD(24),EPS(6),SIGMA(6),                            
     +BT(24,6),FUN(8),KB(IKB1,IKB2),LOADS(ILOADS)                               
      INTEGER G(24),NF(INF,3)                                                   
      DATA IDEE,IBEE,IDBEE,IH/4*6/,IDOF,IBTDB,IBT,IKM/4*24/                     
      DATA IJAC,IJAC1,IT,IDER,IDERIV/5*2/,ICOORD,NOD/2*8/                       
      DATA ISAMP/4/,NODOF/3/                                                    
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NRE,NDE,N,IW,NN,NR,NGP,AA,BB,E,V                                 
      READ(5,*)LTH,IFLAG,CHI                                                    
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      PI=4.*ATAN(1.)                                                            
      CHI=CHI*PI/180.                                                           
      CA=COS(CHI)                                                               
      SA=SIN(CHI)                                                               
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL FORMD3(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NRE                                                            
      DO 10 IQ=1,NDE                                                            
      CALL GENA8X(IP,IQ,NRE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL BNONAX(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,                       
     +            SUM,NOD,IFLAG,LTH)                                            
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=SUM*DET*SAMP(I,2)*SAMP(J,2)                                          
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
C      RECOVER STRAINS AND STRESSES AT ELEMENT CENTRES                          
C                                                                               
      NGP=1                                                                     
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 30 IP=1,NRE                                                            
      DO 30 IQ=1,NDE                                                            
      CALL GENA8X(IP,IQ,NRE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      DO 40 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   40 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      DO 30 I=1,NGP                                                             
      DO 30 J=1,NGP                                                             
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL BNONAX(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,                       
     +            SUM,NOD,IFLAG,LTH)                                            
      DO 50 L=1,IDOF                                                            
      DO 60 K=1,4                                                               
   60 BEE(K,L)=BEE(K,L)*CA                                                      
      DO 50 K=5,6                                                               
   50 BEE(K,L)=BEE(K,L)*SA                                                      
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL PRINTV(SIGMA,IH)                                                     
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
