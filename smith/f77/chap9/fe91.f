      PROGRAM P91                                                               
C                                                                               
C       PROGRAM 9.1 BIOT CONSOLIDATION OF AN ELASTIC SOLID IN PLANE             
C       STRAIN USING 4-NODE QUADRILATERALS                                      
C                                                                               
C                                                                               
C       ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                  
C                                                                               
      PARAMETER(IKV=5000,IPB2=50,ILOADS=200,INF=100,IWID=30,                    
     +          IDEP=30)                                                        
C                                                                               
      REAL DEE(3,3),SAMP(3,2),COORD(4,2),DERIVT(4,2),JAC(2,2),                  
     +JAC1(2,2),KAY(2,2),DER(2,4),DERIV(2,4),KDERIV(2,4),                       
     +BEE(3,8),DBEE(3,8),BT(8,3),BTDB(8,8),KM(8,8),ELD(8),                      
     +EPS(3),SIGMA(3),DTKD(4,4),KP(4,4),KE(12,12),KD(12,12),                    
     +FUN(4),C(8,4),VOLF(8,4),WIDTH(IWID),DEPTH(IDEP),KV(IKV),VOL(8),           
     +PB(ILOADS,IPB2),LOADS(ILOADS),ANS(ILOADS)                                 
      INTEGER G(12),NF(INF,3)                                                   
      DATA IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV,IT/7*2/                           
      DATA IDEE,ISAMP,IBEE,IDBEE,NODOF,IH/6*3/                                  
      DATA ICOORD,IDERVT,IDTKD,IKP,NOD/5*4/                                     
      DATA IBT,IBTDB,IKM,IC,IVOLF,IDOF/6*8/,IKE,IKD,ITOT/3*12/                  
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,PERMX,PERMY,E,V,DTIM,ISTEP,THETA          
      READ(5,*)(WIDTH(I),I=1,NXE+1)                                             
      READ(5,*)(DEPTH(I),I=1,NYE+1)                                             
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IBAND=2*(IW+1)-1                                                          
      IR=N*(IW+1)                                                               
      CALL NULL(PB,ILOADS,N,IBAND)                                              
      CALL NULVEC(LOADS,N)                                                      
      CALL NULVEC(KV,IR)                                                        
      CALL NULL(DEE,IDEE,IH,IH)                                                 
      CALL FMDEPS(DEE,IDEE,E,V)                                                 
      CALL NULL(KAY,IKAY,IT,IT)                                                 
      KAY(1,1)=PERMX                                                            
      KAY(2,2)=PERMY                                                            
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT MATRIX INTEGRATION AND ASSEMBLY                                  
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEV4X3(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)                  
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      CALL NULL(C,IC,IDOF,NOD)                                                  
      CALL NULL(KP,IKP,NOD,NOD)                                                 
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL VOL2D(BEE,IBEE,VOL,NOD)                                              
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
      CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
      CALL MATMUL(KAY,IKAY,DERIV,IDERIV,KDERIV,IKDERV,IT,IT,NOD)                
      CALL MATRAN(DERIVT,IDERVT,DERIV,IDERIV,IT,NOD)                            
      CALL MATMUL(DERIVT,IDERVT,KDERIV,IKDERV,DTKD,IDTKD,NOD,IT,NOD)            
      PROD=QUOT*DTIM                                                            
      CALL MSMULT(DTKD,IDTKD,PROD,NOD,NOD)                                      
      CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
      DO 30 K=1,IDOF                                                            
      DO 30 L=1,NOD                                                             
   30 VOLF(K,L)=VOL(K)*FUN(L)*QUOT                                              
      CALL MATADD(C,IC,VOLF,IVOLF,IDOF,NOD)                                     
   20 CONTINUE                                                                  
      CALL FMKDKE(KM,IKM,KP,IKP,C,IC,KE,IKE,KD,IKD,IDOF,NOD,ITOT,THETA)         
      CALL FORMKV(KV,KE,IKE,G,N,ITOT)                                           
      CALL FORMTB(PB,ILOADS,KD,IKD,G,IW,ITOT)                                   
   10 CONTINUE                                                                  
C                                                                               
C      REDUCE LEFT HAND SIDE                                                    
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
C                                                                               
C      TIME STEPPING LOOP                                                       
C                                                                               
      DO 40 NS=1,ISTEP                                                          
      CALL BANTML(PB,ILOADS,LOADS,ANS,N,IW)                                     
C                                                                               
C      RAMP LOADING                                                             
C                                                                               
      X1=(.1*NS+.1*(THETA-1.))*.5                                               
      X2=X1*2.                                                                  
      IF(NS.GT.10)ANS(1)=ANS(1)-.5                                              
      IF(NS.LE.10)ANS(1)=ANS(1)-X1                                              
      IF(NS.GT.10)ANS(3)=ANS(3)-1.                                              
      IF(NS.LE.10)ANS(3)=ANS(3)-X2                                              
      IF(NS.GT.10)ANS(5)=ANS(5)-1.                                              
      IF(NS.LE.10)ANS(5)=ANS(5)-X2                                              
      IF(NS.GT.10)ANS(6)=ANS(6)-.5                                              
      IF(NS.LE.10)ANS(6)=ANS(6)-X1                                              
      CALL BACSUB(KV,ANS,N,IW)                                                  
      CALL VECCOP(ANS,LOADS,N)                                                  
      CALL PRINTV(ANS,N)                                                        
C                                                                               
C      RECOVER ELEMENT EFFECTIVE STRESSES AT ELEMENT 'CENTRES'                  
C                                                                               
      NGP=1                                                                     
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 50 IP=1,NXE                                                            
      DO 50 IQ=1,NYE                                                            
      CALL GEV4X3(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)                  
      DO 60 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   60 IF(G(M).NE.0)ELD(M)=ANS(G(M))                                             
      DO 70 I=1,NGP                                                             
      DO 70 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      CALL PRINTV(SIGMA,IH)                                                     
   70 CONTINUE                                                                  
   50 CONTINUE                                                                  
   40 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
