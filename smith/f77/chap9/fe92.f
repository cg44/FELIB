      PROGRAM P92                                                               
C                                                                               
C       PROGRAM 9.2 BIOT CONSOLIDATION OF AN ELASTIC SOLID IN PLANE             
C       STRAIN USING 4-NODE QUADRILATERALS FOR FLUID PHASE AND                  
C       8-NODE QUADRILATERALS FOR THE SOLID PHASE                               
C                                                                               
C                                                                               
C       ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                  
C                                                                               
      PARAMETER(IKV=5000,IPB2=50,ILOADS=200,INF=100)                            
C                                                                               
      REAL DEE(3,3),SAMP(3,2),COORD(8,2),DERIVT(8,2),JAC(2,2),                  
     +JAC1(2,2),KAY(2,2),DER(2,8),DERIV(2,8),KDERIV(2,8),                       
     +BEE(3,16),DBEE(3,16),BT(16,3),BTDB(16,16),KM(16,16),ELD(16),              
     +EPS(3),SIGMA(3),DTKD(4,4),KP(4,4),KE(20,20),KD(20,20),                    
     +FUN(8),C(16,4),VOLF(16,4),WIDTH(30),DEPTH(30),KV(IKV),VOL(16),            
     +FUNF(4),COORDF(4,2),DERF(2,4),DERIVF(2,4),                                
     +PB(ILOADS,IPB2),LOADS(ILOADS),ANS(ILOADS)                                 
      INTEGER G(20),NF(INF,3)                                                   
      DATA IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV,IT,IDERF,IDERVF/9*2/              
      DATA IDEE,ISAMP,IBEE,IDBEE,NODOF,IH/6*3/                                  
      DATA ICORDF,IDTKD,IKP,NODF/4*4/                                           
      DATA ICOORD,IDERVT,NOD/3*8/,IBT,IBTDB,IKM,IC,IVOLF,IDOF/6*16/             
      DATA IKE,IKD,ITOT/3*20/                                                   
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
      CALL GEOUVP(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,COORDF,ICORDF,             
     +            G,NF,INF)                                                     
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      CALL NULL(C,IC,IDOF,NODF)                                                 
      CALL NULL(KP,IKP,NODF,NODF)                                               
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
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
C                                                                               
C      FLUID CONTRIBUTION                                                       
C                                                                               
      CALL FORMLN(DERF,IDERF,FUNF,SAMP,ISAMP,I,J)                               
      CALL MATMUL(DERF,IDERF,COORDF,ICORDF,JAC,IJAC,IT,NODF,IT)                 
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DERF,IDERF,DERIVF,IDERVF,IT,IT,NODF)               
      CALL MATMUL(KAY,IKAY,DERIVF,IDERVF,KDERIV,IKDERV,IT,IT,NODF)              
      CALL MATRAN(DERIVT,IDERVT,DERIVF,IDERVF,IT,NODF)                          
      CALL MATMUL(DERIVT,IDERVT,KDERIV,IKDERV,DTKD,IDTKD,NODF,IT,NODF)          
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      PROD=QUOT*DTIM                                                            
      CALL MSMULT(DTKD,IDTKD,PROD,NODF,NODF)                                    
      CALL MATADD(KP,IKP,DTKD,IDTKD,NODF,NODF)                                  
      DO 30 K=1,IDOF                                                            
      DO 30 L=1,NODF                                                            
   30 VOLF(K,L)=VOL(K)*FUNF(L)*QUOT                                             
      CALL MATADD(C,IC,VOLF,IVOLF,IDOF,NODF)                                    
   20 CONTINUE                                                                  
      CALL FMKDKE(KM,IKM,KP,IKP,C,IC,KE,IKE,KD,IKD,IDOF,NODF,ITOT,THETA)        
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
      X1=(.1*NS+.1*(THETA-1.))/6.                                               
      X2=X1*4.                                                                  
      IF(NS.GT.10)ANS(1)=ANS(1)-1./6.                                           
      IF(NS.LE.10)ANS(1)=ANS(1)-X1                                              
      IF(NS.GT.10)ANS(3)=ANS(3)-2./3.                                           
      IF(NS.LE.10)ANS(3)=ANS(3)-X2                                              
      IF(NS.GT.10)ANS(4)=ANS(4)-1./6.                                           
      IF(NS.LE.10)ANS(4)=ANS(4)-X1                                              
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
      CALL GEOUVP(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,COORDF,ICORDF,             
     +            G,NF,INF)                                                     
      DO 60 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   60 IF(G(M).NE.0)ELD(M)=ANS(G(M))                                             
      DO 70 I=1,NGP                                                             
      DO 70 J=1,NGP                                                             
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
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
