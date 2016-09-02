      PROGRAM P90                                                               
C                                                                               
C      PROGRAM 9.0 STEADY STATE NAVIER-STOKES                                   
C      USING 8 NODE VELOCITY ELEMENTS                                           
C      AND   4 NODE PRESSURE ELEMENTS                                           
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IPB1=100,IPB2=50,INF=100,INO=30,                                
     +          IWID=30,IDEP=30)                                                
C                                                                               
      REAL SAMP(3,2),COORD(8,2),DERIVT(8,2),UVEL(8),VVEL(8),                    
     +JAC(2,2),JAC1(2,2),KAY(2,2),DER(2,8),DERIV(2,8),                          
     +DTKD(8,8),KDERIV(2,8),KE(20,20),FUN(8),OLDLDS(IPB1),                      
     +FUNF(4),COORDF(4,2),DERF(2,4),DERIVF(2,4),                                
     +WIDTH(IWID),DEPTH(IDEP),VAL(INO),ROW(8),TEMP(8,8),                        
     +C11(8,8),C12(8,4),C21(4,8),C23(4,8),C32(8,4),                             
     +PB(IPB1,IPB2),WORK(IPB2,IPB1),LOADS(IPB1)                                 
      INTEGER G(20),NO(INO),NF(INF,3)                                           
      DATA IT,IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV,IDERF,IDERVF/9*2/              
      DATA ISAMP,NODOF/2*3/,ICORDF,IC21,IC23,NODF/4*4/,IKE,ITOT/2*20/           
      DATA ICOORD,IDERVT,IDTKD,ITEMP,IC11,IC12,IC32,NOD/8*8/                    
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,VISC,RHO,ITS                              
      READ(5,*)(WIDTH(I),I=1,NXE+1)                                             
      READ(5,*)(DEPTH(I),I=1,NYE+1)                                             
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      READ(5,*)IFIX,(NO(I),VAL(I),I=1,IFIX)                                     
      IWP1=IW+1                                                                 
      IBAND=2*IWP1-1                                                            
      CALL NULVEC(UVEL,NOD)                                                     
      CALL NULVEC(VVEL,NOD)                                                     
      CALL NULVEC(OLDLDS,N)                                                     
      CALL NULVEC(LOADS,N)                                                      
      CALL NULL(KAY,IKAY,IT,IT)                                                 
      KAY(1,1)=VISC/RHO                                                         
      KAY(2,2)=VISC/RHO                                                         
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
      ITERS=0                                                                   
   10 ITERS=ITERS+1                                                             
      CALL NULL(PB,IPB1,N,IBAND)                                                
      CALL NULL(WORK,IPB2,IWP1,N)                                               
      CALL NULL(KE,IKE,ITOT,ITOT)                                               
C                                                                               
C      ELEMENT MATRIX INTEGRATION AND ASSEMBLY                                  
C                                                                               
      DO 20 IP=1,NXE                                                            
      DO 20 IQ=1,NYE                                                            
      CALL GEVUPV(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,COORDF,ICORDF,             
     +            G,NF,INF)                                                     
      DO 30 M=1,NOD                                                             
      IF(G(M).EQ.0)UVEL(M)=0.                                                   
   30 IF(G(M).NE.0)UVEL(M)=(LOADS(G(M))+OLDLDS(G(M)))*.5                        
      DO 40 M=NOD+NODF+1,ITOT                                                   
      IF(G(M).EQ.0)VVEL(M-NOD-NODF)=0.0                                         
   40 IF(G(M).NE.0)VVEL(M-NOD-NODF)=(LOADS(G(M))+ OLDLDS(G(M)))*.5              
      CALL NULL(C11,IC11,NOD,NOD)                                               
      CALL NULL(C12,IC12,NOD,NODF)                                              
      CALL NULL(C21,IC21,NODF,NOD)                                              
      CALL NULL(C23,IC23,NODF,NOD)                                              
      CALL NULL(C32,IC32,NOD,NODF)                                              
      DO 50 I=1,NGP                                                             
      DO 50 J=1,NGP                                                             
C                                                                               
C      VELOCITY CONTRIBUTION                                                    
C                                                                               
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      UBAR=0.                                                                   
      VBAR=0.                                                                   
      DO 60 M=1,NOD                                                             
      UBAR=UBAR+FUN(M)*UVEL(M)                                                  
   60 VBAR=VBAR+FUN(M)*VVEL(M)                                                  
      IF(ITERS.EQ.1)UBAR=1.0                                                   
      IF(ITERS.EQ.1)VBAR=1.0                                                   
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL MATMUL(KAY,IKAY,DERIV,IDERIV,KDERIV,IKDERV,IT,IT,NOD)                
      CALL MATRAN(DERIVT,IDERVT,DERIV,IDERIV,IT,NOD)                            
      CALL MATMUL(DERIVT,IDERVT,KDERIV,IKDERV,DTKD,IDTKD,NOD,IT,NOD)            
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      DO 70 K=1,NOD                                                             
      DO 70 L=1,NOD                                                             
   70 DTKD(K,L)=DTKD(K,L)*QUOT                                                  
      CALL MATADD(C11,IC11,DTKD,IDTKD,NOD,NOD)                                  
      DO 80 K=1,NOD                                                             
   80 ROW(K)=DERIV(1,K)                                                         
      PROD=QUOT*UBAR                                                            
      CALL VVMULT(FUN,ROW,TEMP,ITEMP,NOD,NOD)                                   
      CALL MSMULT(TEMP,ITEMP,PROD,NOD,NOD)                                      
      CALL MATADD(C11,IC11,TEMP,ITEMP,NOD,NOD)                                  
      DO 90 K=1,NOD                                                             
   90 ROW(K)=DERIV(2,K)                                                         
      PROD=QUOT*VBAR                                                            
      CALL VVMULT(FUN,ROW,TEMP,ITEMP,NOD,NOD)                                   
      CALL MSMULT(TEMP,ITEMP,PROD,NOD,NOD)                                      
      CALL MATADD(C11,IC11,TEMP,ITEMP,NOD,NOD)                                  
C                                                                               
C      PRESSURE CONTRIBUTION                                                    
C                                                                               
      CALL FORMLN(DERF,IDERF,FUNF,SAMP,ISAMP,I,J)                               
      CALL MATMUL(DERF,IDERF,COORDF,ICORDF,JAC,IJAC,IT,NODF,IT)                 
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DERF,IDERF,DERIVF,IDERVF,IT,IT,NODF)               
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      PROD=QUOT/RHO                                                             
      DO 100 K=1,NODF                                                           
  100 ROW(K)=DERIVF(1,K)                                                        
      CALL VVMULT(FUN,ROW,TEMP,ITEMP,NOD,NODF)                                  
      CALL MSMULT(TEMP,ITEMP,PROD,NOD,NODF)                                     
      CALL MATADD(C12,IC12,TEMP,ITEMP,NOD,NODF)                                 
      DO 110 K=1,NODF                                                           
  110 ROW(K)=DERIVF(2,K)                                                        
      CALL VVMULT(FUN,ROW,TEMP,ITEMP,NOD,NODF)                                  
      CALL MSMULT(TEMP,ITEMP,PROD,NOD,NODF)                                     
      CALL MATADD(C32,IC32,TEMP,ITEMP,NOD,NODF)                                 
      DO 120 K=1,NOD                                                            
  120 ROW(K)=DERIV(1,K)                                                         
      CALL VVMULT(FUNF,ROW,TEMP,ITEMP,NODF,NOD)                                 
      CALL MSMULT(TEMP,ITEMP,QUOT,NODF,NOD)                                     
      CALL MATADD(C21,IC21,TEMP,ITEMP,NODF,NOD)                                 
      DO 130 K=1,NOD                                                            
  130 ROW(K)=DERIV(2,K)                                                         
      CALL VVMULT(FUNF,ROW,TEMP,ITEMP,NODF,NOD)                                 
      CALL MSMULT(TEMP,ITEMP,QUOT,NODF,NOD)                                     
      CALL MATADD(C23,IC23,TEMP,ITEMP,NODF,NOD)                                 
   50 CONTINUE                                                                  
      CALL FRMUPV(KE,IKE,C11,IC11,C12,IC12,C21,IC21,                            
     +            C23,IC23,C32,IC32,NOD,NODF,ITOT)                              
      CALL FORMTB(PB,IPB1,KE,IKE,G,IW,ITOT)                                     
   20 CONTINUE                                                                  
C                                                                               
C      INSERT PRESCRIBED VALUES OF VELOCITY AND PRESSURE                        
C                                                                               
      CALL NULVEC(LOADS,N)                                                      
      DO 140 I=1,IFIX                                                           
      PB(NO(I),IWP1)=PB(NO(I),IWP1)+1.E20                                       
  140 LOADS(NO(I))=PB(NO(I),IWP1)*VAL(I)                                        
C                                                                               
C      SOLVE SIMULTANEOUS EQUATIONS AND CHECK CONVERGENCE                       
C                                                                               
      CALL GAUSBA(PB,IPB1,WORK,IPB2,N,IW)                                       
      CALL SOLVBA(PB,IPB1,WORK,IPB2,LOADS,N,IW)                                 
      CALL CHECON(LOADS,OLDLDS,N,.001,ICON)                                     
      IF(ICON.EQ.0)GOTO 10                                                      
      CALL PRINTV(LOADS,N)                                                      
      WRITE(6,1000)ITERS                                                        
 1000 FORMAT(I10)                                                               
      STOP                                                                      
      END                                                                       
