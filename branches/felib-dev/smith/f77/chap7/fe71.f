      PROGRAM P71                                                               
C                                                                               
C      PROGRAM 7.1 SOLUTION OF LAPLACE'S EQUATION FOR                           
C      PLANE FREE-SURFACE FLOW USING 4-NODED QUADRILATERALS                     
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=1000,ILOADS=150,INF=70,INO=12,INX=10)                       
C                                                                               
      REAL JAC(2,2),JAC1(2,2),KAY(2,2),SAMP(3,2),DTKD(4,4),KP(4,4),             
     +COORD(4,2),DER(2,4),DERIV(2,4),DERIVT(4,2),KDERIV(2,4),FUN(4),            
     +VAL(INO),KVH(IKV),KV(IKV),LOADS(ILOADS),DISPS(ILOADS),                    
     +OLDPOT(ILOADS),WIDTH(INX),SURF(INX)                                       
      INTEGER G(4),NO(INO),NF(INF,1)                                            
      DATA IT,IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV/7*2/,ISAMP/3/                  
      DATA IDTKD,IKP,ICOORD,IDERVT,NOD/5*4/,NODOF/1/                            
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,ITS,PERMX,PERMY                           
      READ(5,*)(WIDTH(I),I=1,NXE+1)                                             
      READ(5,*)(SURF(I),I=1,NXE+1)                                              
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IR=N*(IW+1)                                                               
      CALL NULVEC(OLDPOT,N)                                                     
      CALL NULL(KAY,IKAY,IT,IT)                                                 
      KAY(1,1)=PERMX                                                            
      KAY(2,2)=PERMY                                                            
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ITERATE FOR POSITION OF FREESURFACE                                      
C                                                                               
      ITERS=0                                                                   
   10 ITERS=ITERS+1                                                             
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C      ELEMENT INTEGRATION AND ASSEMBLY                                         
C                                                                               
      DO 20 IP=1,NXE                                                            
      DO 20 IQ=1,NYE                                                            
      CALL WELGEO(IP,IQ,NXE,NYE,WIDTH,SURF,COORD,ICOORD,G,NF,INF)               
      CALL NULL(KP,IKP,NOD,NOD)                                                 
      DO 30 I=1,NGP                                                             
      DO 30 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL MATMUL(KAY,IKAY,DERIV,IDERIV,KDERIV,IKDERV,IT,IT,NOD)                
      CALL MATRAN(DERIVT,IDERVT,DERIV,IDERIV,IT,NOD)                            
      CALL MATMUL(DERIVT,IDERVT,KDERIV,IKDERV,DTKD,IDTKD,NOD,IT,NOD)            
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      CALL MSMULT(DTKD,IDTKD,QUOT,NOD,NOD)                                      
   30 CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
   20 CALL FORMKV(KV,KP,IKP,G,N,NOD)                                            
      CALL VECCOP(KV,KVH,IR)                                                    
C                                                                               
C      SPECIFY FIXED POTENTIALS AND REDUCE EQUATIONS                            
C                                                                               
      IF(ITERS.EQ.1)READ(5,*)IFIX,(NO(I),I=1,IFIX)                              
      CALL NULVEC(LOADS,N)                                                      
      DO 40 I=1,IFIX                                                            
      KV(NO(I))=KV(NO(I))+1.E20                                                 
   40 LOADS(NO(I))=KV(NO(I))*SURF(NXE+1)                                        
      DO 50 IQ=1,NYE-1                                                          
      J=IQ*(NXE+1)+1                                                            
   50 LOADS(J)=KV(J)*(NYE-IQ)*SURF(1)/NYE                                       
      CALL BANRED(KV,N,IW)                                                      
C                                                                               
C      SOLVE EQUATIONS                                                          
C                                                                               
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL VECCOP(LOADS,SURF,NXE)                                               
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(LOADS,OLDPOT,N,0.001,ICON)                                    
      IF(ITERS.NE.ITS.AND.ICON.EQ.0)GOTO 10                                     
      CALL LINMUL(KVH,LOADS,DISPS,N,IW)                                         
      REACT=0.                                                                  
      DO 60 I=1,NYE                                                             
   60 REACT=REACT+DISPS(I*(NXE+1))                                              
      REACT=REACT+DISPS(N)                                                      
      CALL PRINTV(LOADS,N)                                                      
      WRITE(6,'(E12.4)')REACT                                                   
      WRITE(6,'(I10)')ITERS                                                     
      STOP                                                                      
      END                                                                       
