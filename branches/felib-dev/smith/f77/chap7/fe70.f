      PROGRAM P70                                                               
C                                                                               
C      PROGRAM 7.0 SOLUTION OF LAPLACE'S EQUATION OVER A                        
C      PLANE AREA USING 4-NODE QUADRILATERALS                                   
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=1000,ILOADS=150,INF=50,INO=10)                              
C                                                                               
      REAL JAC(2,2),JAC1(2,2),KAY(2,2),SAMP(3,2),DTKD(4,4),KP(4,4),             
     +COORD(4,2),DER(2,4),DERIV(2,4),DERIVT(4,2),KDERIV(2,4),FUN(4),            
     +VAL(INO),KVH(IKV),KV(IKV),LOADS(ILOADS),DISPS(ILOADS)                     
      INTEGER G(4),NO(INO),NF(INF,1)                                            
      DATA IT,IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV/7*2/,ISAMP/3/                  
      DATA IDTKD,IKP,ICOORD,IDERVT,NOD/5*4/,NODOF/1/                            
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,AA,BB,PERMX,PERMY                         
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL NULL(KAY,IKAY,IT,IT)                                                 
      KAY(1,1)=PERMX                                                            
      KAY(2,2)=PERMY                                                            
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT INTEGRATION AND ASSEMBLY                                         
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KP,IKP,NOD,NOD)                                                 
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL MATMUL(KAY,IKAY,DERIV,IDERIV,KDERIV,IKDERV,IT,IT,NOD)                
      CALL MATRAN(DERIVT,IDERVT,DERIV,IDERIV,IT,NOD)                            
      CALL MATMUL(DERIVT,IDERVT,KDERIV,IKDERV,DTKD,IDTKD,NOD,IT,NOD)            
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      CALL MSMULT(DTKD,IDTKD,QUOT,NOD,NOD)                                      
   20 CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
   10 CALL FORMKV(KV,KP,IKP,G,N,NOD)                                            
      CALL VECCOP(KV,KVH,IR)                                                    
C                                                                               
C      SPECIFY FIXED POTENTIALS AND REDUCE EQUATIONS                            
C                                                                               
      READ(5,*)IFIX,(NO(I),VAL(I),I=1,IFIX)                                     
      CALL NULVEC(LOADS,N)                                                      
      DO 30 I=1,IFIX                                                            
      KV(NO(I))=KV(NO(I))+1.E20                                                 
   30 LOADS(NO(I))=KV(NO(I))*VAL(I)                                             
      CALL BANRED(KV,N,IW)                                                      
C                                                                               
C      SOLVE EQUATIONS AND RETRIEVE FLOW RATE                                   
C                                                                               
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
      CALL LINMUL(KVH,LOADS,DISPS,N,IW)                                         
      REACT=0.                                                                  
      DO 40 I=1,IFIX                                                            
   40 REACT=REACT+DISPS(NO(I))                                                  
      WRITE(6,'(E12.4)')REACT                                                   
      STOP                                                                      
      END                                                                       
