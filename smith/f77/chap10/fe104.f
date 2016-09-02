      PROGRAM P104                                                              
C                                                                               
C      PROGRAM 10.4 EIGENVALUES AND EIGENVECTORS OF A                           
C      RECTANGULAR SOLID IN PLANE STRAIN USING 4-NODE                           
C      QUADRILATERAL ELEMENTS : CONSISTENT MASS,LANCZOS METHOD                  
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=200,IKB2=20,LALFA=500,LEIG=20,LX=80,LY=200,LZ=500,         
     +          INF=110)                                                        
C                                                                               
      REAL DEE(3,3),SAMP(3,2),COORD(4,2),FUN(4),JAC(2,2),JAC1(2,2),             
     +DER(2,4),DERIV(2,4),BEE(3,8),DBEE(3,8),BTDB(8,8),KM(8,8),                 
     +BT(8,3),VOL(8),EMM(8,8),ECM(8,8),TN(8,8),NT(8,2),                         
     +UA(IKB1),VA(IKB1),EIG(LEIG),X(LX),DEL(LX),UDIAG(IKB1),                    
     +ALFA(LALFA),BETA(LALFA),W1(IKB1),Y(LY,LEIG),Z(LZ,LEIG),                   
     +KB(IKB1,IKB2),MB(IKB1,IKB2),DIAG(IKB1)                                    
      INTEGER G(8),NF(INF,2),NU(LX),JEIG(2,LEIG)                                
      DATA NODOF,IT,IJAC,IJAC1,IDER,IDERIV/6*2/,ICOORD,NOD/2*4/                 
      DATA IDEE,IBEE,IDBEE,IH,ISAMP/5*3/                                        
      DATA IBTDB,IKM,IBT,IEMM,IECM,ITN,INT,IDOF/8*8/                            
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      DATA EL/0./,ER/20./,ACC/1.E-6/,LP/6/,ITAPE/1/,IFLAG/-1/                   
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,AA,BB,RHO,E,V,NMODES                      
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL NULL(MB,IKB1,N,IWP1)                                                 
      CALL NULL(DEE,IDEE,IH,IH)                                                 
      CALL FMDEPS(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT STIFFNESS AND MASS INTEGRATION AND ASSEMBLY                      
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEOM4Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      CALL NULL(EMM,IEMM,IDOF,IDOF)                                             
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      CALL ECMAT(ECM,IECM,TN,ITN,NT,INT,FUN,NOD,NODOF)                          
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      DO 30 K=1,IDOF                                                            
      DO 30 L=1,IDOF                                                            
      ECM(K,L)=ECM(K,L)*QUOT*RHO                                                
   30 BTDB(K,L)=BTDB(K,L)*QUOT                                                  
      CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   20 CALL MATADD(EMM,IEMM,ECM,IECM,IDOF,IDOF)                                  
      CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
   10 CALL FORMKB(MB,IKB1,EMM,IEMM,G,IW,IDOF)                                   
C                                                                               
C      FIND EIGENVALUES                                                         
C                                                                               
      CALL CHOLIN(MB,IKB1,N,IW)                                                 
      DO 40 ITERS=1,LALFA                                                       
      CALL LANCZ1(N,EL,ER,ACC,LEIG,LX,LALFA,LP,ITAPE,                           
     +            IFLAG,UA,VA,EIG,JEIG,NEIG,X,DEL,NU,                           
     +            ALFA,BETA)                                                    
      IF(IFLAG.EQ.0)GOTO 50                                                     
      IF(IFLAG.GT.1)GOTO 60                                                     
C                                                                               
C      IFLAG=1   FORM U+AV                                                      
C                                                                               
      CALL VECCOP(VA,UDIAG,N)                                                   
      CALL CHOBK2(MB,IKB1,UDIAG,N,IW)                                           
      CALL BANMUL(KB,IKB1,UDIAG,DIAG,N,IW)                                      
      CALL CHOBK1(MB,IKB1,DIAG,N,IW)                                            
   40 CALL VECADD(UA,DIAG,UA,N)                                                 
      GOTO 60                                                                   
C                                                                               
C      WRITE OUT SPECTRUM FOUND                                                 
C                                                                               
   50 WRITE(6,1000)ITERS,EL,ER                                                  
      CALL PRINTV(EIG,NEIG)                                                     
C                                                                               
C      CALCULATE EIGENVECTORS                                                   
C                                                                               
      NEXQT=1                                                                   
      IF(NEXQT.EQ.0)GOTO 70                                                     
      IF(NEIG.GT.10)NEIG=10                                                     
      CALL LANCZ2(N,LALFA,LP,ITAPE,EIG,JEIG,NEIG,ALFA,BETA,                     
     +            LY,LZ,JFLAG,Y,W1,Z)                                           
      IF(JFLAG.NE.0)GOTO 80                                                     
C                                                                               
C      JFLAG=0  CALCULATE EIGENVECTORS                                          
C                                                                               
      DO 90 I=1,NMODES                                                          
      DO 100 J=1,N                                                              
  100 UDIAG(J)=Y(J,I)                                                           
      CALL CHOBK2(MB,IKB1,UDIAG,N,IW)                                           
      CALL PRINTV(UDIAG,N)                                                      
   90 CONTINUE                                                                  
      GOTO 70                                                                   
C                                                                               
C      LANCZ1 IS SIGNALLING FAILURE                                             
C                                                                               
   60 WRITE(6,2000)IFLAG                                                        
      GOTO 70                                                                   
C                                                                               
C      EAI5ED IS SIGNALLING FAILURE                                             
C                                                                               
   80 WRITE(6,3000)JFLAG                                                        
   70 CONTINUE                                                                  
 1000 FORMAT(I10,2E12.4)                                                        
 2000 FORMAT(26X,'LANCZ1 HAS FAILED. IFLAG=',I2)                                
 3000 FORMAT(26X,'LANCZ2 HAS FAILED. JFLAG=',I2)                                
      STOP                                                                      
      END                                                                       
