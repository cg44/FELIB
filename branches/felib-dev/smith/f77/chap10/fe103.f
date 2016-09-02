      PROGRAM P103                                                              
C                                                                               
C      PROGRAM 10.3 EIGENVALUES AND EIGENVECTORS OF A RECTANGULAR               
C      SOLID IN PLANE STRESS USING 8-NODE QUADRILATERALS                        
C      CONSISTENT MASS                                                          
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=103,IKB2=16,INF=85)                                        
C                                                                               
      REAL DEE(3,3),SAMP(3,2),COORD(8,2),FUN(8),JAC(2,2),JAC1(2,2),             
     +DER(2,8),DERIV(2,8),BEE(3,16),DBEE(3,16),BTDB(16,16),BT(16,3),            
     +EMM(16,16),ECM(16,16),KM(16,16),TN(16,16),NT(16,2),                       
     +BIGK(IKB1,IKB1),DIAG(IKB1),UDIAG(IKB1),KB(IKB1,IKB2),MB(IKB1,IKB2)        
      INTEGER NF(INF,2),G(16)                                                   
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/                                 
      DATA IH,ISAMP,IDEE,IBEE,IDBEE/5*3/                                        
      DATA ICOORD,NOD/2*8/,IBTDB,IKM,IBT,IEMM,IECM,IDOF,ITN,INT/8*16/           
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGPK,NGPM,AA,BB,RHO,E,V,NMODES                
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL NULL(MB,IKB1,N,IWP1)                                                 
      CALL NULL(BIGK,IKB1,N,N)                                                  
      CALL FMDSIG(DEE,IDEE,E,V)                                                 
C                                                                               
C      ELEMENT STIFFNESS AND MASS INTEGRATION AND ASSEMBLY                      
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEOM8Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      CALL NULL(EMM,IEMM,IDOF,IDOF)                                             
      CALL GAUSS(SAMP,ISAMP,NGPK)                                               
      DO 20 I=1,NGPK                                                            
      DO 20 J=1,NGPK                                                            
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
   20 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
      CALL GAUSS(SAMP,ISAMP,NGPM)                                               
      DO 30 I=1,NGPM                                                            
      DO 30 J=1,NGPM                                                            
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      QUOT=DET*SAMP(I,2)*SAMP(J,2)*RHO                                          
      CALL ECMAT(ECM,IECM,TN,ITN,NT,INT,FUN,NOD,NODOF)                          
      CALL MSMULT(ECM,IECM,QUOT,IDOF,IDOF)                                      
   30 CALL MATADD(EMM,IEMM,ECM,IECM,IDOF,IDOF)                                  
      CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
   10 CALL FORMKB(MB,IKB1,EMM,IEMM,G,IW,IDOF)                                   
C                                                                               
C      REDUCE TO STANDARD EIGENVALUE PROBLEM                                    
C                                                                               
      CALL CHOLIN(MB,IKB1,N,IW)                                                 
      CALL LBKBAN(MB,BIGK,KB,IKB1,IW,N)                                         
      DO 40 I=1,N                                                               
      DO 40 J=I+1,N                                                             
   40 BIGK(I,J)=BIGK(J,I)                                                       
      CALL LBBT(MB,BIGK,IKB1,IW,N)                                              
C                                                                               
C      EXTRACT EIGENVALUES AND EIGENVECTORS                                     
C                                                                               
      CALL TRIDIA(N,1.E-20,BIGK,DIAG,UDIAG,BIGK,IKB1)                           
      IFAIL=1                                                                   
      CALL EVECTS(N,1.E-20,DIAG,UDIAG,BIGK,IKB1,IFAIL)                          
      CALL PRINTV(DIAG,N)                                                       
      DO 50 J=1,NMODES                                                          
      DO 60 I=1,N                                                               
   60 UDIAG(I)=BIGK(I,J)                                                        
      CALL CHOBK2(MB,IKB1,UDIAG,N,IW)                                           
   50 CALL PRINTV(UDIAG,N)                                                      
      STOP                                                                      
      END                                                                       
