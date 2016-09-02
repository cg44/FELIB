      PROGRAM P102                                                              
C                                                                               
C      PROGRAM 10.2 EIGENVALUES AND EIGENVECTORS OF A RECTANGULAR               
C      SOLID IN PLANE STRAIN USING 8-NODE QUADRILATERALS                        
C      LUMPED MASS                                                              
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IBIGK=103,INF=85)                                               
C                                                                               
      REAL DEE(3,3),SAMP(3,2),COORD(8,2),FUN(8),JAC(2,2),JAC1(2,2),             
     +DER(2,8),DERIV(2,8),BEE(3,16),DBEE(3,16),BTDB(16,16),BT(16,3),            
     +BIGK(IBIGK,IBIGK),LOADS(IBIGK),DIAG(IBIGK),UDIAG(IBIGK),                  
     +KM(16,16),EMM(16,16)                                                      
      INTEGER NF(INF,2),G(16)                                                   
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/                                 
      DATA IH,ISAMP,IDEE,IBEE,IDBEE/5*3/                                        
      DATA ICOORD,NOD/2*8/,IBTDB,IKM,IBT,IEMM,IDOF/5*16/                        
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,NN,NR,NGP,AA,BB,RHO,E,V,NMODES                         
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      CALL NULL(BIGK,IBIGK,N,N)                                                 
      CALL NULL(DEE,IDEE,IH,IH)                                                 
      CALL FMDEPS(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      CALL NULL(EMM,IEMM,IDOF,IDOF)                                             
      CALL NULVEC(DIAG,N)                                                       
C                                                                               
C      FORM LUMPED MASS MATRIX                                                  
C                                                                               
      DO 30 I=1,IDOF                                                            
   30 EMM(I,I)=AA*BB*RHO*.2                                                     
      DO 40 I=1,13,4                                                            
   40 EMM(I,I)=EMM(3,3)*.25                                                     
      DO 50 I=2,14,4                                                            
   50 EMM(I,I)=EMM(3,3)*.25                                                     
C                                                                               
C      ELEMENT STIFFNESS AND MASS INTEGRATION AND ASSEMBLY                      
C                                                                               
      DO 60 IP=1,NXE                                                            
      DO 60 IQ=1,NYE                                                            
      CALL GEOM8Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      DO 70 I=1,NGP                                                             
      DO 70 J=1,NGP                                                             
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
   70 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
      CALL FMLUMP(DIAG,EMM,IEMM,G,IDOF)                                         
   60 CALL FMBIGK(BIGK,IBIGK,KM,IKM,G,IDOF)                                     
      CALL PRINTV(DIAG,N)                                                       
C                                                                               
C      REDUCE TO STANDARD EIGENVALUE PROBLEM                                    
C                                                                               
      DO 80 I=1,N                                                               
      DIAG(I)=1./SQRT(DIAG(I))                                                  
   80 LOADS(I)=DIAG(I)                                                          
      DO 90 I=1,N                                                               
      DO 90 J=1,N                                                               
   90 BIGK(I,J)=BIGK(I,J)*DIAG(I)*DIAG(J)                                       
C                                                                               
C      EXTRACT EIGENVALUES                                                      
C                                                                               
      CALL TRIDIA(N,1.E-20,BIGK,DIAG,UDIAG,BIGK,IBIGK)                          
      IFAIL=1                                                                   
      CALL EVECTS(N,1.E-20,DIAG,UDIAG,BIGK,IBIGK,IFAIL)                         
      CALL PRINTV(DIAG,N)                                                       
      DO 100 J=1,NMODES                                                         
      DO 110 I=1,N                                                              
  110 DIAG(I)=BIGK(I,J)*LOADS(I)                                                
  100 CALL PRINTV(DIAG,N)                                                       
      STOP                                                                      
      END                                                                       

