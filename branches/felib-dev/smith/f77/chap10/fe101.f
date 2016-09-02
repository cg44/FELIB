      PROGRAM P101                                                              
C                                                                               
C      PROGRAM 10.1 EIGENVALUES OF A RECTANGULAR                                
C      SOLID IN PLANE STRAIN USING 4-NODE QUADRILATERALS                        
C      LUMPED MASS                                                              
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKU1=200,IKU2=50,INF=100)                                       
C                                                                               
      REAL DEE(3,3),SAMP(3,2),COORD(4,2),FUN(4),JAC(2,2),JAC1(2,2),             
     +DER(2,4),DERIV(2,4),BEE(3,8),DBEE(3,8),BTDB(8,8),BT(8,3),                 
     +KM(8,8),KU(IKU1,IKU2),LOADS(IKU1),DIAG(IKU1),UDIAG(IKU1)                  
      INTEGER NF(INF,2),G(8)                                                    
      DATA IJAC,IJAC1,IDER,IDERIV,NODOF,IT/6*2/                                 
      DATA IH,ISAMP,IDEE,IBEE,IDBEE/5*3/                                        
      DATA ICOORD,NOD/2*4/,IBTDB,IKM,IBT,IDOF/4*8/                              
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,AA,BB,RHO,E,V                             
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IWP1=IW+1                                                                 
      CALL NULL(KU,IKU1,N,IWP1)                                                 
      CALL FMDEPS(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT STIFFNESS AND MASS INTEGRATION AND ASSEMBLY                      
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEOM4Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
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
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
   20 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   10 CALL FORMKU(KU,IKU1,KM,IKM,G,IW,IDOF)                                     
C                                                                               
C      FORM LUMPED MASS MATRIX                                                  
C                                                                               
      DO 30 I=1,N                                                               
   30 DIAG(I)=.5*AA*BB*RHO                                                      
      DO 40 I=9,12                                                              
   40 DIAG(I)=DIAG(I)*.5                                                        
      CALL PRINTV(DIAG,N)                                                       
C                                                                               
C      REDUCE TO STANDARD EIGENVALUE PROBLEM                                    
C                                                                               
      DO 50 I=1,N                                                               
   50 DIAG(I)=1./SQRT(DIAG(I))                                                  
      DO 60 I=1,N                                                               
      IF(I.LE.N-IW)K=IWP1                                                       
      IF(I.GT.N-IW)K=N-I+1                                                      
      DO 60 J=1,K                                                               
   60 KU(I,J)=KU(I,J)*DIAG(I)*DIAG(I+J-1)                                       
C                                                                               
C      EXTRACT EIGENVALUES                                                      
C                                                                               
      CALL BANDRD(N,IW,KU,IKU1,DIAG,UDIAG,LOADS)                                
      IFAIL=1                                                                   
      CALL BISECT(N,1.E-20,DIAG,UDIAG,IFAIL)                                    
      CALL PRINTV(DIAG,N)                                                       
      STOP                                                                      
      END                                                                       
