      PROGRAM P49                                                               
C                                                                               
C      PROGRAM 4.9 EQUILIBRIUM OF TRANSVERSLY LOADED                            
C      RECTANGULAR PLATES                                                       
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,INF=100)                                     
C                                                                               
      REAL KV(IKV),LOADS(ILOADS),KM(16,16),DTD(16,16),FUN(16),                  
     +D1X(16),D2X(16),D1Y(16),D2Y(16),D2XY(16),SAMP(7,2)                        
      INTEGER NF(INF,4),G(16)                                                   
      DATA NODOF/4/,ISAMP/7/,IDOF,IKM,IDTD/3*16/                                
C                                                                               
C      INPUT SECTION                                                            
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,E,V,TH,AA,BB                              
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      D=E*TH**3/(12.*(1.-V*V))                                                  
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      FORM ELEMENT STIFFNESS MATRIX BY NUMERICAL INTEGRATION                   
C                                                                               
      DO 10 I=1,NGP                                                             
      DO 10 J=1,NGP                                                             
      CALL FMPLAT(FUN,D1X,D1Y,D2X,D2Y,D2XY,SAMP,ISAMP,AA,BB,I,J)                
      DO 20 K=1,IDOF                                                            
      DO 20 L=K,IDOF                                                            
      DTD(K,L)=4.*AA*BB*D*SAMP(I,2)*SAMP(J,2)*(D2X(K)*D2X(L)/(AA**4)+           
     +         D2Y(K)*D2Y(L)/(BB**4)+(V*D2X(K)*D2Y(L)+                          
     +         V*D2X(L)*D2Y(K)+2.*(1.-V)*D2XY(K)*D2XY(L))/(AA**2*BB**2))        
   20 DTD(L,K)=DTD(K,L)                                                         
   10 CALL MATADD(KM,IKM,DTD,IDTD,IDOF,IDOF)                                    
C                                                                               
C      GLOBAL STIFFNESS MATRIX ASSEMBLY                                         
C                                                                               
      DO 30 IP=1,NXE                                                            
      DO 30 IQ=1,NYE                                                            
      CALL FORMGP(IP,IQ,NYE,G,NF,INF)                                           
   30 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
      STOP                                                                      
      END                                                                       
