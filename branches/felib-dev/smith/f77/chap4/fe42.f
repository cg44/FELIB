      PROGRAM P42                                                               
C                                                                               
C      PROGRAM 4.2 NUMERICALLY INTEGRATED BEAM ON ELASTIC FOUNDATION            
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,INF=100,IPROP=15)                            
C                                                                               
      REAL KM(4,4),ELD(4),ACTION(4),KV(IKV),LOADS(ILOADS),MM(4,4),              
     +DTD(4,4),FTF(4,4),DER2(4),FUN(4),SAMP(7,2),MOM(INF),                      
     +STORKM(IPROP,4,4)                                                         
      INTEGER G(4),NF(INF,2)                                                    
      DATA IKM,IDOF/2*4/,NODOF/2/,IW/3/,ISAMP/7/                                
C                                                                               
C      INPUT SECTION                                                            
C                                                                               
      READ(5,*)NXE,N,NN,NR,NGP,EI0,EI1,FS0,FS1,ELL                              
      IR=(IW+1)*N                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      NODE FREEDOM DATA                                                        
C                                                                               
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
C                                                                               
C      GLOBAL STIFFNESS AND MASS MATRIX ASSEMBLY                                
C                                                                               
      X=0.                                                                      
      DO 10 IP=1,NXE                                                            
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      CALL NULL(MM,IKM,IDOF,IDOF)                                               
      DO 20 I=1,NGP                                                             
      SP=X+ELL*.5*(SAMP(I,1)+1.)                                                
      EI=SP/(ELL*NXE)*(EI1-EI0)+EI0                                             
      FS=SP/(ELL*NXE)*(FS1-FS0)+FS0                                             
      WT=SAMP(I,2)                                                              
      CALL FMBEAM(DER2,FUN,SAMP,ISAMP,ELL,I)                                    
      DO 30 K=1,IDOF                                                            
      DO 30 L=1,IDOF                                                            
      FTF(K,L)=FUN(K)*FUN(L)*WT*.5*ELL*FS                                       
   30 DTD(K,L)=DER2(K)*DER2(L)*WT*8.*EI/(ELL*ELL*ELL)                           
      CALL MATADD(MM,IKM,FTF,IKM,IDOF,IDOF)                                     
   20 CALL MATADD(KM,IKM,DTD,IKM,IDOF,IDOF)                                     
      CALL MATADD(KM,IKM,MM,IKM,IDOF,IDOF)                                      
      DO 40 I=1,IDOF                                                            
      DO 40 J=1,IDOF                                                            
   40 STORKM(IP,I,J)=KM(I,J)                                                    
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      X=X+ELL                                                                   
   10 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL BACSUB(KV,LOADS,N,IW)                                                
C                                                                               
C      RETRIEVE ELEMENT END FORCES AND MOMENTS                                  
C                                                                               
      DO 50 IP=1,NXE                                                            
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      DO 60 I=1,IDOF                                                            
      DO 60 J=1,IDOF                                                            
   60 KM(I,J)=STORKM(IP,I,J)                                                    
      DO 70 I=1,IDOF                                                            
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   70 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
      MOM(IP)=ACTION(2)                                                         
      IF(IP.EQ.NXE)MOM(IP+1)=-ACTION(4)                                         
   50 CONTINUE                                                                  
      DO 80 I=1,NN                                                              
      WRITE(6,'(2E12.4)')LOADS(2*I-1),MOM(I)                                    
   80 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
