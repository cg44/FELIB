      PROGRAM P40                                                               
C                                                                               
C      PROGRAM 4.0 EQUILIBRIUM OF UNIFORM BEAMS                                 
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,INF=100)                                     
C                                                                               
      REAL KM(4,4),ELD(4),ACTION(4),KV(IKV),LOADS(ILOADS)                       
      INTEGER G(4),NF(INF,2)                                                    
      DATA IKM,IDOF/2*4/,NODOF/2/,IW/3/                                         
C                                                                               
C      INPUT SECTION                                                            
C                                                                               
C     open(5,file='book40.d')
C     open(6,file='book40.r')
      READ(5,*)NXE,N,NN,NR,EI,ELL                                               
      IR=(IW+1)*N                                                               
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C      NODE FREEDOM DATA                                                        
C                                                                               
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
C                                                                               
C      ELEMENT STIFFNESS MATRIX                                                 
C                                                                               
      CALL BEAMKM(KM,EI,ELL)                                                    
C                                                                               
C      GLOBAL STIFFNESS MATRIX ASSEMBLY                                         
C                                                                               
      DO 10 IP=1,NXE                                                            
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
   10 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
C                                                                               
C      RETRIEVE ELEMENT END FORCES AND MOMENTS                                  
C                                                                               
      DO 20 IP=1,NXE                                                            
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      DO 30 I=1,IDOF                                                            
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   30 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
   20 CALL PRINTV(ACTION,IDOF)                                                  
      STOP                                                                      
      END                                                                       
