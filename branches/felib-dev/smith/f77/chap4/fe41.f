      PROGRAM P41                                                               
C                                                                               
C      PROGRAM 4.1  EQUILIBRIUM OF STEPPED BEAMS                                
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,IPROP=20,INO=20,INF=100)                     
C                                                                               
      REAL KM(4,4),ELD(4),ACTION(4),KV(IKV),LOADS(ILOADS),VAL(INO),             
     +PROP(IPROP,2)                                                             
      INTEGER G(4),NF(INF,2),NO(INO)                                            
      DATA IKM,IDOF/2*4/,NODOF/2/,IW/3/                                         
C                                                                               
C      INPUT SECTION                                                            
C                                                                               
      READ(5,*)NXE,N,NN,NR                                                      
      IR=(IW+1)*N                                                               
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C      NODE FREEDOM DATA                                                        
C                                                                               
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
C                                                                               
C      GLOBAL STIFFNESS MATRIX ASSEMBLY                                         
C                                                                               
      DO 10 IP=1,NXE                                                            
      READ(5,*)EI,ELL                                                           
      PROP(IP,1)=EI                                                             
      PROP(IP,2)=ELL                                                            
      CALL BEAMKM(KM,EI,ELL)                                                    
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
   10 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)IFIX,(NO(I),VAL(I),I=1,IFIX)                                     
      DO 20 I=1,IFIX                                                            
      KV(NO(I))=KV(NO(I))+1.E20                                                 
   20 LOADS(NO(I))=KV(NO(I))*VAL(I)                                             
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL BANRED(KV,N,IW)                                                      
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
C                                                                               
C      RETRIEVE ELEMENT END FORCES AND MOMENTS                                  
C                                                                               
      DO 30 IP=1,NXE                                                            
      EI=PROP(IP,1)                                                             
      ELL=PROP(IP,2)                                                            
      CALL BEAMKM(KM,EI,ELL)                                                    
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      DO 40 I=1,IDOF                                                            
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   40 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
   30 CALL PRINTV(ACTION,IDOF)                                                  
      STOP                                                                      
      END                                                                       
