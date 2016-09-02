      PROGRAM P45                                                               
C                                                                               
C     PROGRAM 4.5  ANALYSIS OF 2-D TRUSSES                                      
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,IPROP=20)                                    
C                                                                               
      REAL KM(4,4),ELD(4),ACTION(4),                                            
     +KV(IKV),LOADS(ILOADS),PROP(IPROP),COORD(IPROP,4)                          
      INTEGER G(4),STOREG(IPROP,4)                                              
      DATA IKM,IDOF/2*4/,NODOF/2/                                               
C                                                                               
C      INPUT SECTION                                                            
C                                                                               
      READ(5,*)NXE,N,IW                                                         
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C      GLOBAL STIFFNESS MATRIX ASSEMBLY                                         
C                                                                               
      DO 10 IP=1,NXE                                                            
      READ(5,*)EA,(COORD(IP,I),I=1,4),(G(I),I=1,IDOF)                           
      PROP(IP)=EA                                                               
      CALL PINJ2(KM,EA,IP,COORD,IPROP)                                          
      DO 20 I=1,IDOF                                                            
   20 STOREG(IP,I)=G(I)                                                         
   10 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C       EQUATION SOLUTION                                                       
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
      CALL NULVEC(LOADS,N)                                                      
      READ(5,*)NL,(K,LOADS(K),I=1,NL)                                           
      CALL BACSUB(KV,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
C                                                                               
C      RETRIEVE ELEMENT AXIAL LOADS                                             
C                                                                               
      DO 30 IP=1,NXE                                                            
      EA=PROP(IP)                                                               
      CALL PINJ2(KM,EA,IP,COORD,IPROP)                                          
      DO 40 I=1,IDOF                                                            
      G(I)=STOREG(IP,I)                                                         
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   40 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
      CALL PRINTV(ACTION,IDOF)                                                  
      CALL LOC2T(AXIAL,ACTION,IP,COORD,IPROP)                                   
   30 WRITE(6,'(E12.4)')AXIAL                                                   
      STOP                                                                      
      END                                                                       
