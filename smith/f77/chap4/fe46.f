      PROGRAM P46                                                               
C                                                                               
C     PROGRAM 4.6  ANALYSIS OF 3-D TRUSSES                                      
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,IPROP=20)                                    
C                                                                               
      REAL KM(6,6),ELD(6),ACTION(6),                                            
     +KV(IKV),LOADS(ILOADS),PROP(IPROP),COORD(IPROP,6)                          
      INTEGER G(6),STOREG(IPROP,6)                                              
      DATA IKM,IDOF/2*6/,NODOF/3/                                               
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
      READ(5,*)EA,(COORD(IP,I),I=1,6),(G(I),I=1,IDOF)                           
      PROP(IP)=EA                                                               
      CALL PINJ3(KM,EA,IP,COORD,IPROP)                                          
      DO 20 I=1,IDOF                                                            
   20 STOREG(IP,I)=G(I)                                                         
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
C      RETRIEVE ELEMENT AXIAL LOADS                                             
C                                                                               
      DO 30 IP=1,NXE                                                            
      EA=PROP(IP)                                                               
      CALL PINJ3(KM,EA,IP,COORD,IPROP)                                          
      DO 40 I=1,IDOF                                                            
      G(I)=STOREG(IP,I)                                                         
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   40 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
      CALL PRINTV(ACTION,IDOF)                                                  
      CALL LOC3T(AXIAL,ACTION,IP,COORD,IPROP)                                   
   30 WRITE(6,'(E12.4)')AXIAL                                                   
      STOP                                                                      
      END                                                                       
