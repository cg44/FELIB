      PROGRAM P43                                                               
C                                                                               
C       PROGRAM 4.3 ANALYSIS OF 2-D FRAMES                                      
C                                                                               
C                                                                               
C       ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                  
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,IPROP=20)                                    
C                                                                               
      REAL ACTION(6),LOCAL(6),ELD(6),KM(6,6),                                   
     +LOADS(ILOADS),KV(IKV),PROP(IPROP,2),COORD(IPROP,4)                        
      INTEGER G(6),STOREG(IPROP,6)                                              
      DATA IKM,IDOF/2*6/,NODOF/3/                                               
C                                                                               
C       INPUT SECTION                                                           
C                                                                               
      READ(5,*)NXE,N,IW                                                         
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C       GLOBAL STIFFNESS MATRIX ASSEMBLY                                        
C                                                                               
      DO 10 IP=1,NXE                                                            
      READ(5,*)EA,EI,(COORD(IP,I),I=1,4),(G(I),I=1,IDOF)                        
      PROP(IP,1)=EA                                                             
      PROP(IP,2)=EI                                                             
      CALL BMCOL2(KM,EA,EI,IP,COORD,IPROP)                                      
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
C       RETRIEVE ELEMENT END FORCES AND MOMENTS                                 
C                                                                               
      DO 30 IP=1,NXE                                                            
      EA=PROP(IP,1)                                                             
      EI=PROP(IP,2)                                                             
      CALL BMCOL2(KM,EA,EI,IP,COORD,IPROP)                                      
      DO 40 I=1,IDOF                                                            
      G(I)=STOREG(IP,I)                                                         
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   40 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
      CALL LOC2F(LOCAL,ACTION,IP,COORD,IPROP)                                   
      CALL PRINTV(ACTION,IDOF)                                                  
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
