      PROGRAM P44                                                               
C                                                                               
C      PROGRAM 4.4 ANALYSIS OF 3-D FRAMES                                       
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,IPROP=20)                                    
C                                                                               
      REAL ACTION(12),LOCAL(12),ELD(12),KM(12,12),                              
     +KV(IKV),LOADS(ILOADS),PROP(IPROP,4),COORD(IPROP,7)                        
      INTEGER G(12),STOREG(IPROP,12)                                            
      DATA IKM,IDOF/2*12/,NODOF/6/                                              
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
      READ(5,*)EA,EIY,EIZ,GJ,(COORD(IP,I),I=1,7),(G(I),I=1,IDOF)                
      PROP(IP,1)=EA                                                             
      PROP(IP,2)=EIY                                                            
      PROP(IP,3)=EIZ                                                            
      PROP(IP,4)=GJ                                                             
      CALL BMCOL3(KM,EA,EIY,EIZ,GJ,IP,COORD,IPROP)                              
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
C      RETRIEVE ELEMENT END FORCES AND MOMENTS                                  
C                                                                               
      DO 30 IP=1,NXE                                                            
      EA=PROP(IP,1)                                                             
      EIY=PROP(IP,2)                                                            
      EIZ=PROP(IP,3)                                                            
      GJ=PROP(IP,4)                                                             
      CALL BMCOL3(KM,EA,EIY,EIZ,GJ,IP,COORD,IPROP)                              
      DO 40 I=1,IDOF                                                            
      G(I)=STOREG(IP,I)                                                         
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   40 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
      CALL LOC3F(LOCAL,ACTION,IP,COORD,IPROP)                                   
      CALL PRINTV(ACTION,IDOF)                                                  
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
