      PROGRAM P48                                                               
C                                                                               
C      PROGRAM 4.8 STABILITY ANALYSIS OF PLANE FRAMES                           
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,IPROP=20,INO=20,ISTEP=20)                    
C                                                                               
      REAL OLDSPS(ILOADS),ACTION(6),LOCAL(6),ELD(6),KM(6,6),KP(6,6),            
     +LOADS(ILOADS),KV(IKV),PROP(IPROP,2),VAL(INO),COORD(IPROP,4),              
     +DISPS(ILOADS),ELDTOT(ILOADS),AXIF(IPROP),AXIP(IPROP),KCOP(IKV),           
     +DLOAD(ISTEP)                                                              
      INTEGER G(6),STOREG(IPROP,6),NO(INO)                                      
      DATA IKM,IDOF/2*6/,NODOF/3/                                               
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,N,IW,ITS,TOL                                                 
      READ(5,*)((PROP(IP,I),I=1,2),(COORD(IP,I),I=1,4),                         
     +          (STOREG(IP,I),I=1,6),IP=1,NXE)                                  
      READ(5,*)NL,(NO(I),VAL(I),I=1,NL)                                         
      READ(5,*)INCS,(DLOAD(I),I=1,INCS)                                         
      IR=N*(IW+1)                                                               
      CALL NULVEC(AXIP,NXE)                                                     
      CALL NULVEC(AXIF,NXE)                                                     
      CALL NULVEC(ELDTOT,N)                                                     
C                                                                               
C      LOAD INCREMENT LOOP                                                      
C                                                                               
      TOTLO=0.                                                                  
      DO 10 IY=1,INCS                                                           
      TOTLO=TOTLO+DLOAD(IY)                                                     
      CALL NULVEC(LOADS,N)                                                      
      DO 20 I=1,NL                                                              
   20 LOADS(NO(I))=DLOAD(IY)*VAL(I)                                             
      CALL NULVEC(OLDSPS,N)                                                     
      ITERS=0                                                                   
   30 ITERS=ITERS+1                                                             
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C       GLOBAL STIFFNESS MATRIX ASSEMBLY                                        
C                                                                               
      DO 40 IP=1,NXE                                                            
      DO 50 I=1,6                                                               
   50 G(I)=STOREG(IP,I)                                                         
      EA=PROP(IP,1)                                                             
      EI=PROP(IP,2)                                                             
      PAX=AXIF(IP)                                                              
      CALL STAB2D(KM,EA,EI,IP,COORD,IPROP,PAX)                                  
   40 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
      CALL VECCOP(KV,KCOP,IR)                                                   
      CALL KVDET(KCOP,N,IW,DET,KSC)                                             
C                                                                               
C      EQUATION SOLUTION                                                        
C                                                                               
      CALL VECCOP(LOADS,DISPS,N)                                                
      CALL BANRED(KV,N,IW)                                                      
      CALL BACSUB(KV,DISPS,N,IW)                                                
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(DISPS,OLDSPS,N,TOL,ICON)                                      
C                                                                               
C      RETRIEVE ELEMENT END FORCES AND MOMENTS                                  
C                                                                               
      DO 60 IP=1,NXE                                                            
      DO 70 I=1,6                                                               
   70 G(I)=STOREG(IP,I)                                                         
      EA=PROP(IP,1)                                                             
      EI=PROP(IP,2)                                                             
      PAX=AXIF(IP)                                                              
      CALL STAB2D(KM,EA,EI,IP,COORD,IPROP,PAX)                                  
      DO 80 I=1,IDOF                                                            
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   80 IF(G(I).NE.0)ELD(I)=DISPS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
      CALL LOC2F(LOCAL,ACTION,IP,COORD,IPROP)                                   
      AXIF(IP)=AXIP(IP)+LOCAL(4)                                                
   60 CONTINUE                                                                  
      IF(ITERS.NE.ITS.AND.ICON.EQ.0)GOTO 30                                     
C                                                                               
C      AT CONVERGENCE UPDATE DISPLACEMENTS AND AXIAL FORCES                     
C                                                                               
      CALL VECCOP(AXIF,AXIP,NXE)                                                
      CALL VECADD(ELDTOT,DISPS,ELDTOT,N)                                        
      WRITE(6,'(/,E12.4)')TOTLO                                                 
      WRITE(6,'(10E12.4)')(ELDTOT(NO(I)),I=1,NL)                                
      WRITE(6,'(I10,E12.4)')ITERS,DET                                           
   10 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
