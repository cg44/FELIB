      PROGRAM P47                                                               
C                                                                               
C      PROGRAM 4.7 PLASTIC HINGE ANALYSIS OF RIGID JOINTED FRAMES               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=400,ILOADS=100,IPROP=20,INO=20,ISTEP=20)                    
C                                                                               
      REAL OLDSPS(ILOADS),LOADS(ILOADS),BDYLDS(ILOADS),ELDTOT(ILOADS),          
     +KV(IKV),HOLDR(IPROP,6),PROP(IPROP,3),VAL(INO),                            
     +ACTION(6),REACT(6),ELD(6),KM(6,6),COORD(IPROP,4),DLOAD(ISTEP)             
      INTEGER G(6),STOREG(IPROP,6),NO(INO)                                      
      DATA IKM,IDOF/2*6/,NODOF/3/                                               
C                                                                               
C      INPUT SECTION                                                            
C                                                                               
      READ(5,*)NXE,N,IW,ITS,TOL                                                 
      IR=N*(IW+1)                                                               
      CALL NULVEC(KV,IR)                                                        
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(ELDTOT,N)                                                     
      CALL NULL(HOLDR,IPROP,NXE,IDOF)                                           
      READ(5,*)EA,EI                                                            
C                                                                               
C      GLOBAL STIFFNESS MATRIX ASSEMBLY                                         
C                                                                               
      DO 10 IP=1,NXE                                                            
      READ(5,*)PM,(COORD(IP,I),I=1,4),(G(I),I=1,6)                              
      PROP(IP,1)=EA                                                             
      PROP(IP,2)=EI                                                             
      PROP(IP,3)=PM                                                             
      DO 20 I=1,IDOF                                                            
   20 STOREG(IP,I)=G(I)                                                         
      CALL BMCOL2(KM,EA,EI,IP,COORD,IPROP)                                      
   10 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      NODAL LOADING AND STIFFNESS MATRIX REDUCTION                             
C                                                                               
      READ(5,*)NL,(NO(I),VAL(I),I=1,NL)                                         
      READ(5,*)INCS,(DLOAD(I),I=1,INCS)                                         
      CALL BANRED(KV,N,IW)                                                      
C                                                                               
C      LOAD INCREMENT LOOP                                                      
C                                                                               
      TOTLO=0.                                                                  
      DO 30 IY=1,INCS                                                           
      TOTLO=TOTLO+DLOAD(IY)                                                     
      CALL NULVEC(OLDSPS,N)                                                     
      ITERS=0                                                                   
   40 ITERS=ITERS+1                                                             
      CALL NULVEC(LOADS,N)                                                      
      DO 50 I=1,NL                                                              
   50 LOADS(NO(I))=DLOAD(IY)*VAL(I)                                             
      CALL VECADD(LOADS,BDYLDS,LOADS,N)                                         
      CALL NULVEC(BDYLDS,N)                                                     
      CALL BACSUB(KV,LOADS,N,IW)                                                
C                                                                               
C     CHECK CONVERGENCE                                                         
C                                                                               
      CALL CHECON(LOADS,OLDSPS,N,TOL,ICON)                                      
C                                                                               
C     INSPECT MOMENTS IN ALL ELEMENTS                                           
C                                                                               
      DO 60 IP=1,NXE                                                            
      EA=PROP(IP,1)                                                             
      EI=PROP(IP,2)                                                             
      PM=PROP(IP,3)                                                             
      CALL BMCOL2(KM,EA,EI,IP,COORD,IPROP)                                      
      DO 70 I=1,IDOF                                                            
      G(I)=STOREG(IP,I)                                                         
      IF(G(I).EQ.0)ELD(I)=0.                                                    
   70 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(KM,IKM,ELD,IDOF,IDOF,ACTION)                                  
      CALL NULVEC(REACT,IDOF)                                                   
      IF(ITS.EQ.1)GOTO 80                                                       
C                                                                               
C     IF PM EXCEEDED GENERATE SELF-EQUILIBRATING VECTOR 'REACT'                 
C     TO SUBTRACT FROM LOADS VECTOR                                             
C                                                                               
      CALL HING2(IP,HOLDR,COORD,IPROP,ACTION,REACT,PM)                          
      DO 90 M=1,IDOF                                                            
      IF(G(M).EQ.0)GOTO 90                                                      
      BDYLDS(G(M))=BDYLDS(G(M))-REACT(M)                                        
   90 CONTINUE                                                                  
C                                                                               
C     AT CONVERGENCE UPDATE ELEMENT REACTIONS, PRINT RESULTS,                   
C     AND MOVE ON TO NEXT LOAD INCREMENT                                        
C                                                                               
   80 IF(ITERS.NE.ITS.AND.ICON.NE.1)GOTO 60                                     
      DO 100 M=1,IDOF                                                           
  100 HOLDR(IP,M)=HOLDR(IP,M)+REACT(M)+ACTION(M)                                
   60 CONTINUE                                                                  
      IF(ITERS.NE.ITS.AND.ICON.NE.1)GOTO 40                                     
      CALL VECADD(LOADS,ELDTOT,ELDTOT,N)                                        
      WRITE(6,'(/,E12.4)')TOTLO                                                 
      WRITE(6,'(10E12.4)')(ELDTOT(NO(I)),I=1,NL)                                
      WRITE(6,'(10I10)')ITERS                                                   
   30 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
