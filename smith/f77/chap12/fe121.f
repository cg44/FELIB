      PROGRAM P121                                                              
C                                                                               
C      PROGRAM 12.1 P-Y ANALYSIS OF LATERALLY LOADED PILES                      
C      USING 2-NODE BEAM ELEMENTS                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=500,ILOADS=200,INO=10,ISEG=12,INF=100)                      
C                                                                               
      REAL KV(IKV),LOADS(ILOADS),BDYLDS(ILOADS),DISPS(ILOADS),                  
     +OLDIS(ILOADS),RU(INF,ISEG),QU(INF,ISEG),F(INF),DF(INF),                   
     +MOM(INF),VAL(INO),STORE(INO),KM(4,4),KP(4,4),ACTION(4),ELD(4)             
      INTEGER G(4),NF(INF,2),NO(INO)                                            
      DATA IKM,IDOF/2*4/,NODOF/2/,IW/3/                                         
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,N,NN,NR,PA,NP,EI,ELL,INCS,ITS                                
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      READ(5,*)((RU(I,J),J=1,NP+1),I=1,NN)                                      
      READ(5,*)((QU(I,J),J=1,NP+1),I=1,NN)                                      
      READ(5,*)NL,(NO(I),VAL(I),I=1,NL)                                         
      IR=N*(IW+1)                                                               
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(OLDIS,N)                                                      
      CALL NULVEC(DISPS,N)                                                      
      CALL NULVEC(F,NN)                                                         
      CALL NULVEC(MOM,NN)                                                       
      CALL BEAMKM(KM,EI,ELL)                                                    
C                                                                               
C      LOAD INCREMENT LOOP                                                      
C                                                                               
      DO 10 L=1,INCS                                                            
      CALL BEAMKP(KP,ELL)                                                       
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C      MODIFY ELEMENT STIFFNESS FOR AXIAL LOADING                               
C                                                                               
      DO 20 I=1,IDOF                                                            
      DO 20 J=1,IDOF                                                            
   20 KP(I,J)=KM(I,J)-L*PA*KP(I,J)                                              
C                                                                               
C      ASSEMBLE GLOBAL STIFFNESS MATRIX                                         
C                                                                               
      DO 30 IP=1,NXE                                                            
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
   30 CALL FORMKV(KV,KP,IKM,G,N,IDOF)                                           
C                                                                               
C      ADD P-Y SPRINGS                                                          
C                                                                               
      DO 40 I=1,NN                                                              
      II=2*I-1                                                                  
   40 KV(II)=KV(II)+RU(I,2)/QU(I,2)                                             
C                                                                               
C      ADD 'BIG SPRINGS' FOR PRESCRIBED DISPLACEMENTS                           
C      AND REDUCE EQUATIONS                                                     
C                                                                               
      DO 50 I=1,NL                                                              
      KV(NO(I))=KV(NO(I))+1.E20                                                 
   50 STORE(I)=KV(NO(I))                                                        
      CALL BANRED(KV,N,IW)                                                      
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
      ITERS=0                                                                   
   60 ITERS=ITERS+1                                                             
      CALL NULVEC(LOADS,N)                                                      
      DO 70 I=1,NL                                                              
   70 LOADS(NO(I))=STORE(I)*VAL(I)                                              
      CALL VECADD(LOADS,BDYLDS,LOADS,N)                                         
      CALL NULVEC(BDYLDS,N)                                                     
C                                                                               
C      SOLVE EQUATIONS                                                          
C                                                                               
      CALL BACSUB(KV,LOADS,N,IW)                                                
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(LOADS,OLDIS,N,.00001,ICON)                                    
      IF(ITERS.EQ.1)ICON=0                                                      
C                                                                               
C      REDISTRIBUTE EXCESS SPRING FORCES                                        
C                                                                               
      DO 80 I=1,NN                                                              
      II=2*I-1                                                                  
      DF(I)=-LOADS(II)*RU(I,2)/QU(I,2)                                          
      J=1                                                                       
      DO 90 IP=2,NP+1                                                           
   90 IF(ABS(DISPS(II)+LOADS(II)).GT.QU(I,IP))J=J+1                             
      IF(J.EQ.1.AND.DISPS(II)*LOADS(II).GT.0.)GOTO 100                          
      FORCE=(ABS(DISPS(II)+LOADS(II))-QU(I,J))*(RU(I,J+1)-RU(I,J))              
     +/(QU(I,J+1)-QU(I,J))+RU(I,J)                                              
      IF(DISPS(II).LT.0.)BDYLDS(II)=-F(I)-DF(I)+FORCE                           
      IF(DISPS(II).GT.0.)BDYLDS(II)=-F(I)-DF(I)-FORCE                           
  100 DF(I)=BDYLDS(II)+DF(I)                                                    
   80 CONTINUE                                                                  
      IF(ICON.EQ.0.AND.ITERS.NE.ITS)GOTO 60                                     
C                                                                               
C      COMPUTE ELEMENT MOMENTS AND SHEARS                                       
C                                                                               
      DO 110 IP=1,NXE                                                           
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      DO 120 M=1,IDOF                                                           
      IF(G(M).EQ.0)ELD(M)=.0                                                    
  120 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      CALL MVMULT(KP,IKM,ELD,IDOF,IDOF,ACTION)                                  
      MOM(IP)=MOM(IP)+ACTION(2)                                                 
      IF(IP.EQ.NXE)MOM(IP+1)=MOM(IP+1)+ACTION(4)                                
  110 CONTINUE                                                                  
      CALL VECADD(DISPS,LOADS,DISPS,N)                                          
      SUM=0.0                                                                   
      DO 130 I=1,NN                                                             
      F(I)=F(I)+DF(I)                                                           
  130 SUM=SUM+F(I)                                                              
      WRITE(6,1000)SUM,DISPS(1),ITERS                                           
      IF(ITERS.EQ.ITS)GOTO 140                                                  
   10 CONTINUE                                                                  
  140 CONTINUE                                                                  
 1000 FORMAT(2E12.4,I10)                                                        
      STOP                                                                      
      END                                                                       
