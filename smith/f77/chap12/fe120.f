      PROGRAM P120                                                              
C                                                                               
C      PROGRAM 12.0 T-Z ANALYSIS OF AXIALLY LOADED PILES                        
C      USING 2-NODE LINE ELEMENTS                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=200,ILOADS=100,INO=10,INF=100,IQINC=50)                     
C                                                                               
      REAL KV(IKV),LOADS(ILOADS),DISPS(ILOADS),BDYLDS(ILOADS),                  
     +OLDIS(ILOADS),RU(INF),QU(INF),F(INF),DF(INF),SX(INF),                     
     +QINC(IQINC),VAL(INO),ELD(2),KM(2,2)                                       
      INTEGER G(2),NF(INF,1),NO(INO)                                            
      DATA IKM,IDOF/2*2/,IW,NODOF/2*1/                                          
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,N,NN,NR,CSA,E,ELL,ITS                                        
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      READ(5,*)(RU(I),I=1,NXE+1)                                                
      READ(5,*)(QU(I),I=1,NXE+1)                                                
      READ(5,*)NL,(NO(I),I=1,NL)                                                
      READ(5,*)INCS,(QINC(I),I=1,INCS)                                          
      IR=N*(IW+1)                                                               
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(OLDIS,N)                                                      
      CALL NULVEC(DISPS,N)                                                      
      CALL NULVEC(SX,NXE)                                                       
      CALL NULVEC(F,NN)                                                         
      CALL AXIKM(KM,CSA,E,ELL)                                                  
C                                                                               
C      LOAD INCREMENT LOOP                                                      
C                                                                               
      DO 10 L=1,INCS                                                            
      CALL NULVEC(KV,IR)                                                        
C                                                                               
C      ASSEMBLE GLOBAL STIFFNESS MATRIX                                         
C                                                                               
      DO 20 IP=1,NXE                                                            
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
   20 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
C                                                                               
C      ADD T-Z SPRINGS                                                          
C                                                                               
      DO 30 I=1,NN                                                              
      IF(QU(I).LT.1.E-6)GOTO 30                                                 
      IF(ABS(DISPS(I)).LE..75*QU(I))KV(I)=KV(I)+RU(I)/QU(I)                     
      IF(ABS(DISPS(I)).GT..75*QU(I))KV(I)=KV(I)+.1*RU(I)/QU(I)                  
   30 CONTINUE                                                                  
C                                                                               
C      REDUCE EQUATIONS                                                         
C                                                                               
      CALL BANRED(KV,N,IW)                                                      
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
      ITERS=0                                                                   
   40 ITERS=ITERS+1                                                             
      CALL NULVEC(LOADS,N)                                                      
      DO 50 I=1,NL                                                              
   50 LOADS(NO(I))=QINC(L)                                                      
      CALL VECADD(LOADS,BDYLDS,LOADS,N)                                         
      CALL NULVEC(BDYLDS,N)                                                     
C                                                                               
C      SOLVE EQUATIONS                                                          
C                                                                               
      CALL BACSUB(KV,LOADS,N,IW)                                                
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      IF(ITS.EQ.1)GOTO 60                                                       
      CALL CHECON(LOADS,OLDIS,N,0.00001,ICON)                                   
      IF(ITERS.EQ.1)ICON=0                                                      
C                                                                               
C      REDISTRIBUTE EXCESS SPRING FORCES                                        
C                                                                               
      DO 70 I=1,NN                                                              
      IF(QU(I).LT.1.E-6)GOTO 70                                                 
      DF(I)=LOADS(I)*RU(I)/QU(I)                                                
      IF(ABS(DISPS(I)).GT..75*QU(I))DF(I)=.1*DF(I)                              
      IF(ABS(F(I)+DF(I)).GT.RU(I))THEN                                          
      BDYLDS(I)=F(I)+DF(I)+RU(I)                                                
      DF(I)=-RU(I)-F(I)                                                         
      END IF                                                                    
   70 CONTINUE                                                                  
      IF(ICON.EQ.0.AND.ITERS.NE.ITS)GOTO 40                                     
C                                                                               
C      COMPUTE ELEMENT STRESSES AND STRAINS                                     
C                                                                               
   60 DO 80 IP=1,NXE                                                            
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      DO 90 J=1,2                                                               
      IF(G(J).EQ.0)ELD(J)=0.                                                    
   90 IF(G(J).NE.0)ELD(J)=LOADS(G(J))                                           
      EPS=(ELD(2)-ELD(1))/ELL                                                   
      SIGMA=E*EPS                                                               
      SX(IP)=SX(IP)+SIGMA                                                       
   80 CONTINUE                                                                  
      CALL VECADD(DISPS,LOADS,DISPS,N)                                          
C                                                                               
C      CHECK SUM OF SPRING FORCES                                               
C                                                                               
      SUM=0.                                                                    
      DO 100 I=1,NN                                                             
      F(I)=F(I)+DF(I)                                                           
  100 SUM=SUM+F(I)                                                              
      WRITE(6,1000)SUM,DISPS(1),ITERS                                           
      IF(ITERS.EQ.ITS)GOTO 110                                                  
   10 CONTINUE                                                                  
  110 CONTINUE                                                                  
 1000 FORMAT(2E12.4,I10)                                                        
      STOP                                                                      
      END                                                                       
