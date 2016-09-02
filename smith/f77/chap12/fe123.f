      PROGRAM P123                                                              
C                                                                               
C      PROGRAM 12.3 PILE DRIVABILITY BY THE WAVE EQUATION                       
C      2-NODE LINE ELEMENTS, NEWMARK IMPLICIT INTEGRATION                       
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=1000,ILOADS=200,INF=101)                                    
C                                                                               
      REAL E(INF),RHO(INF),SX(INF),ELL(INF),CSA(INF),KM(2,2),EMM(2,2),          
     +RU(INF),QU(INF),JJ(INF),F(INF),DF(INF),KV(IKV),F1(IKV),MM(IKV),           
     +LOADS(ILOADS),X0(ILOADS),D1X0(ILOADS),D2X0(ILOADS),ELD(2),                
     +BDYLDS(ILOADS),X1(ILOADS),D1X1(ILOADS),D2X1(ILOADS)                       
      INTEGER NF(INF,1),G(2)                                                    
      DATA IKM,IEMM,IDOF/3*2/,NODOF,IW/2*1/                                     
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,N,NN,NR,ALPHA,BETA,ITS,DTIM,ISTEP,THETA,ITYPE,               
     +         RAMVEL,STATLD                                                    
      IR=N*(IW+1)                                                               
      CALL NULVEC(LOADS,N)                                                      
      CALL NULVEC(BDYLDS,N)                                                     
      CALL NULVEC(SX,NXE)                                                       
      CALL NULVEC(F,NN)                                                         
      CALL NULVEC(DF,NN)                                                        
      CALL NULVEC(KV,IR)                                                        
      CALL NULVEC(MM,IR)                                                        
      READ(5,*)(E(I),I=1,NXE)                                                   
      READ(5,*)(RHO(I),I=1,NXE)                                                 
      READ(5,*)(CSA(I),I=1,NXE)                                                 
      READ(5,*)(ELL(I),I=1,NXE)                                                 
      READ(5,*)(RU(I),I=1,NN)                                                   
      READ(5,*)(QU(I),I=1,NN)                                                   
      READ(5,*)(JJ(I),I=1,NN)                                                   
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
C                                                                               
C      ELEMENT MASS AND STIFFNESS ASSEMBLY                                      
C                                                                               
      DO 10 IP=1,NXE                                                            
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      CALL AXIKM(KM,CSA(IP),E(IP),ELL(IP))                                      
      IF(ITYPE.EQ.1)GO TO 20                                                    
      EMM(1,1)=RHO(IP)*CSA(IP)*ELL(IP)/3.                                       
      EMM(2,2)=EMM(1,1)                                                         
      EMM(2,1)=EMM(1,1)*.5                                                      
      EMM(1,2)=EMM(2,1)                                                         
      GOTO 30                                                                   
   20 EMM(1,1)=RHO(IP)*CSA(IP)*ELL(IP)*.5                                       
      EMM(2,2)=EMM(1,1)                                                         
      EMM(2,1)=.0                                                               
      EMM(1,2)=EMM(2,1)                                                         
   30 CALL FORMKV(KV,KM,IKM,G,N,IDOF)                                           
   10 CALL FORMKV(MM,EMM,IEMM,G,N,IDOF)                                         
      CALL VECCOP(KV,F1,IR)                                                     
C                                                                               
C      PRELIMINARY STATIC SOLUTION                                              
C                                                                               
      F1(N)=F1(N)+1.E6                                                          
      CALL BANRED(F1,N,IW)                                                      
      LOADS(1)=STATLD                                                           
      CALL BACSUB(F1,LOADS,N,IW)                                                
      CALL PRINTV(LOADS,N)                                                      
      CALL NULVEC(LOADS,N)                                                      
C                                                                               
C      INITIAL CONDITIONS                                                       
C                                                                               
      CALL NULVEC(X0,N)                                                         
      CALL NULVEC(D1X0,N)                                                       
      CALL NULVEC(D2X0,N)                                                       
      CALL NULVEC(D1X1,N)                                                       
      D1X0(2)=RAMVEL                                                            
      D1X0(1)=D1X0(2)                                                           
      C1=(1.-THETA)*DTIM                                                        
      C2=BETA-C1                                                                
      C3=ALPHA+1./(THETA*DTIM)                                                  
      C4=BETA+THETA*DTIM                                                        
C                                                                               
C      NEWMARK THETA METHOD FOR TIME STEPPING                                   
C                                                                               
      TIM=0.0                                                                   
      DO 40 J=1,ISTEP                                                           
      TIM=TIM+DTIM                                                              
      DO 50 I=1,IR                                                              
   50 F1(I)=C3*MM(I)+C4*KV(I)                                                   
      DO 60 I=1,NN                                                              
      IF(QU(I).LT.1.E-6)GOTO 60                                                 
      F1(I)=F1(I)+C4*RU(I)*(1.+JJ(I)*D1X1(I))/QU(I)                             
   60 CONTINUE                                                                  
      CALL BANRED(F1,N,IW)                                                      
C                                                                               
C      REDISTRIBUTE EXCESS SPRING FORCES                                        
C                                                                               
      DO 70 K=1,ITS                                                             
      DO 80 I=1,N                                                               
      LOADS(I)=DTIM*BDYLDS(I)                                                   
      BDYLDS(I)=0.0                                                             
   80 X1(I)=C3*X0(I)+D1X0(I)/THETA                                              
      CALL LINMUL(MM,X1,D1X1,N,IW)                                              
      CALL VECADD(D1X1,LOADS,D1X1,N)                                            
      DO 90 I=1,N                                                               
   90 LOADS(I)=C2*X0(I)                                                         
      CALL LINMUL(KV,LOADS,X1,N,IW)                                             
      CALL VECADD(X1,D1X1,X1,N)                                                 
      CALL BACSUB(F1,X1,N,IW)                                                   
      DO 100 I=1,N                                                              
      D1X1(I)=(X1(I)-X0(I))/(THETA*DTIM)-D1X0(I)*(1.-THETA)/THETA               
  100 D2X1(I)=(D1X1(I)-D1X0(I))/(THETA*DTIM)-D2X0(I)*(1.-THETA)/THETA           
      DO 110 I=1,NN                                                             
      IF(QU(I).LT.1.E-6)GO TO 110                                               
      W1=1.+JJ(I)*D1X1(I)                                                       
      DF(I)=-(X1(I)-X0(I))*W1*RU(I)/QU(I)                                       
      IF(ABS(F(I)+DF(I)).LT.W1*RU(I))GO TO 120                                  
      IF(DF(I).GT..0)BDYLDS(I)=-F(I)-DF(I)+W1*RU(I)                             
      IF(DF(I).LE..0)BDYLDS(I)=-F(I)-DF(I)-W1*RU(I)                             
  120 DF(I)=-(X1(I)-X0(I))*W1*RU(I)/QU(I)+BDYLDS(I)                             
  110 CONTINUE                                                                  
   70 CONTINUE                                                                  
      DO 130 IP=1,NXE                                                           
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      DO 140 M=1,2                                                              
      IF(G(M).EQ.0)ELD(M)=.0                                                    
  140 IF(G(M).NE.0)ELD(M)=X1(G(M))                                              
      EPS=(ELD(2)-ELD(1))/ELL(IP)                                               
      SIGMA=E(IP)*EPS                                                           
  130 SX(IP)=SIGMA                                                              
      CALL VECADD(F,DF,F,NN)                                                    
      WRITE(6,1000)TIM                                                          
      CALL PRINTV(X1,N)                                                         
      CALL PRINTV(D1X1,N)                                                       
      CALL PRINTV(F,NN)                                                         
      CALL PRINTV(SX,NXE)                                                       
C                                                                               
C      UPDATE DISPLACEMENTS,VELOCITIES,ACCELERATIONS                            
C                                                                               
      CALL VECCOP(X1,X0,N)                                                      
      CALL VECCOP(D1X1,D1X0,N)                                                  
      CALL VECCOP(D2X1,D2X0,N)                                                  
   40 CONTINUE                                                                  
 1000 FORMAT(F8.5)                                                              
      STOP                                                                      
      END                                                                       
