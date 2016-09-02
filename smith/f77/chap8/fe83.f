      PROGRAM P83                                                               
C                                                                               
C      PROGRAM 8.3 SOLUTION OF THE CONDUCTION EQUATION                          
C      OVER A RECTANGULAR AREA USING 4-NODED QUADRILATERALS                     
C      ELEMENT BY ELEMENT PRODUCT INTEGRATION IN TIME                           
C      TWO PASS ALGORITHM                                                       
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IEL=100,ILOADS=150,INF=70)                                      
C                                                                               
      REAL JAC(2,2),JAC1(2,2),KAY(2,2),SAMP(3,2),FTF(4,4),DTKD(4,4),            
     +COORD(4,2),DER(2,4),DERIV(2,4),DERIVT(4,2),KDERIV(2,4),FUN(4),            
     +PM(4,4),KP(4,4),LOADS(ILOADS),MASS(4),GLOBMA(ILOADS),                     
     +STKP(4,4,IEL),STKA(4,4,IEL)                                               
      INTEGER G(4),NF(INF,1)                                                    
      DATA IT,IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV/7*2/,ISAMP/3/                  
      DATA IDTKD,IKP,ICOORD,IDERVT,NOD,IFTF,IPM/7*4/,NODOF/1/                   
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,NN,NR,NGP,AA,BB,PERMX,PERMY                            
      READ(5,*)DTIM,ISTEP,THETA,NPRI,NRES                                       
      CALL NULVEC(GLOBMA,NN)                                                    
      CALL NULL(KAY,IKAY,IT,IT)                                                 
      KAY(1,1)=PERMX                                                            
      KAY(2,2)=PERMY                                                            
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 10 I=1,NN                                                              
   10 NF(I,1)=I                                                                 
C                                                                               
C      ELEMENT INTEGRATION AND LUMPED MASS ASSEMBLY                             
C                                                                               
      NM=0                                                                      
      DO 20 IP=1,NXE                                                            
      DO 20 IQ=1,NYE                                                            
      NM=NM+1                                                                   
      CALL GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KP,IKP,NOD,NOD)                                                 
      CALL NULL(PM,IPM,NOD,NOD)                                                 
      DO 30 I=1,NGP                                                             
      DO 30 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL MATMUL(KAY,IKAY,DERIV,IDERIV,KDERIV,IKDERV,IT,IT,NOD)                
      CALL MATRAN(DERIVT,IDERVT,DERIV,IDERIV,IT,NOD)                            
      CALL MATMUL(DERIVT,IDERVT,KDERIV,IKDERV,DTKD,IDTKD,NOD,IT,NOD)            
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      DO 40 K=1,NOD                                                             
      DO 40 L=1,NOD                                                             
      FTF(K,L)=FUN(K)*FUN(L)*QUOT                                               
   40 DTKD(K,L)=DTKD(K,L)*QUOT                                                  
      CALL MATADD(PM,IPM,FTF,IFTF,NOD,NOD)                                      
   30 CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
      DO 50 I=1,NOD                                                             
      DO 50 J=1,NOD                                                             
   50 STKP(I,J,NM)=KP(I,J)                                                      
      DO 60 I=1,NOD                                                             
      QUOT=0.                                                                   
      DO 70 J=1,NOD                                                             
   70 QUOT=QUOT+PM(I,J)                                                         
   60 GLOBMA(G(I))=GLOBMA(G(I))+QUOT                                            
   20 CONTINUE                                                                  
C                                                                               
C      CALCULATE AND STORE ELEMENT A AND B MATRICES                             
C                                                                               
      NM=0                                                                      
      DO 80 IP=1,NXE                                                            
      DO 80 IQ=1,NYE                                                            
      NM=NM+1                                                                   
      CALL GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      DO 90 I=1,NOD                                                             
      DO 90 J=1,NOD                                                             
      KP(I,J)=-STKP(I,J,NM)*(1.-THETA)*DTIM*.5                                  
   90 FTF(I,J)=STKP(I,J,NM)*THETA*DTIM*.5                                       
      DO 100 I=1,NOD                                                            
      FTF(I,I)=FTF(I,I)+GLOBMA(G(I))                                            
  100 KP(I,I)=KP(I,I)+GLOBMA(G(I))                                              
      CALL MATINV(FTF,IFTF,NOD)                                                 
      CALL MATMUL(FTF,IFTF,KP,IKP,PM,IPM,NOD,NOD,NOD)                           
      DO 110 I=1,NOD                                                            
      DO 110 J=1,NOD                                                            
  110 STKP(I,J,NM)=PM(I,J)                                                      
   80 CONTINUE                                                                  
C                                                                               
C      TAKE ACCOUNT OF INITIAL AND BOUNDARY CONDITIONS                          
C                                                                               
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      READ(5,*)VAL0                                                             
      DO 120 I=1,N                                                              
  120 LOADS(I)=VAL0                                                             
C                                                                               
C      TIME STEPPING LOOP                                                       
C                                                                               
      DO 130 J=1,ISTEP                                                          
      TIME=DTIM*J                                                               
      DO 140 IPASS=1,2                                                          
      IF(IPASS.EQ.1)THEN                                                        
      ILP=1                                                                     
      IHP=NXE                                                                   
      INCP=1                                                                    
      ILQ=1                                                                     
      IHQ=NYE                                                                   
      NM=0                                                                      
      ELSE                                                                      
      ILP=NXE                                                                   
      IHP=1                                                                     
      INCP=-1                                                                   
      ILQ=NYE                                                                   
      IHQ=1                                                                     
      NM=NXE*NYE+1                                                              
      END IF                                                                    
      DO 150 IP=ILP,IHP,INCP                                                    
      DO 150 IQ=ILQ,IHQ,INCP                                                    
      NM=NM+INCP                                                                
      CALL GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      DO 160 K=1,NOD                                                            
      DO 160 L=1,NOD                                                            
  160 FTF(K,L)=STKP(K,L,NM)                                                     
      DO 170 K=1,NOD                                                            
      IF(G(K).EQ.0)MASS(K)=0.                                                   
  170 IF(G(K).NE.0)MASS(K)=LOADS(G(K))                                          
      CALL MVMULT(FTF,IFTF,MASS,NOD,NOD,FUN)                                    
      DO 180 K=1,NOD                                                            
  180 IF(G(K).NE.0)LOADS(G(K))=FUN(K)                                           
  150 CONTINUE                                                                  
  140 CONTINUE                                                                  
      IF(J/NPRI*NPRI.EQ.J)WRITE(6,'(2E12.4)')TIME,LOADS(NRES)                   
  130 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
