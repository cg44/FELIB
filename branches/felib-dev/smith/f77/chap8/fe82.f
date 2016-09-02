      PROGRAM P82                                                               
C                                                                               
C      PROGRAM 8.2 SOLUTION OF THE CONDUCTION EQUATION                          
C      OVER A RECTANGULAR AREA USING 4-NODED QUADRILATERALS                     
C      SIMPLE EXPLICIT INTEGRATION IN TIME USING                                
C      ELEMENT BY ELEMENT SUMMATION                                             
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IEL=30,ILOADS=50,INF=70)                                        
C                                                                               
      REAL JAC(2,2),JAC1(2,2),KAY(2,2),SAMP(3,2),FTF(4,4),DTKD(4,4),            
     +COORD(4,2),DER(2,4),DERIV(2,4),DERIVT(4,2),KDERIV(2,4),FUN(4),            
     +PM(4,4),KP(4,4),NEWLO(ILOADS),LOADS(ILOADS),ELD(4),                       
     +STPM(4,4,IEL),GLOBMA(INF),MASS(4)                                         
      INTEGER G(4),NF(INF,1)                                                    
      DATA IT,IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV/7*2/,ISAMP/3/                  
      DATA IDTKD,IKP,ICOORD,IDERVT,NOD,IFTF,IPM/7*4/,NODOF/1/                   
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,NN,NR,NGP,AA,BB,PERMX,PERMY                            
      READ(5,*)DTIM,ISTEP,NPRI,NRES                                             
      CALL NULVEC(GLOBMA,NN)                                                    
      CALL NULL(KAY,IKAY,IT,IT)                                                 
      KAY(1,1)=PERMX                                                            
      KAY(2,2)=PERMY                                                            
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 10 I=1,NN                                                              
   10 NF(I,1)=I                                                                 
C                                                                               
C      ELEMENT INTEGRATION AND STORAGE                                          
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
   40 DTKD(K,L)=DTKD(K,L)*QUOT*DTIM                                             
      CALL MATADD(PM,IPM,FTF,IFTF,NOD,NOD)                                      
   30 CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
      DO 50 I=1,NOD                                                             
      QUOT=0.                                                                   
      DO 60 J=1,NOD                                                             
   60 QUOT=QUOT+PM(I,J)                                                         
      MASS(I)=QUOT                                                              
   50 GLOBMA(G(I))=GLOBMA(G(I))+MASS(I)                                         
      CALL NULL(PM,IPM,NOD,NOD)                                                 
      DO 70 I=1,NOD                                                             
   70 PM(I,I)=MASS(I)                                                           
      DO 80 I=1,NOD                                                             
      DO 80 J=1,NOD                                                             
   80 STPM(I,J,NM)=PM(I,J)-KP(I,J)                                              
   20 CONTINUE                                                                  
C                                                                               
C      INITIAL CONDITIONS AND REARRANGE GLOBMA INTO M-1                         
C                                                                               
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      J=0                                                                       
      DO 90 I=1,NN                                                              
      IF(NF(I,1).EQ.0)GO TO 90                                                  
      J=J+1                                                                     
      GLOBMA(J)=1./GLOBMA(I)                                                    
   90 CONTINUE                                                                  
      READ(5,*)VAL0                                                             
      DO 100 I=1,N                                                              
  100 LOADS(I)=VAL0                                                             
C                                                                               
C      TIME STEPPING RECURSION                                                  
C                                                                               
      DO 110 J=1,ISTEP                                                          
      CALL NULVEC(NEWLO,N)                                                      
      TIME=J*DTIM                                                               
      NM=0                                                                      
      DO 120 IP=1,NXE                                                           
      DO 120 IQ=1,NYE                                                           
      NM=NM+1                                                                   
      CALL GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      DO 130 I=1,NOD                                                            
      DO 130 K=1,NOD                                                            
  130 PM(I,K)=STPM(I,K,NM)                                                      
      DO 140 I=1,NOD                                                            
      IF(G(I).EQ.0)ELD(I)=0.                                                    
  140 IF(G(I).NE.0)ELD(I)=LOADS(G(I))                                           
      CALL MVMULT(PM,IPM,ELD,NOD,NOD,FUN)                                       
      DO 150 I=1,NOD                                                            
  150 IF(G(I).NE.0)NEWLO(G(I))=NEWLO(G(I))+FUN(I)                               
  120 CONTINUE                                                                  
      DO 160 I=1,N                                                              
  160 LOADS(I)=NEWLO(I)*GLOBMA(I)                                               
      IF(J/NPRI*NPRI.EQ.J)WRITE(6,'(2E12.4)')TIME,LOADS(NRES)                   
  110 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
