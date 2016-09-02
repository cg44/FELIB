      PROGRAM P81                                                               
C                                                                               
C      PROGRAM 8.1 SOLUTION OF THE CONDUCTION EQUATION                          
C      OVER A CYLINDER USING 4-NODED QUADRILATERALS                             
C      IMPLICIT INTEGRATION IN TIME BY 'THETA' METHOD                           
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=1000,ILOADS=150,INF=70)                                     
C                                                                               
      REAL JAC(2,2),JAC1(2,2),KAY(2,2),SAMP(3,2),FTF(4,4),DTKD(4,4),            
     +COORD(4,2),DER(2,4),DERIV(2,4),DERIVT(4,2),KDERIV(2,4),FUN(4),            
     +PM(4,4),KP(4,4),NEWLO(ILOADS),LOADS(ILOADS),BP(IKV),BK(IKV)               
      INTEGER G(4),NF(INF,1)                                                    
      DATA IT,IJAC,IJAC1,IKAY,IDER,IDERIV,IKDERV/7*2/,ISAMP/3/                  
      DATA IDTKD,IKP,ICOORD,IDERVT,NOD,IFTF,IPM/7*4/,NODOF/1/                   
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NR,NGP,AA,BB,PERMX,PERMY                         
      READ(5,*)DTIM,ISTEP,THETA,NPRI,NRES                                       
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IR=N*(IW+1)                                                               
      CALL NULVEC(BK,IR)                                                        
      CALL NULVEC(BP,IR)                                                        
      CALL NULL(KAY,IKAY,IT,IT)                                                 
      KAY(1,1)=PERMX                                                            
      KAY(2,2)=PERMY                                                            
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT INTEGRATION AND ASSEMBLY                                         
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KP,IKP,NOD,NOD)                                                 
      CALL NULL(PM,IPM,NOD,NOD)                                                 
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL MATMUL(KAY,IKAY,DERIV,IDERIV,KDERIV,IKDERV,IT,IT,NOD)                
      CALL MATRAN(DERIVT,IDERVT,DERIV,IDERIV,IT,NOD)                            
      CALL MATMUL(DERIVT,IDERVT,KDERIV,IKDERV,DTKD,IDTKD,NOD,IT,NOD)            
      SUM=0.                                                                    
      DO 30 K=1,NOD                                                             
   30 SUM=SUM+FUN(K)*COORD(K,1)                                                 
      QUOT=SUM*DET*SAMP(I,2)*SAMP(J,2)                                          
      DO 40 K=1,NOD                                                             
      DO 40 L=1,NOD                                                             
      FTF(K,L)=FUN(K)*FUN(L)*QUOT                                               
   40 DTKD(K,L)=DTKD(K,L)*QUOT*THETA*DTIM                                       
      CALL MATADD(PM,IPM,FTF,IFTF,NOD,NOD)                                      
   20 CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
      CALL FORMKV(BP,PM,IPM,G,N,NOD)                                            
   10 CALL FORMKV(BK,KP,IKP,G,N,NOD)                                            
C                                                                               
C      REDUCTION OF LEFT HAND SIDE                                              
C                                                                               
      DO 50 I=1,IR                                                              
      BP(I)=BP(I)+BK(I)                                                         
   50 BK(I)=BP(I)-BK(I)/THETA                                                   
      CALL BANRED(BP,N,IW)                                                      
      READ(5,*)VAL0                                                             
      DO 60 I=1,N                                                               
   60 LOADS(I)=VAL0                                                             
C                                                                               
C      TIME STEPPING RECURSION                                                  
C                                                                               
      DO 70 J=1,ISTEP                                                           
      TIME=J*DTIM                                                               
      CALL LINMUL(BK,LOADS,NEWLO,N,IW)                                          
      CALL BACSUB(BP,NEWLO,N,IW)                                                
      CALL VECCOP(NEWLO,LOADS,N)                                                
      IF(J/NPRI*NPRI.EQ.J)WRITE(6,'(2E12.4)')TIME,LOADS(NRES)                   
   70 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
