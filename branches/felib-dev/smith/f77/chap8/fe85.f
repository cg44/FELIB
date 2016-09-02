      PROGRAM P85                                                               
C                                                                               
C      PROGRAM 8.5  DIFFUSION CONVECTION EQUATION ON A RECTANGLE                
C      USING 4-NODED QUADRILATERAL ELEMENTS. UNTRANSFORMED                      
C      EQUATION BY GALERKIN'S METHOD; IMPLICIT INTEGRATION                      
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=250,IKB2=10,INF=250,INO=20,IWORK=10)                       
C                                                                               
      REAL FUN(4),COORD(4,2),DERIVT(4,2),DER(2,4),DERIV(2,4),KDERIV(2,4)        
     +,                                                                         
     +JAC(2,2),JAC1(2,2),SAMP(3,2),DTKD(4,4),KP(4,4),PM(4,4),FTF(4,4),          
     +KB(IKB1,IKB2),PB(IKB1,IKB2),LOADS(IKB1),ANS(IKB1),STORPB(INO),            
     +WORK(IWORK,IKB1),COPY(IWORK,IKB1)                                         
      INTEGER NO(INO),G(4),NF(INF,1)                                            
      DATA ICOORD,IDERVT,IDTKD,IKP,IPM,IFTF,NOD/7*4/,ISAMP/3/                   
      DATA IDER,IDERIV,IKDERV,IJAC,IJAC1,IT/6*2/,NODOF/1/                       
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NGP,AA,BB,PERMX,PERMY,UX,UY,                     
     +         DTIM,ISTEP,THETA                                                 
      IWP1=IW+1                                                                 
      IBAND=2*IWP1-1                                                            
      CALL NULL(KB,IKB1,N,IBAND)                                                
      CALL NULL(PB,IKB1,N,IBAND)                                                
      CALL NULL(WORK,IKB2,IWP1,N)                                               
      CALL NULVEC(LOADS,N)                                                      
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
      DO 10 I=1,NN                                                              
   10 NF(I,1)=I                                                                 
C                                                                               
C      ELEMENT INTEGRATION AND ASSEMBLY                                         
C                                                                               
      DO 20 IP=1,NXE                                                            
      DO 20 IQ=1,NYE                                                            
      CALL GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)                        
      CALL NULL(KP,IKP,NOD,NOD)                                                 
      CALL NULL(PM,IPM,NOD,NOD)                                                 
      DO 30 I=1,NGP                                                             
      DO 30 J=1,NGP                                                             
      CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      DO 40 K=1,NOD                                                             
      DO 40 L=1,NOD                                                             
      PART1=PERMX*DERIV(1,K)*DERIV(1,L)+PERMY*DERIV(2,K)*DERIV(2,L)             
      PART2=UX*FUN(K)*DERIV(1,L)+UY*FUN(K)*DERIV(2,L)                           
      DTKD(K,L)=QUOT*(PART1-PART2)                                              
      FTF(K,L)=FUN(K)*FUN(L)*QUOT/(THETA*DTIM)                                  
   40 CONTINUE                                                                  
      CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
      CALL MATADD(PM,IPM,FTF,IFTF,NOD,NOD)                                      
   30 CONTINUE                                                                  
      CALL FORMTB(KB,IKB1,KP,IKP,G,IW,NOD)                                      
      CALL FORMTB(PB,IKB1,PM,IPM,G,IW,NOD)                                      
   20 CONTINUE                                                                  
C                                                                               
C      SPECIFY FIXED NODAL VALUES                                               
C                                                                               
      CALL MATADD(PB,IKB1,KB,IKB1,N,IBAND)                                      
      DO 50 I=1,N                                                               
      DO 50 J=1,IBAND                                                           
   50 KB(I,J)=PB(I,J)-KB(I,J)/THETA                                             
      READ(5,*)IFIX,(NO(I),I=1,IFIX)                                            
      DO 60 I=1,IFIX                                                            
      PB(NO(I),IWP1)=PB(NO(I),IWP1)+1.E20                                       
   60 STORPB(I)=PB(NO(I),IWP1)                                                  
C                                                                               
C      REDUCTION OF LEFT HAND SIDE                                              
C                                                                               
      CALL GAUSBA(PB,IKB1,WORK,IKB2,N,IW)                                       
C                                                                               
C      TIME STEPPING RECURSION                                                  
C                                                                               
      DO 80 J=1,ISTEP                                                           
      CALL MATCOP(WORK,IKB2,COPY,IKB2,IWP1,N)                                   
      CALL BANTML(KB,IKB1,LOADS,ANS,N,IW)                                       
      DO 90 I=1,IFIX                                                            
      IF(J*DTIM.LE..2)ANS(NO(I))=STORPB(I)                                      
      IF(J*DTIM.GT..2)ANS(NO(I))=0.                                             
   90 CONTINUE                                                                  
      CALL SOLVBA(PB,IKB1,COPY,IKB2,ANS,N,IW)                                   
      DO 100 I=1,N                                                              
  100 LOADS(I)=ANS(I)                                                           
      CALL PRINTV(ANS,N)                                                        
   80 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
