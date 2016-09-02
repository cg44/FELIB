      PROGRAM P84                                                               
C                                                                               
C      PROGRAM 8.4  DIFFUSION CONVECTION EQUATION ON A RECTANGLE                
C      USING 4-NODED QUADRILATERAL ELEMENTS. SELF-ADJOINT                       
C      TRANSFORMATION; IMPLICIT INTEGRATION                                     
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=100,IKB2=10,INF=200)                                       
C                                                                               
      REAL FUN(4),COORD(4,2),DERIVT(4,2),DER(2,4),DERIV(2,4),KDERIV(2,4)        
     +,                                                                         
     +JAC(2,2),JAC1(2,2),SAMP(3,2),DTKD(4,4),KP(4,4),PM(4,4),FTF(4,4),          
     +KB(IKB1,IKB2),PB(IKB1,IKB2),LOADS(IKB1),ANS(IKB1)                         
      INTEGER G(4),NF(INF,1)                                                    
      DATA ICOORD,IDERVT,IDTKD,IKP,IPM,IFTF,NOD/7*4/,ISAMP/3/                   
      DATA IDER,IDERIV,IKDERV,IJAC,IJAC1,IT/6*2/,NODOF/1/                       
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,NYE,N,IW,NN,NGP,AA,BB,PERMX,PERMY,UX,UY,                     
     +         DTIM,ISTEP,THETA                                                 
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL NULL(PB,IKB1,N,IWP1)                                                 
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
      PART2=(UX*UX/PERMX+UY*UY/PERMY)*FUN(K)*FUN(L)*.25                         
      DTKD(K,L)=QUOT*(PART1+PART2)                                              
      FTF(K,L)=FUN(K)*FUN(L)*QUOT/(THETA*DTIM)                                  
   40 CONTINUE                                                                  
      CALL MATADD(KP,IKP,DTKD,IDTKD,NOD,NOD)                                    
      CALL MATADD(PM,IPM,FTF,IFTF,NOD,NOD)                                      
   30 CONTINUE                                                                  
C                                                                               
C      INSERT DERIVATIVE BOUNDARY CONDITIONS                                    
C                                                                               
      IF(IQ.NE.1)GOTO 50                                                        
      KP(2,2)=KP(2,2)+UY*AA/6.                                                  
      KP(2,3)=KP(2,3)+UY*AA/12.                                                 
      KP(3,2)=KP(3,2)+UY*AA/12.                                                 
      KP(3,3)=KP(3,3)+UY*AA/6.                                                  
   50 CONTINUE                                                                  
      IF(IQ.NE.NYE)GOTO 60                                                      
      KP(1,1)=KP(1,1)+UY*AA/6.                                                  
      KP(1,4)=KP(1,4)+UY*AA/12.                                                 
      KP(4,1)=KP(4,1)+UY*AA/12.                                                 
      KP(4,4)=KP(4,4)+UY*AA/6.                                                  
   60 CONTINUE                                                                  
      CALL FORMKB(KB,IKB1,KP,IKP,G,IW,NOD)                                      
      CALL FORMKB(PB,IKB1,PM,IPM,G,IW,NOD)                                      
   20 CONTINUE                                                                  
C                                                                               
C      REDUCTION OF LEFT HAND SIDE                                              
C                                                                               
      F1=UY*AA/(2.*THETA)                                                       
      F2=F1                                                                     
      CALL MATADD(PB,IKB1,KB,IKB1,N,IWP1)                                       
      DO 70 I=1,N                                                               
      DO 70 J=1,IWP1                                                            
   70 KB(I,J)=PB(I,J)-KB(I,J)/THETA                                             
      CALL CHOLIN(PB,IKB1,N,IW)                                                 
C                                                                               
C      TIME STEPPING RECURSION                                                  
C                                                                               
      DO 80 J=1,ISTEP                                                           
      CALL BANMUL(KB,IKB1,LOADS,ANS,N,IW)                                       
      ANS(N)=ANS(N)+F1                                                          
      ANS(N-1)=ANS(N-1)+F2                                                      
      CALL CHOBAC(PB,IKB1,ANS,N,IW)                                             
      DO 90 I=1,N                                                               
   90 LOADS(I)=ANS(I)                                                           
      CALL PRINTV(LOADS,N)                                                      
   80 CONTINUE                                                                  
      STOP                                                                      
      END                                                                       
