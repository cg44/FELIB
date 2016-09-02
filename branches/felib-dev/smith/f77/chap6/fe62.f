      PROGRAM P62                                                               
C                                                                               
C      PROGRAM 6.2 PLANE STRAIN OF AN ELASTIC-PLASTIC                           
C      (MOHR-COULOMB) SOLID USING 8-NODE QUADRILATERAL ELEMENTS                 
C      INITIAL STRESS METHOD                                                    
C                                                                               
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKB1=294,IKB2=48,ILOADS=294,INF=176,                    
     +          INX=20,INY=20,INO=10)                                           
C                                                                               
      REAL DEE(4,4),SAMP(4,2),COORD(8,2),JAC(2,2),JAC1(2,2),                    
     +DER(2,8),DERIV(2,8),BEE(4,16),DBEE(4,16),WIDTH(INX),DEPTH(INY),           
     +BTDB(16,16),KM(16,16),ELD(16),EPS(4),SIGMA(4),BLOAD(16),                  
     +BT(16,4),FUN(8),KB(IKB1,IKB2),LOADS(ILOADS),ELOAD(16),                    
     +TOTD(ILOADS),BDYLDS(ILOADS),OLDIS(ILOADS),ELSO(4),                        
     +SX(INX,INY,4),SY(INX,INY,4),TXY(INX,INY,4),SZ(INX,INY,4),                 
     +STORKB(INO),STRESS(4),PL(4,4),GC(2)                                       
      INTEGER NF(INF,2),G(16),NO(INO)                               
      DATA IDEE,IBEE,IDBEE,IH,IFLOW,IPL/6*4/,IDOF,IBTDB,IBT,IKM/4*16/           
      DATA IJAC,IJAC1,NODOF,IT,IDER,IDERIV/6*2/,ICOORD,NOD/2*8/                 
      DATA ISAMP/4/                                                             
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)PHI,C,PSI,GAMA,EPK0,E,V,NXE,NYE,N,IW,NN,NR,NGP                   
      CALL READNF(NF,INF,NN,NODOF,NR)                                  
      READ(5,*)(WIDTH(I),I=1,NXE+1)                                             
      READ(5,*)(DEPTH(I),I=1,NYE+1)                                             
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL NULVEC(OLDIS,N)                                                      
      CALL NULVEC(TOTD,N)                                                       
      CALL NULVEC(BDYLDS,N)                                                     
      CALL FMDRAD(DEE,IDEE,E,V)                                                 
      CALL GAUSS(SAMP,ISAMP,NGP)                                                
C                                                                               
C      ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY                               
C                                                                               
      DO 10 IP=1,NXE                                                            
      DO 10 IQ=1,NYE                                                            
      CALL GEOV8Y(IP,IQ,NYE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)                  
      CALL NULL(KM,IKM,IDOF,IDOF)                                               
      IG=0                                                                      
      DO 20 I=1,NGP                                                             
      DO 20 J=1,NGP                                                             
      IG=IG+1                                                                   
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL GCOORD(FUN,COORD,ICOORD,NOD,IT,GC)                                   
      SY(IP,IQ,IG)=GC(2)*GAMA                                                   
      SX(IP,IQ,IG)=GC(2)*GAMA*EPK0                                              
      SZ(IP,IQ,IG)=GC(2)*GAMA*EPK0                                              
      TXY(IP,IQ,IG)=0.                                                          
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)                      
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)                                    
   20 CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)                                  
   10 CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
C                                                                               
C      READ PRESCRIBED FREEDOMS                                                 
C      AUGMENT AND REDUCE STIFFNESS MATRIX                                      
C                                                                               
      READ(5,*)NL,(NO(I),I=1,NL),PRESC,INCS,ITS                                 
      DO 30 I=1,NL                                                              
      KB(NO(I),IWP1)=KB(NO(I),IWP1)+1.E20                                       
   30 STORKB(I)=KB(NO(I),IWP1)                                                  
      CALL CHOLIN(KB,IKB1,N,IW)                                                 
C                                                                               
C      DISPLACEMENT INCREMENT LOOP                                              
C                                                                               
      DO 40 IY=1,INCS                                                           
      PTOT=PRESC*IY                                                             
      ITERS=0                                                                   
C                                                                               
C      ITERATION LOOP                                                           
C                                                                               
   50 ITERS=ITERS+1                                                             
      CALL NULVEC(LOADS,N)                                                      
      DO 60 I=1,NL                                                              
   60 LOADS(NO(I))=STORKB(I)*PRESC                                              
      CALL VECADD(LOADS,BDYLDS,LOADS,N)                                         
      CALL CHOBAC(KB,IKB1,LOADS,N,IW)                                           
      CALL NULVEC(BDYLDS,N)                                                     
C                                                                               
C      CHECK CONVERGENCE                                                        
C                                                                               
      CALL CHECON(LOADS,OLDIS,N,0.001,ICON)                                     
      IF(ITERS.EQ.1)ICON=0                                                      
C                                                                               
C      INSPECT ALL GAUSS POINTS                                                 
C                                                                               
      DO 70 IP=1,NXE                                                            
      DO 70 IQ=1,NYE                                                            
      CALL NULVEC(BLOAD,IDOF)                                                   
      CALL GEOV8Y(IP,IQ,NYE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)                  
      DO 80 M=1,IDOF                                                            
      IF(G(M).EQ.0)ELD(M)=0.                                                    
   80 IF(G(M).NE.0)ELD(M)=LOADS(G(M))                                           
      IG=0                                                                      
      DO 90 I=1,NGP                                                             
      DO 90 J=1,NGP                                                             
      IG=IG+1                                                                   
      CALL NULVEC(ELSO,IH)                                                      
      CALL FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)                                  
      CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)                     
      CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)                                      
      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)                   
      CALL NULL(BEE,IBEE,IH,IDOF)                                               
      CALL FORMB(BEE,IBEE,DERIV,IDERIV,NOD)                                     
      CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)                                      
      CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)                                     
      CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)                                     
      STRESS(1)=SIGMA(1)+SX(IP,IQ,IG)                                           
      STRESS(2)=SIGMA(2)+SY(IP,IQ,IG)                                           
      STRESS(3)=SIGMA(3)+TXY(IP,IQ,IG)                                          
      STRESS(4)=SIGMA(4)+SZ(IP,IQ,IG)                                           
      CALL INVAR(STRESS,SIGM,DSBAR,THETA)                                       
C                                                                               
C      CHECK WHETHER YIELD IS VIOLATED                                          
C                                                                               
      CALL MOCOUF(PHI,C,SIGM,DSBAR,THETA,FNEW)                                  
      IF(FNEW.LT.0.)GOTO 100                                                    
      STRESS(1)=SX(IP,IQ,IG)                                                    
      STRESS(2)=SY(IP,IQ,IG)                                                    
      STRESS(3)=TXY(IP,IQ,IG)                                                   
      STRESS(4)=SZ(IP,IQ,IG)                                                    
      CALL INVAR(STRESS,SIGM,DSBAR,THETA)                                       
      CALL MOCOUF(PHI,C,SIGM,DSBAR,THETA,F)                                     
      FAC=FNEW/(FNEW-F)                                                         
      STRESS(1)=SX(IP,IQ,IG)+(1.-FAC)*SIGMA(1)                                  
      STRESS(2)=SY(IP,IQ,IG)+(1.-FAC)*SIGMA(2)                                  
      STRESS(3)=TXY(IP,IQ,IG)+(1.-FAC)*SIGMA(3)                                 
      STRESS(4)=SZ(IP,IQ,IG)+(1.-FAC)*SIGMA(4)                                  
      CALL MOCOPL(PHI,PSI,E,V,STRESS,PL)                                        
      CALL MSMULT(PL,IPL,FAC,IH,IH)                                             
      CALL MVMULT(PL,IPL,EPS,IH,IH,ELSO)                                        
      CALL MVMULT(BT,IBT,ELSO,IDOF,IH,ELOAD)                                    
      QUOT=DET*SAMP(I,2)*SAMP(J,2)                                              
      DO 110 K=1,IDOF                                                           
  110 BLOAD(K)=BLOAD(K)+ELOAD(K)*QUOT                                           
  100 IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 90                                     
C                                                                               
C      UPDATE GAUSS POINT STRESSES                                              
C                                                                               
      SX(IP,IQ,IG)=SX(IP,IQ,IG)+SIGMA(1)-ELSO(1)                                
      SY(IP,IQ,IG)=SY(IP,IQ,IG)+SIGMA(2)-ELSO(2)                                
      TXY(IP,IQ,IG)=TXY(IP,IQ,IG)+SIGMA(3)-ELSO(3)                              
      SZ(IP,IQ,IG)=SZ(IP,IQ,IG)+SIGMA(4)-ELSO(4)                                
   90 CONTINUE                                                                  
C                                                                               
C      COMPUTE TOTAL BODYLOADS VECTOR                                           
C                                                                               
      DO 120 M=1,IDOF                                                           
      IF(G(M).EQ.0)GOTO 120                                                     
      BDYLDS(G(M))=BDYLDS(G(M))+BLOAD(M)                                        
  120 CONTINUE                                                                  
   70 CONTINUE                                                                  
      IF(ICON.NE.1.AND.ITERS.NE.ITS)GOTO 50                                     
      CALL VECADD(TOTD,LOADS,TOTD,N)                                            
      PAV=.5*((DEPTH(1)-DEPTH(2))*(SX(1,1,2)+SX(1,1,4))                         
     +       +(DEPTH(2)-DEPTH(3))*(SX(1,2,2)+SX(1,2,4))                         
     +       +(DEPTH(3)-DEPTH(4))*(SX(1,3,2)+SX(1,3,4))                         
     +       +(DEPTH(4)-DEPTH(5))*(SX(1,4,2)+SX(1,4,4)))                        
      WRITE(6,1000)PTOT,PAV,ITERS                                               
      IF(ITERS.EQ.ITS)GOTO 130                                                  
   40 CONTINUE                                                                  
  130 CONTINUE                                                                  
 1000 FORMAT(2E12.4,I12)                                                        
      STOP                                                                      
      END                                                                       
